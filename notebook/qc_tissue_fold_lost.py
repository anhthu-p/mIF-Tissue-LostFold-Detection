from pathlib import Path
import numpy as np
from skimage import filters
import load_assess_image
import napari
from importlib import reload
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import load_assess_image


def norm(x):
    x = x.astype('float32')
    lo, hi = np.percentile(x, [1, 99])
    return np.clip((x - lo) / (hi - lo + 1e-9), 0, 1)

class QC_Fold_Lost:
    def __init__(self, main_dir: str, filename: str):
        self.main_dir = main_dir
        self.filename = filename
        self.full_image = load_assess_image.AssessImage(main_dir=main_dir, filename=filename)
    

    # ---------- 0. shared image loading + whole-image agreement ----------
    GLOBAL_DISAGREE_THRESHOLD = 0.7  # placeholder -- calibrate on known good/bad pairs

    def _load_round(self, round_dna: str):
        """
        Load (and cache) img_a ('dna', round 0) + its Otsu threshold once per
        instance, and load img_b for the requested round_dna. Shared by
        global_disagreement() and tile_disagreement() so either can be called
        on its own, in any order, without redoing work.
        """
        if not hasattr(self, 'image_a'):
            self.image_a = self.full_image.assess_single_image('dna', 0)
            self._thr_a = filters.threshold_otsu(self.image_a)

        if getattr(self, 'round_dna', None) != round_dna:
            self.image_b = self.full_image.assess_single_image(round_dna, 0)
            self.round_dna = round_dna

    def global_disagreement(self, round_dna: str = None) -> float:
        """
        1 - NCC computed over ALL foreground pixels at once (not per tile).
        Restricting to foreground (image_a > thr) avoids background pixels
        (agreeing trivially near zero) diluting the score.

        round_dna : which round to compare against 'dna' (round 0). Can be
                    called standalone (loads/caches images itself) or omitted
                    to reuse whatever round is already loaded (e.g. right
                    after tile_disagreement()).
        """
        if round_dna is not None:
            self._load_round(round_dna)

        mask = self.image_a > self._thr_a
        a = self.image_a[mask].astype('float64')
        b = self.image_b[mask].astype('float64')
        if a.size == 0:
            return np.nan

        a = a - a.mean()
        b = b - b.mean()
        denom = np.sqrt((a*a).sum() * (b*b).sum()) + 1e-9
        ncc = float((a*b).sum() / denom)
        self.global_disagree = 1.0 - ncc
        return 1.0 - ncc


    # ---------- 1. per-tile disagreement between two registered rounds ----------
    def tile_disagreement(self, round_dna: str, # dna(1), dna(2), etc.
                           tile: int, min_foreground: float = 0.05, verbose: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        img_a, img_b : 2D float arrays, same shape, already co-registered.
        tile         : tile size in pixels. # eg 128 pixels
        min_foreground: skip near-empty tiles (fraction of pixels above threshold).
        verbose      : print per-round progress messages.

        Returns three (rows, cols) grids:
        disagree  : 1 - NCC per tile (NaN where tile skipped for low content)
        dmean     : mean(img_b) - mean(img_a) per tile  (sign of intensity change)
        fg        : foreground fraction per tile

        If whole-image disagreement is already very high (gross misregistration /
        little-to-no overlap), the per-tile median+k*MAD logic in flag_folds() is
        unreliable (its "typical" baseline gets contaminated by the majority-bad
        tiles). In that case we skip per-tile scoring entirely and flag the whole
        round as failed registration via self.registration_failed.
        """
        self._load_round(round_dna)
        img_a = self.image_a
        img_b = self.image_b
        thr = self._thr_a
        self.tile = tile

        if verbose:
            print(f"Round {round_dna}")

        H, W = img_a.shape
        rows, cols = H // tile, W // tile

        self.global_disagree = self.global_disagreement()
        self.registration_failed = self.global_disagree > self.GLOBAL_DISAGREE_THRESHOLD

        if self.registration_failed:
            if verbose:
                print(f"Round {round_dna}: global_disagree={self.global_disagree:.3f} "
                      f"> {self.GLOBAL_DISAGREE_THRESHOLD} -> flagging as failed registration, skipping per-tile scoring")
            disagree = np.full((rows, cols), np.nan)
            dmean    = np.full((rows, cols), np.nan)
            fg       = np.zeros((rows, cols))
            self.tile_disagreement_array = disagree
            self.tile_dmean_array = dmean
            self.tile_foreground_array = fg
            self.fold_tile = []
            self.lose_tile = []
            return disagree, dmean, fg

        # Compute disagreement and mean difference per tile
        disagree = np.full((rows, cols), np.nan)
        dmean    = np.full((rows, cols), np.nan)
        fg       = np.zeros((rows, cols))
        real_signal = np.zeros((rows, cols))

        if verbose:
            print(f"Processing {rows} x {cols} tiles of size {tile}x{tile}")


        negative_intensity_tiles = []
        for r in range(rows):
            for c in range(cols):
                ta = img_a[r*tile:(r+1)*tile, c*tile:(c+1)*tile].ravel()
                tb = img_b[r*tile:(r+1)*tile, c*tile:(c+1)*tile].ravel()

                frac = np.mean(ta > thr) # find mean foreground in a specific tile that pass otsu threshold 
                fg[r, c] = frac # save fg for this specific tile
                
                if frac < min_foreground: 
                    continue  # too little tissue to judge -> leave NaN

                # NCC: normalized cross-correlation (brightness/contrast invariant)
                a = ta - ta.mean()
                b = tb - tb.mean()
                denom = np.sqrt((a*a).sum() * (b*b).sum()) + 1e-9
                ncc = float((a*b).sum() / denom)


                disagree[r, c] = 1.0 - ncc
                dmean[r, c]    = tb.mean() - ta.mean()
                # print(f"Tile {r},{c}: NCC={ncc:.3f}, disagree={disagree[r,c]:.3f}, dmean={dmean[r,c]:.3f}, fg={frac:.3f}")
                self.tile_disagreement_array = disagree
                self.tile_dmean_array = dmean
                self.tile_foreground_array = fg

        return disagree, dmean, fg


    # ---------- 2. flag fold-candidate tiles ----------
    def flag_folds(self, k: float = 3.0) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
        """
        Flag tiles that are BOTH:
        - spatial outliers in disagreement (robust: median + k*MAD), AND
        - positive intensity change (dmean > 0, i.e. doubled DNA).

        Returns boolean grid of fold-candidate tiles.
        """
        vals = self.tile_disagreement_array[~np.isnan(self.tile_disagreement_array)]
        if vals.size == 0:
            self.fold_tile = []
            self.lose_tile = []
            return self.fold_tile, self.lose_tile

        med = np.median(vals)
        mad = np.median(np.abs(vals - med))
        # 1.4826 scales MAD to be comparable to a standard deviation for
        # normally distributed data -- the adjustment I mentioned earlier.
        robust_sigma = 1.4826 * mad + 1e-9
        threshold = med + k * robust_sigma

        is_spatial_outlier = self.tile_disagreement_array > threshold
        # print(f"Number of spatial outliers: {np.sum(is_spatial_outlier)}")
        is_losing_tissue = self.tile_dmean_array < 0
        is_fold_candidate = self.tile_dmean_array > 0

        array_lose = np.nan_to_num(is_spatial_outlier & is_losing_tissue).astype(bool)
        array_fold = np.nan_to_num(is_spatial_outlier & is_fold_candidate).astype(bool)

        fold_tile = [tuple(rc) for rc in np.argwhere(array_fold)]
        lose_tile = [tuple(rc) for rc in np.argwhere(array_lose)]

        self.fold_tile = fold_tile
        self.lose_tile = lose_tile

        return fold_tile, lose_tile


    # ---------- 2b. pixel bounding boxes of flagged tiles, for downstream cell exclusion ----------
    def get_flagged_tile_bounds(self) -> list[dict]:
        """
        Must be called after tile_disagreement + flag_folds for the current round.
        Returns one record per flagged tile with its pixel bounding box, so a
        later step can test cell centroids (x, y) against (x0 <= x < x1 and
        y0 <= y < y1) to exclude cells sitting in fold/lost-tissue tiles.
        (x0,y0)---(x1,y0)
        |          |
        (x0,y1)---(x1,y1)

        Row/col -> pixel mapping matches plot_flagged_tiles(): row -> y, col -> x.
        """
        records = []
        for tile_type, tiles in (('fold', self.fold_tile), ('lose', self.lose_tile)):
            for (r, c) in tiles:
                x0, y0 = c * self.tile, r * self.tile
                x1, y1 = x0 + self.tile, y0 + self.tile
                records.append({
                    'filename': self.filename,
                    'round_dna': self.round_dna,
                    'tile_type': tile_type,
                    'row': r,
                    'col': c,
                    'x0': x0, 'y0': y0,
                    'x1': x1, 'y1': y1,
                })
        return records


    # ---------- 3. per-round QC metric for CSV export ----------
    WARN_RATIO_THRESHOLD = 0.10

    def get_qc_record(self) -> dict:
        """
        Summarize the current round (must be called after tile_disagreement + flag_folds)
        as a dict suitable for appending to a QC summary table.
        """
        n_flagged = len(self.fold_tile) + len(self.lose_tile)
        # tiles below min_foreground were left as NaN in tile_disagreement -> exclude them
        # so the ratio is computed over tissue-containing tiles only, not background
        total_tiles = int(np.sum(~np.isnan(self.tile_disagreement_array)))
        ratio = n_flagged / total_tiles if total_tiles else 0.0

        # Consider gloabl disagreement and flagged ratio to determine warning level:

        if self.registration_failed:
            warning = 'Majority is lost/fold - Manual review'
        elif ratio > self.WARN_RATIO_THRESHOLD:
            warning = 'Flag manually check'
        else:
            warning = 'PASS'

        return {
            'filename': self.filename,
            'round_dna': self.round_dna,
            'global_disagree': self.global_disagree,
            'registration_failed': self.registration_failed,
            'n_flagged_tiles': n_flagged,
            'total_tiles': total_tiles,
            'flagged_ratio': ratio,
            'warning': warning,
        }
    

    #################################################
    # PLOT
    def plot_flagged_tiles(self, figsize=(10, 5), ax=None, show=False, save_path=None):
        """
        image_a, image_b : 2D float arrays, same shape, already co-registered.
        lose_tile     : list of (row, col) tuples for losing tissue tiles
        fold_tile     : list of (row, col) tuples for fold candidate tiles
        tile          : tile size in pixels
        ax            : optional pair of matplotlib Axes (ax_merge, ax_flagged) to draw into;
                        if None, a new 1x2 figure is created
        save_path     : optional path (file or folder) to save the figure to.
                        If a folder is given, the filename is derived from self.filename.
        """

        external_ax = ax is not None
        if not external_ax:
            fig, (ax_merge, ax_flagged) = plt.subplots(1, 2, figsize=figsize)
        else:
            ax_merge, ax_flagged = ax

        a, b = norm(self.image_a), norm(self.image_b)   # both normalized to [0,1]
        overlay = np.zeros((*a.shape, 3), dtype='float32')
        overlay[..., 0] = a      # red   channel = round A
        overlay[..., 2] = a      # blue  channel = round A   -> A shows as magenta
        overlay[..., 1] = b

        ax_merge.imshow(overlay)
        ax_merge.set_title(f"{self.round_dna}: merged channels")
        ax_merge.axis('off')

        ax_flagged.imshow(overlay)

        for (r, c) in self.lose_tile:
            rect = patches.Rectangle(
                (c * self.tile, r * self.tile),    # (x, y) top-left corner — note x=col, y=row
                self.tile, self.tile,              # width, height
                linewidth=1.5, edgecolor='red', facecolor='none')
            ax_flagged.add_patch(rect)
            # ax_flagged.text(c * self.tile + 4, r * self.tile + 16, f"{r},{c}",
            #         color='red', fontsize=8)

        for (r, c) in self.fold_tile:
            rect = patches.Rectangle(
                (c * self.tile, r * self.tile),    # (x, y) top-left corner — note x=col, y=row
                self.tile, self.tile,              # width, height
                linewidth=1.5, edgecolor='white', facecolor='none')
            ax_flagged.add_patch(rect)
            # ax_flagged.text(c * self.tile + 4, r * self.tile + 16, f"{r},{c}",
            #         color='white', fontsize=8)

        ax_flagged.set_title(f"{self.round_dna}: {len(self.lose_tile)+ len(self.fold_tile)} flagged tiles")
        ax_flagged.axis('off')

        legend_handles = [
            patches.Patch(edgecolor='red', facecolor='none', label='lost tissue'),
            patches.Patch(edgecolor='white', facecolor='none', label='fold tissue'),
        ]
        ax_flagged.legend(handles=legend_handles, loc='lower right', fontsize=8, framealpha=0.7)

        if save_path is not None:
            save_path = Path(save_path)
            if save_path.suffix == '':
                save_path.mkdir(parents=True, exist_ok=True)
                save_path = save_path / f"{Path(self.filename).stem}_flagged_tiles_{self.round_dna}.png"
            else:
                save_path.parent.mkdir(parents=True, exist_ok=True)
            ax_flagged.figure.savefig(save_path, dpi=150, bbox_inches='tight')

        if show:
            plt.tight_layout()
            plt.show()
        elif not external_ax:
            plt.close(fig)



