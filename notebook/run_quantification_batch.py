"""
Batch cell quantification (shape + per-biomarker intensity stats) across all
segmentation .npz files in mask_dir, parallelized across CPU cores.

Moved out of 5_Quantification.ipynb because Windows' spawn-based
multiprocessing needs the worker function importable from a real module --
a function defined in a notebook cell can't be pickled to worker processes.

Usage:
    python run_quantification_batch.py
"""
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table
import load_assess_image as load_assess


def single_sample_quantification(main_dir, mask_dir, save_dir, file_name):
    stem = file_name.replace('.npz', '')
    npz_path = os.path.join(mask_dir, file_name)
    data = np.load(npz_path)
    nuclei_mask = data['nuclei']

    image_test = load_assess.AssessImage(main_dir, stem + '.ome.tif')
    shape_props = regionprops_table(
        nuclei_mask,
        properties=['label', 'area', 'centroid', 'equivalent_diameter_area',
                    'solidity', 'eccentricity', 'axis_major_length', 'axis_minor_length', 'perimeter']
    )
    df = pd.DataFrame(shape_props)
    df['sample_id'] = stem

    for biomarker, channel_index in image_test.biomarker_channel.items():
        image_biomarker = np.array(image_test.assess_single_image(biomarker=biomarker, resolution_level=0))
        intensity_props = regionprops_table(
            nuclei_mask,
            intensity_image=image_biomarker,
            properties=['label', 'intensity_mean', 'intensity_median', 'intensity_std']
        )
        intensity_df = pd.DataFrame(intensity_props).rename(columns={
            'intensity_mean':   f'{biomarker}_mean',
            'intensity_median': f'{biomarker}_median',
            'intensity_std':    f'{biomarker}_std',
        })
        df = df.merge(intensity_df, on='label')

    output_path = os.path.join(save_dir, stem + '.csv')
    df.to_csv(output_path, index=True)
    return f"Saved {len(df)} cells -> {output_path}"


if __name__ == "__main__":
    main_dir = r'D:\thu\HTAN\images\ILC\level_2'
    mask_dir = r'D:\thu\HTAN\images\mask\test\ILC\qc'
    save_dir = r'D:\thu\HTAN\images\quantification\ILC\qc'
    os.makedirs(save_dir, exist_ok=True)

    npz_files = [f for f in os.listdir(mask_dir) if f.endswith('.npz')]

    max_workers = max(1, os.cpu_count() - 2)
    print(f"Processing {len(npz_files)} files with {max_workers} workers")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(single_sample_quantification, main_dir, mask_dir, save_dir, f): f
            for f in npz_files
        }
        for future in as_completed(futures):
            f = futures[future]
            try:
                print(future.result())
            except Exception as e:
                print(f"{f} failed: {e}")

    print("Done.")
