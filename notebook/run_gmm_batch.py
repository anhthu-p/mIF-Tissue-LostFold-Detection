"""
Batch GMM clustering of a biomarker across all quantification .csv files in
quant_dir, parallelized with a thread pool.

Moved out of 7_analysis.ipynb so the worker function is a real, importable
module-level function.

A thread pool (not a process pool) is used here because each file's work is
mostly numpy/scikit-learn linear algebra (GMM.fit releases the GIL during
its BLAS calls) plus CSV I/O -- both benefit from threads without the
per-worker process/import overhead ProcessPoolExecutor would add for a
comparatively lightweight per-file job.

Usage:
    python run_gmm_batch.py
"""
import os
import sys


import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.append(r"R:\Thu\From_2026\1. Projects\TP121_BRC962_Reorganized")
import gmm_func

import re

def get_dna_round_order(columns):
    """
    Ordered list of DNA round names as they appear in df's columns
    (e.g. ['dna', 'dna(2)', 'dna(3)', ...]), matching the imaging-cycle order.
    """
    pattern = re.compile(r'^dna(\(\d+\))?_mean$')
    return [col[:-len('_mean')] for col in columns if pattern.match(col)]


def exclude_flagged_cells(df, flagged_bounds_df, x_col='centroid-1', y_col='centroid-0'):
    """
    For each round_dna present in flagged_bounds_df, find cells whose centroid
    falls inside any flagged tile for that round, and set NaN across the
    imaging-cycle column block from that round's _mean through the NEXT
    round's _std (inclusive) -- i.e. the DNA channel + markers stained in
    that same physical cycle, since a fold/lost-tissue artifact detected at
    a given DNA round corrupts everything imaged in that cycle.

    Assumes centroid-0 = row = y, centroid-1 = col = x (skimage regionprops
    convention), matching how x0/y0/x1/y1 were built by
    QC_Fold_Lost.get_flagged_tile_bounds().
    """
    df = df.copy()
    dna_round_order = get_dna_round_order(df.columns)

    for round_dna, group in flagged_bounds_df.groupby('round_dna'):
        if round_dna not in dna_round_order:
            print(f"Skipping {round_dna}: not found in df columns")
            continue

        start_col = f'{round_dna}_mean'
        idx = dna_round_order.index(round_dna)
        end_col = f'{dna_round_order[idx + 1]}_std' if idx + 1 < len(dna_round_order) else df.columns[-1]

        x = df[x_col].to_numpy()
        y = df[y_col].to_numpy()
        in_any_flagged_tile = np.zeros(len(df), dtype=bool)
        for _, tile in group.iterrows():
            in_any_flagged_tile |= (
                (x >= tile['x0']) & (x < tile['x1']) &
                (y >= tile['y0']) & (y < tile['y1'])
            )

        df.loc[in_any_flagged_tile, start_col:end_col] = np.nan
        print(f"{round_dna}: {in_any_flagged_tile.sum()} cells NaN'd across {start_col}..{end_col}")

    return df

def single_file_gmm(file_name, quant_dir, save_dir, biomarker='pRB_median', cluster_range=range(1, 10), log_or_raw = 'log', show_plot=True):
    stem = file_name.replace('.csv', '')
    df = pd.read_csv(os.path.join(quant_dir, file_name), index_col=0)
    os.makedirs(save_dir, exist_ok=True)
    qc_tissues = pd.read_csv(rf"D:\thu\HTAN\images\QC_Fold_Lost\{subtype}\{file_name.split('.')[0]}\{file_name.split('.')[0]}_flagged_bounds.csv")
    df_qc = exclude_flagged_cells(df, qc_tissues, x_col='centroid-1', y_col='centroid-0')
    df_qc[f'log_{biomarker}'] = np.log1p(df_qc[biomarker])
    if log_or_raw == 'raw':     
        obj = gmm_func.FindClusterGMM(biomarker=f'{biomarker}', cluster_range=cluster_range, df_or_path=df_qc)
    else:
        obj = gmm_func.FindClusterGMM(biomarker=f'log_{biomarker}', cluster_range=cluster_range, df_or_path=df_qc)
    obj.calculate_bic()
    obj.find_optimal_n()
    obj.predict_label()
    labeled_df = obj.assign_labels_to_df()

    if show_plot:
        fig = obj.plot_histogram_and_bic(show=False)
        out_path = os.path.join(save_dir, f"{file_name.split('.')[0]}_{biomarker}_gmm_hist_bic.png")
        fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"Saved plot to {out_path}")

    output_path = os.path.join(save_dir, stem + '.csv')
    if log_or_raw == 'raw':
        labeled_df[[f'cluster_{biomarker}','label']].to_csv(output_path, index=True)
    else:
        labeled_df[[f'cluster_log_{biomarker}','label']].to_csv(output_path, index=True)
    return f"{stem}: optimal_n={obj.optimal_n} -> Saved {len(labeled_df)} cells -> {output_path}"


if __name__ == "__main__":
    
    # Load er_pos samples in each subtype
    er_pos_combined_pRB_ki67 = pd.read_csv(rf"R:\Thu\From_2026\6. Other\DL_Summer_2026\result\er_pos_combined_pRB_ki67_metadata.csv",index_col=0)
    er_pos_ILC=er_pos_combined_pRB_ki67[er_pos_combined_pRB_ki67['Diagnosis']=='Lobular carcinoma NOS']['Filename_simple'].tolist()
    er_pos_NST = er_pos_combined_pRB_ki67[er_pos_combined_pRB_ki67['Diagnosis']=='Ductal carcinoma NOS']['Filename_simple'].tolist()

    subtype = 'ILC'
    quant_dir = rf"D:\thu\HTAN\images\quantification\{subtype}\qc"

    #datatype to use for GMM: 'log' or 'raw'
    log_or_raw = 'log'

    biomarker_list = ['aSMA']

    
    # Find all CSV files in quant_dir that have a corresponding .ome.tif file in er_pos_ILC
    csv_files = []
    
    for biomarker in biomarker_list:
        # if biomarker is not in the quantification CSV, skip it
        quant_dir = rf"D:\thu\HTAN\images\quantification\{subtype}\qc"
        sample_csv = next((f for f in os.listdir(quant_dir) if f.endswith('.csv')), None)
        if subtype == 'ILC':    
            main_dir = rf"D:\thu\HTAN\images\{subtype}\level_2"
            for i in range(len(er_pos_ILC)):
                file_name = er_pos_ILC[i]
                image_test = AssessImage(main_dir, file_name)
                if biomarker not in image_test.biomarker_channel:
                    print(f"Biomarker {biomarker} not found in {file_name}. Skipping.")
                    continue
            for f in os.listdir(quant_dir):
                if f.endswith('.csv') and f.split('.')[0] + '.ome.tif' in er_pos_ILC:
                    csv_files.append(f)



        elif subtype == 'NST':
            main_dir = rf"D:\thu\HTAN\images\{subtype}"
            for i in range(len(er_pos_NST)):
                file_name = er_pos_NST[i]
                image_test = AssessImage(main_dir, file_name)
                if biomarker not in image_test.biomarker_channel:
                    print(f"Biomarker {biomarker} not found in {file_name}. Skipping.")
                    continue
            for f in os.listdir(quant_dir):
                if f.endswith('.csv') and f.split('.')[0] + '.ome.tif' in er_pos_NST:
                    csv_files.append(f)

        if sample_csv is None:
            print(f"No CSV files found in {quant_dir}. Skipping biomarker {biomarker}.")
            continue

    

        # Create save directory for the biomarker
        save_dir = rf"D:\thu\HTAN\images\analysis\{subtype}\gmm\test\{log_or_raw}_{biomarker}"
        os.makedirs(save_dir, exist_ok=True)
        
        # calculate GMM for the biomarker across all CSV files in quant_dir
        # set cluster_range to the range of cluster numbers to test for BIC
        cluster_range = range(1, 10)


       
        # parallelize the processing of CSV files using ThreadPoolExecutor        
        max_workers = min(26, len(csv_files)) or 1
        print(f"Processing {len(csv_files)} files with {max_workers} threads")

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(single_file_gmm, f, quant_dir, save_dir, biomarker, cluster_range, log_or_raw): f
                for f in csv_files
            }

            for future in as_completed(futures):
                f = futures[future]
                try:
                    print(future.result())
                except Exception as e:
                    print(f"{f} failed: {e}")

        print("Done.")
