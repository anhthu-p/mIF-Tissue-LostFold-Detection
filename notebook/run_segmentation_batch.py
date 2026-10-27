import numpy as np
import pandas as pd
import skimage as ski
import matplotlib.pyplot as plt
import tifffile
import zarr
import os
import xml.etree.ElementTree as ET
import napari
from importlib import reload
from instanseg import InstanSeg
import load_assess_image as load_assess
import qc_widget
from concurrent.futures import ThreadPoolExecutor   # was ProcessPoolExecutor
from skimage.measure import regionprops_table
import torch
print(torch.cuda.is_available())
print(torch.cuda.get_device_name(0))  # shows your GPU name



# --- parallel MIP (CPU/IO bound) ---
def process_MIP(main_dir, file_name, cytomarker:list=['CD45','aSMA','panCK','Vim','PDPN','CK7']) -> tuple[str, np.ndarray, np.ndarray, float]:
    image = load_assess.AssessImage(main_dir, file_name)
    image_full = zarr.open(image.loader.image_store, mode='r')
    dna_idx = 0
    dna = np.array(image_full['0'][dna_idx])
    cyto = image.max_projection(cytomarker, resolution_level=0)
    pixel_size = 0.0
    for element in image.loader.root.iter():
        if element.tag.endswith('Pixels'):
            pixel_size = float(element.attrib['PhysicalSizeX'])

    image.loader.tif.close()
    image.loader.image_store.close()
    return file_name, dna, cyto, pixel_size

def get_qc_warning(df, round_dna):
    match = df.loc[df['round_dna'] == f'dna({round_dna})', 'warning']
    if match.empty:
        return None
    return match.iloc[0]

def save_mask(filename, labeled_output, save_dir):
    stem = filename.replace('.ome.tif', '')
    output_path = os.path.join(save_dir, stem + '.npz')
    nuclei_mask = labeled_output[0][0, 0].cpu().numpy().astype(int)
    cells_mask  = labeled_output[0][0, 1].cpu().numpy().astype(int)
    np.savez_compressed(output_path, nuclei=nuclei_mask, cells=cells_mask)
    return stem


if __name__ == "__main__":
    root_dir=rf'D:\thu\HTAN\images'
    subtype = 'NST'
    main_dir = os.path.join(root_dir, subtype)
    save_dir = os.path.join(root_dir, 'mask', 'test', subtype, 'qc')

    # Find files with both pRB and Ki67 biomarkers
    list_file = []
    for i in os.listdir(main_dir):
        if i.endswith('.ome.tif'):
            image_test = load_assess.AssessImage(main_dir, i)
            if 'pRB' in image_test.biomarker_channel:
                if 'Ki67' in image_test.biomarker_channel:
                    list_file.append(i)

    #Filter files that pass QC of tissue lost and fold 
    qc_file = []
    for file in list_file:
        image_test = load_assess.AssessImage(main_dir, file)
        if 'pRB' in image_test.biomarker_channel and 'Ki67' in image_test.biomarker_channel:
            pRB_DNA_round = int(int(image_test.biomarker_channel['pRB']) / 4)
            Ki67_DNA_round = int(int(image_test.biomarker_channel['Ki67']) / 4)

            df = pd.read_csv(rf"{root_dir}\QC_Fold_Lost\{subtype}\{file.split('.')[0]}\{file.split('.')[0]}_qc.csv")

            if get_qc_warning(df, pRB_DNA_round) == 'PASS':
                if get_qc_warning(df, Ki67_DNA_round) == 'PASS':
                    qc_file.append(file)

    # --- parallel MIP (CPU/IO bound) ---
    with ThreadPoolExecutor(max_workers=20) as pool:
        mip_results = list(pool.map(lambda i: process_MIP(main_dir, i), qc_file))

    # --- sequential segmentation (GPU bound) ---
    model = InstanSeg('fluorescence_nuclei_and_cells', device='cuda')

    with ThreadPoolExecutor(max_workers=20) as save_pool:   # separate pool from the MIP one
        futures = []
        for file_name, dna, cyto, pixel_size in mip_results:
            stack = np.stack([dna, cyto])
            labeled_output = model.eval_medium_image(stack, pixel_size=pixel_size)

            future = save_pool.submit(save_mask, file_name, labeled_output, save_dir)
            futures.append(future)

        # block here until all saves are done, and surface any exceptions
        for f in futures:
            f.result()

    print("DONE: Segmentation and saving of masks completed for all files.")
