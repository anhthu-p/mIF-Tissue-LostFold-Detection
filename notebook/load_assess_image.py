import os
import re
from tifffile import imread
import tifffile
import xml.etree.ElementTree as ET
import pandas as pd
import zarr
from itertools import cycle
import numpy as np
from tqdm import tqdm 
import matplotlib.pyplot as plt  # ADDED: needed by get_histogram
import sys
sys.path.append('unify_biomarker_name')
import unify_biomarker_name
from unify_biomarker_name import canonicalize, validate_panel

import importlib
from importlib import reload
reload(unify_biomarker_name)


class FindImage:
    def __init__(self,directory_path,meta_data_path):
        self.dir = directory_path
        self.metadata=pd.read_csv(meta_data_path,sep='\t',index_col=0)

class OmeTifLoader:
    def __init__(self,main_dir, filename):
        self.file_name = filename
        self.filepath = os.path.join(main_dir,filename)
        self.tif = tifffile.TiffFile(self.filepath)
        self.image_store = self.tif.aszarr()
        self.root = ET.fromstring(self.tif.ome_metadata)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.image_store.close()
        self.tif.close()

    def __repr__(self):
        return f"({self.image_store}, {self.root})"

class AssessImage:
    def __init__(self,main_dir, filename):
        self.loader=OmeTifLoader(main_dir,filename)
        self.unify_biomarker_name()
    
    def unify_biomarker_name(self) -> dict: # To see all the biomarkers in this file
        dict_of_biomarker={}
        for element in self.loader.root.iter():
            if element.tag.endswith('Channel'):
                name = canonicalize(element.attrib['Name']) # unify biomarker names to 
                dict_of_biomarker[name]=element.attrib['ID'].split(":")[-1]
        self.biomarker_channel = dict_of_biomarker
        return dict_of_biomarker
    
    def find_biomarkers_of_interest(self,
                                  biomarkers_of_interest:list) -> dict:
        '''
        Find the channels that overlapped with channels of interests in each image
        input: images, list of channels
        Output: Dict{"Channel_name": "Channel index"]}
        '''

        print(f"Processing File:",self.loader.file_name)
        channel_data={}
        for i in biomarkers_of_interest:
            if i in list(self.biomarker_channel):
                channel_index = self.biomarker_channel[i]
                channel_data[i] = int(channel_index)
        missing = [c for c in biomarkers_of_interest if c not in channel_data]
        if missing:
            print(f"Warning: channels not found in { self.loader.file_name}: {missing}")
        self.channel_of_interest_data=channel_data
        return channel_data
    

    def view_images(self,
                    biomarkers_of_interest:list,
                    resolution_level:int,
                    low_high_auto_contrast=None,
                    colors_list=None):
        # put library here because it's heavy instead of putting it on top of the file:
        import napari
        from napari.settings import get_settings
        settings=get_settings()
        settings.application.ipy_interactive=True


        if colors_list is None:
            colors_list = ['cyan', 'magenta', 'yellow', 'red', 'green', 'blue']
        
        viewer=napari.Viewer()
        img_full=zarr.open(self.loader.image_store,mode='r') # put it outside of loop to reduce cost


        for biomarker,channel_index in tqdm(self.biomarker_channel.items(),desc="Finding biomarkers"):
            if biomarker in biomarkers_of_interest:
                img=img_full[str(resolution_level)][int(channel_index)]
                if low_high_auto_contrast:
                    low, high = np.percentile(img, (2.5, 99.8))
                    viewer.add_image(img, name=f'{self.loader.file_name}___{biomarker}', blending='additive',
                                    contrast_limits=[low, high], colormap=colors_list[biomarkers_of_interest.index(biomarker) % len(colors_list)])
                else:
                    viewer.add_image(img, name=f'{self.loader.file_name}___{biomarker}', blending='additive', colormap=colors_list[biomarkers_of_interest.index(biomarker) % len(colors_list)])
        return viewer
    
    def assess_single_image(self, biomarker: str, resolution_level: int):
        image_full = zarr.open(self.loader.image_store, mode='r')
        single_image = image_full[str(resolution_level)][int(self.biomarker_channel[biomarker])]
        
        return single_image


    ################################################____________STASTS_____________#########################################################################
    '''
    Metric-------------------------	Starting threshold--------------------Reasoning
    5th percentile (background)	    > 5% of max bit depth	              For 16-bit images (max 65535), flag if 5th pct > ~1000
    SNR	                            < 3 → poor, > 10 → good	              Physics-based: 3-sigma rule is the minimum to distinguish signal from noise
    Dynamic range	                < 10% of bit depth → flat	          For 16-bit, flag if 99th - 5th pct < ~2000
    CV of background	            < 0.1 → uniform, > 0.3 → uneven	      Rough split; calibrate on your data
    '''

    def get_stats(self, resolution_level) -> pd.DataFrame:
        image_full = zarr.open(self.loader.image_store, mode='r')
        rows = []

        for biomarker, channel_index in self.biomarker_channel.items():
            image = image_full[str(resolution_level)][int(channel_index)]
            bg = image[image < np.percentile(image, 25)]

            fifth_percentile = np.percentile(image, 5)
            cv_background    = np.std(bg) / np.mean(bg)
            snr              = (np.percentile(image, 99) - np.percentile(image, 5)) / np.std(bg)
            dynamic_range    = np.percentile(image, 99) - np.percentile(image, 5)

            rows.append({
                'biomarker'    : biomarker,
                'channel_index': int(channel_index), 
                'fifth_pct'    : fifth_percentile,
                'cv_background': cv_background,
                'snr'          : snr,
                'dynamic_range': dynamic_range,
                'flag'         : self._generate_flag(fifth_percentile, cv_background, snr, dynamic_range)
            })

        return pd.DataFrame(rows)
    
    def _generate_flag(self, fifth_pct, cv_bg, snr, dyn_range) -> str:
        high_bg = fifth_pct > 500
        high_cv = cv_bg > 0.3
        low_cv  = cv_bg < 0.1
        low_snr = snr < 1
        flat    = dyn_range < 200

        if high_bg and high_cv:
            return "CRITICAL: uneven illumination → BaSiCPy"
        elif high_bg and low_snr:
            return "CRITICAL: background masking signal → background subtraction first"
        elif high_bg and low_cv:
            return "WARNING: uniform background → rolling ball subtraction"
        elif low_snr and flat:
            return "WARNING: dim signal → CLAHE"
        elif low_snr:
            return "WARNING: poor SNR → check staining quality"
        elif flat:
            return "WARNING: narrow dynamic range → CLAHE or normalization"
        else:
            return "OK"

    def get_histogram(self, biomarker, resolution_level):
        image_full = zarr.open(self.loader.image_store, mode='r')  # FIXED: was self.assessor.loader (old train wreck)
        if biomarker in self.biomarker_channel:  # FIXED: was self.assessor.biomarker_channel().keys — dict, not a callable; use `in dict` directly
              channel_index = self.biomarker_channel[biomarker]  # FIXED: was ['biomarker'] string literal instead of variable
              image = image_full[str(resolution_level)][int(channel_index)]
              plt.hist(image.ravel(), bins=256)


    #### SET UP FOR INSTANSEG#####
    def max_projection(self, markers: list, resolution_level: int) -> np.ndarray:
        image_full = zarr.open(self.loader.image_store, mode='r')
        stack = []

        for marker in markers:
            if marker not in self.biomarker_channel:
                print(f"Warning: {marker} not found, skipping")
                continue
            image = np.array(image_full[str(resolution_level)][int(self.biomarker_channel[marker])])

            p_low, p_high = np.percentile(image, [0.5, 99.5])
            normalized = (image - p_low) / (p_high - p_low)  # scale each channel to its own [0,1] range
            stack.append(np.clip(normalized, 0, 1))

        return np.max(np.stack(stack), axis=0)

    



       
