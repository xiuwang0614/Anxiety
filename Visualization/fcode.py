# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 13:38:02 2024

@author: Lenovo
"""

from neurosynth.base.dataset import Dataset
from neurosynth import decode


import os

import nibabel as nib
import numpy as np
from nilearn.plotting import plot_roi

from nimare.dataset import Dataset
from nimare.decode import discrete
from nimare.utils import get_resource_path
import nimare 

from pprint import pprint

from nimare.extract import download_abstracts, fetch_neuroquery, fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset

# biopython is unnecessary here, but is required by download_abstracts.
# We import it here only to document the dependency and cause an early failure if it's missing.
import Bio  # pip install biopython

#https://nimare.readthedocs.io/en/stable/decoding.html
nimare.extract.fetch_neurosynth(path="C:\\Users\\Lenovo\\.nimare\\neurosynth")

out_dir = os.path.abspath("E:\\2DCM\\Neurosynth\\resultLDA")
os.makedirs(out_dir, exist_ok=True)

files = fetch_neurosynth(
    data_dir=out_dir,
    version="7",
    overwrite=False,
    source="abstract",
    vocab="LDA50",
)

# Note that the files are saved to a new folder within "out_dir" named "neurosynth".
pprint(files)
neurosynth_db = files[0]

neurosynth_dset = convert_neurosynth_to_dataset(
    coordinates_file=neurosynth_db["coordinates"],
    metadata_file=neurosynth_db["metadata"],
    annotations_files=neurosynth_db["features"],
)
dset2 = neurosynth_dset
from nilearn import datasets, image, plotting
#from nimare.utils import get_template
#from nilearn.image import resample_to_img
mask_img= image.load_img("I:\\lda\\0803\\neurosynth\\your_file_binary.nii")
from nimare.decode.discrete import ROIAssociationDecoder

decoder = ROIAssociationDecoder(
   mask_img,
   u=0.05,
   correction="fdr_bh",
)
decoder.fit(dset2)
decoding_results = decoder.transform()
decoding_results.to_csv(r'I:\\lda\\0803\\neurosynth\0831\\r_caudate_400_th0.01.csv', sep=',', encoding='utf-8', header='true')
















