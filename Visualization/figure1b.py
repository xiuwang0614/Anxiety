# -*- coding: utf-8 -*-
"""
Created on Tue May 14 18:07:05 2024

@author: Lenovo
"""

# First we fetch the Yeo atlas
from nilearn import datasets
import nibabel as nib
import os

# Specify the directory path you want to set as the working directory
directory_path = 'I:\\lda\\0419\\0mask'

# Set the working directory
os.chdir(directory_path)

#mne.datasets.sample.data_path('C:\Users\Eric\data')
#atlas_yeo_2011 = datasets.fetch_atlas_yeo_2011()
#atlas_yeo = atlas_yeo_2011.thick_7

sub = nib.load('figure1b.nii')


# Let's now plot it
from nilearn import plotting
p=plotting.plot_roi(
    sub,
    #title="Original Yeo atlas",
    cut_coords=(-17, 0, -3),
    colorbar=False,
    cmap="Paired",
    draw_cross=False,
    #colorbar=False
)


fig = plotting.show()
save_path = "I:\\lda\\0419\\brain\\figure1b.tiff"
p.savefig(save_path, dpi=600)
