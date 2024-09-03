# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 00:44:05 2024

@author: Lenovo
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# Install and load required libraries 
# e.g. pip install numpy brainspace surfplot neuromaps nibabel 

from surfplot import Plot
from brainspace.datasets import load_parcellation
from brainspace.mesh.mesh_io import read_surface
from neuromaps.datasets import fetch_fslr
import numpy as np
import pandas as pd
import os

# Specify the directory path you want to set as the working directory
directory_path = 'I:\\lda\\0803\\CONJ'

# Set the working directory
os.chdir(directory_path)
#from numpy import unique
# Read the CSV file into a DataFrame


df = pd.read_csv("dataall.csv", header=None)

# Convert the DataFrame to a NumPy array
csv_values = df.values
# Load the surface we want to use as the background
# Read in your own background surface {'.ply', '.vtp', '.vtk', '.fs', '.asc', '.gii'} 
#surfaces = read_surface('./path/to/surface/file.gii.gz') 

# Or use one of the surface files included in the neuromaps package 
# Here we are using the 32k FsLR (a symmetric version 
# fsaverage space template with ~32k verticies in each hemisphere
surfaces = fetch_fslr()
lh, rh = surfaces['inflated']

# Next we want to load the parcellation/atlas we want to plot
# on the background surface. A parcellation is a array or surface

# Or use one of the surface files included in the brainspace package
atlas = load_parcellation('schaefer', 400, join = True)
# Convert the modified atlas array to a DataFrame
atlas_df = pd.DataFrame(atlas)

# Save the DataFrame as a CSV file
atlas_df.to_csv('atlas_index.csv', index=False, header=False)


# You can either plot this atlas directly, or assign new values 
# to each parcel to demonstrate an statistical effect. Here we assign a 

unique = np.unique(atlas[1:])
 # Extract unique elements excluding the first one

for i in range(1,401):
   # i = 10
    #rd = np.random.uniform(low=0.0, high=1.0, size=1).round(3)
    rd = csv_values[i-1]
    atlas = np.where(atlas == unique[i], rd, atlas)
    
max(atlas)
from matplotlib.colors import LinearSegmentedColormap

colors = ['#fef1f3','#992224']
cmap = LinearSegmentedColormap.from_list("camp",colors)
cmap
color_range = (0, 0.016)
atlas[atlas < 2e-3] = 0
# Generate plot
p = Plot(lh, rh, views=['lateral','medial'], zoom=1.2)
lh_parc, rh_parc = load_parcellation('schaefer')
p.add_layer({'left': lh_parc, 'right': rh_parc}, cmap='gray',
            as_outline=True, cbar=False)
p.add_layer(atlas, cbar=True, cmap=cmap,color_range=color_range)
p.add_layer(atlas, cmap='gray', as_outline=True, cbar=False)
kws = {'location': 'right', 'label_direction': 45, 'decimals': 10,
       'fontsize': 8, 'n_ticks':2,'shrink': .15, 'aspect':5,
       'draw_border': False}
fig = p.build(cbar_kws=kws)
save_path ='I:\\lda\\0611\\figure\\conj'
fig.savefig(save_path, format='tiff', dpi=600)
fig.show()
#p.build()

# Convert the modified atlas array to a DataFrame
atlas_df = pd.DataFrame(atlas)

# Save the DataFrame as a CSV file
atlas_df.to_csv('oatlas_modified.csv', index=False, header=False)





