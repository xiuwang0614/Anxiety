# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 15:04:02 2024

@author: Lenovo
"""

from surfplot import Plot
from brainspace.datasets import load_parcellation
from brainspace.mesh.mesh_io import read_surface
from neuromaps.datasets import fetch_fslr
import numpy as np
import pandas as pd
import os
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Specify the directory path you want to set as the working directory
directory_path = 'I:\\lda\\mask'

# Set the working directory
os.chdir(directory_path)

# Read the CSV file into a DataFrame
df = pd.read_csv("data.csv", header=None)

# Convert the DataFrame to a NumPy array
csv_values = df.values

# Load the surface we want to use as the background
surfaces = fetch_fslr()
lh, rh = surfaces['inflated']

# Load the parcellation/atlas we want to plot on the background surface
atlas = load_parcellation('schaefer', 400, join=True)

# Convert the modified atlas array to a DataFrame
atlas_df = pd.DataFrame(atlas)

# Save the DataFrame as a CSV file
atlas_df.to_csv('atlas_index.csv', index=False, header=False)

# Assign new values to each parcel
unique = np.unique(atlas[1:])
for i in range(1, 401):
    rd = csv_values[i-1]
    atlas = np.where(atlas == unique[i], rd, atlas)

# Define a custom colormap
colors = ['#e7e4f8','#a6badf','#6892c5','#4675b0','#28437d']
cmap = LinearSegmentedColormap.from_list("camp", colors)

# Define the color range (1 to 28)
color_range = (1, 28)

# Define corresponding labels
labels = ['Vis','SMN',
          'DorsAttnPost','DorsAttnFEF','DorsAttnPrCv',
          'SalVentAttnParOper','SalVentAttnTempOcc','SalVentAttnFrOperIns','SalVentAttnPFCl','SalVentAttnMed','SalVentAttnPrC','SalVentAttnTempOccPar',
          'LimbicOFC','LimbicTempPole',
          'ContPar','ContTemp','ContOFC','ContPFCl','ContPFCv','ContpCun','ContCing','ContPFCmp',
          'DefaultTemp','DefaultPar','DefaultPFC','DefaultPFCv','DefaultPFCdPFCm','DefaultpCunPCC']

# Generate plot
p = Plot(lh, rh, views=['lateral','medial'], zoom=1.2)

# Add layers to the plot
lh_parc, rh_parc = load_parcellation('schaefer')
p.add_layer({'left': lh_parc, 'right': rh_parc}, cmap='gray',
            as_outline=True, cbar=False)
p.add_layer(atlas, cbar=False, cmap=cmap, color_range=color_range)
p.add_layer(atlas, cmap='gray', as_outline=True, cbar=False)

# Customize colorbar
kws = {
    'location': 'right', 
    'label_direction': 45, 
    'decimals': 0,
    'fontsize': 8, 
    'n_ticks': 28,
    'shrink': .9,  # Adjusted to make the colorbar longer
    'aspect': 20,  # Adjusted to make the colorbar thinner
    'draw_border': False
}
fig = p.build(cbar_kws=kws)

# Create legend entries (colored squares and corresponding labels)
legend_patches = []
for i, label in enumerate(labels):
    color = cmap(i / (len(labels) - 1))  # Get the color from the colormap
    patch = mpatches.Patch(color=color, label=label)
    legend_patches.append(patch)

# Add legend to the plot
fig.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1.05, 0.5),
           title="Parcellations", title_fontsize='13', fontsize='10', frameon=False)

# Show the plot
fig.show()

fig.savefig("D:\\Desktop", format='tiff', dpi=600)
fig.show()