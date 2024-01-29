# Summary

This tutorial demonstrates the workflow for generating the semi-synthetic multi-sample ST data from the STARmap, Smart-seq mouse visual cortex data.

## Dataset

The STARmap mouse visual cortex is publicly available but for sake of convenience, all the required data for semi-synthetic ST data generation can be downloaded from the link below: 

- [STARmap mouse visual cortex data](https://drive.google.com/file/d/1DUoTw2VSxT5a5uOLSoyzk5bTkdT3Cuwm/view?usp=sharing)


## Semi-synthetic Data generation

Here, we describe how one can generate the simulated low resolution semi-synthetic ST data from the STARmap mouse visual cortex data. As an example here we generate one sample with the window size of 600 pixels:

```python
"""
Semi-synthetic Data generation for MUSTANG Benchmarking.

Acknowledgements:
Based on the descriptions and codes from the 
https://github.com/QuKunLab/SpatialBenchmarking
"""


import numpy as np
import pandas as pd
import sys
import pickle
import os
import time as tm
from functools import partial
import scipy.stats as st
from scipy.stats import wasserstein_distance
import scipy.stats
import copy
from sklearn.model_selection import KFold
import pandas as pd
import multiprocessing
import matplotlib as mpl 
import matplotlib.pyplot as plt
import scanpy as sc
import warnings
import subprocess
import seaborn as sns
import os

warnings.filterwarnings('ignore')
import time
from scipy.spatial import distance_matrix
from sklearn.metrics import matthews_corrcoef
from scipy import stats

from scipy.spatial.distance import cdist
import h5py

import anndata
import torch
import sys

import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.serif'] = ['Arial']
matplotlib.rcParams['figure.dpi'] = 60

def Simulated(spatial_rna, spatial_meta, spatial_loc, CoordinateXlable, CoordinateYlable, window, outdir):
    if os.path.exists(outdir):
        print ('The output file is in ' + outdir)
    else:
        os.mkdir(outdir)
    combined_spot = []
    combined_spot_loc = []
    window=window
    c = 0
    for x in np.arange((spatial_loc[CoordinateXlable].min()//window),spatial_loc[CoordinateXlable].max()//window+1):
        for y in np.arange((spatial_loc[CoordinateYlable].min()//window),spatial_loc[CoordinateYlable].max()//window+1):
            tmp_loc = spatial_loc[(x*window < spatial_loc[CoordinateXlable]) & (spatial_loc[CoordinateXlable] < (x+1)*window) & (y*window < spatial_loc[CoordinateYlable]) & (spatial_loc[CoordinateYlable] < (y+1)*window)]
            if len(tmp_loc) > 0:
                c += 1
                combined_spot_loc.append([x,y])
                combined_spot.append(tmp_loc.index.to_list())
            
    combined_cell_counts = pd.DataFrame([len(s) for s in combined_spot],columns=['cell_count'])
    combined_cell_counts.to_csv(outdir + '/combined_cell_counts.txt',sep='\t')
    combined_cell_counts = pd.read_csv(outdir + '/combined_cell_counts.txt',sep='\t',index_col=0)
    print ('The simulated spot has cells with ' + str(combined_cell_counts.min()[0]) + ' to ' + str(combined_cell_counts.max()[0]))
    combined_spot_loc = pd.DataFrame(combined_spot_loc, columns=['x','y'])
    combined_spot_loc.to_csv(outdir + '/combined_Locations.txt',sep='\t',index=False)


    combined_spot_exp = []
    for s in combined_spot:
        combined_spot_exp.append(spatial_rna.loc[s,:].sum(axis=0).values)
    combined_spot_exp = pd.DataFrame(combined_spot_exp, columns=spatial_rna.columns)
    combined_spot_exp.to_csv(outdir + '/combined_spatial_count.txt',sep='\t',index=False)

    combined_spot_clusters = pd.DataFrame(np.zeros((len(combined_spot_loc.index),len(np.unique(spatial_meta['celltype'])))),columns=np.unique(spatial_meta['celltype']))
    for i,c in enumerate(combined_spot):
        for clt in spatial_meta.loc[c,'celltype']:
            combined_spot_clusters.loc[i,clt] += 1
    combined_spot_clusters.to_csv(outdir + '/combined_spot_clusters.txt',sep='\t')
    print ('The simulated spot has size ' + str(combined_spot_clusters.shape[0]))


time_start=time.time()
PATH = 'Dataset10/'
sc_rna = pd.read_csv(PATH + 'scRNA_count.txt', sep='\t',index_col=0,engine='c',low_memory=False)
spatial_rna = pd.read_csv(PATH + 'Spatial_count.txt',sep='\t')
spatial_rna = spatial_rna.loc[:,spatial_rna.columns[spatial_rna.columns.isin(sc_rna.index)]]
spatial_meta = pd.read_csv(PATH + 'Spatial_annotate.txt',sep='\t',index_col=0)
spatial_loc = pd.read_csv(PATH + 'Locations.txt',sep='\t')

window = 600
CoordinateX = 'X'
CoordinateY = 'Y'
outdir = PATH + 'Simulated_STARmap/'

Simulated(spatial_rna, spatial_meta, spatial_loc, CoordinateX, CoordinateY, window, outdir) 
time_end=time.time()
print('STARmap dataset simulated step costs',time_end-time_start,'s')


spatial_meta = pd.read_csv(PATH + 'Spatial_annotate.txt',sep='\t')
spatial_loc = pd.read_csv(PATH + 'Locations.txt',sep='\t')
fig,ax = plt.subplots(figsize=(4,7))
cmap = sns.color_palette('tab20',n_colors=len(np.unique(spatial_meta.celltype)))
for i,c in enumerate(np.unique(spatial_meta.celltype)):
    ax.scatter(x=spatial_loc[spatial_meta['celltype'] == c]['X'],y=spatial_loc[spatial_meta['celltype'] == c]['Y'],c=matplotlib.colors.to_hex(cmap[i]),label=c,s=30,marker='o',edgecolors='none')
ax.set_title('gd')
ax.legend(bbox_to_anchor=(1,1))
plt.xticks(np.arange(0,7500,600))
plt.yticks(np.arange(0,14250,600))
plt.xlim(0,7500)
plt.ylim(0,14250)
plt.grid()
plt.show()

```

<img src="https://github.com/namini94/MUSTANG/blob/main/Miscel/Mouse_Brain_Markdown_Figs/MB_NoBatch.png" width="50%" height="50%">
