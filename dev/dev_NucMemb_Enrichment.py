# Update system path to find pyfishquant
import sys
import os
import importlib
sys.path.insert(0,os.path.abspath(os.path.join('..','rna-loc')))
from rnaloc import LOCtoolbox

#%% Test function with entire analysis workflow 
importlib.reload(LOCtoolbox)
file_load =  'D:\\ongoing\\rna-loc\\NucMembEnrichment\\C1-N2_imb-2-670_NG-610_04_R3D__spots.txt'

results_all = LOCtoolbox.process_file(FQ_file=file_load, 
                        img_size=(533,541),
                        bin_prop=(0,90,20),
                        channels={'cells':'C4-'},
                        data_category={'roi':''},
                        annotation_extension='__RoiSet.zip',
                        img_extension='.tif',
                        show_plots = False,
                        Zrange=(0,0),
                        dZ = 0,
                        plot_callback=None,
                        progress_callback = None)

#%% BUILD ENTIRE WORKFLOW 
import sys
import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
import re
import os

# My Stuff
from rnaloc import annotationImporter, maskGenerator, FQtoolbox

from read_roi import read_roi_zip # https://github.com/hadim/read-roi
from read_roi import read_roi_file

#%% Read roi file
annote_name_full = 'C:\\Users\\muell\\Dropbox (Imod Pasteur)\\ImJoy\\DATA_SYNCED\\ImJoy\\rna-loc\\NucMembEnrichment\\zstack_segmentation\\C4-N2_imb-2-670_NG-610_04_R3D_Z64__RoiSet.zip'
roi_dict_complete = read_roi_zip(annote_name_full)


#%% Test workflow
importlib.reload(annotationImporter)
importlib.reload(maskGenerator)


## Open FQ results
file_load =  'D:\\ongoing\\rna-loc\\NucMembEnrichment\\C1-N2_imb-2-670_NG-610_04_R3D__spots.txt'


Zrange = (30,80)
dZ =3


fq_dict = FQtoolbox.read_FQ_matlab(file_load)
spots_all = FQtoolbox.get_rna(fq_dict)
Zrna = spots_all[:,[18]]

## Get folders
drive, path_and_file = os.path.splitdrive(file_load)
path_results, file_results = os.path.split(path_and_file)

## Import annotations
path_annot = os.path.join(path_results,'zstack_segmentation')
folderImporter = annotationImporter.FolderImporter(channels={'nuclei':'C4-'},
                                                   data_category={'roi':''},
                                                   annot_ext='__RoiSet.zip',
                                                   progress_callback=None)
annotDict = folderImporter.load(path_annot)
print('average roi size:', annotDict['roi_size'])


# Generate binary masks for a selected data-set
binaryGen = maskGenerator.BinaryMaskGenerator(erose_size=5,
                                              obj_size_rem=500,
                                              save_indiv=True,
                                              progress_callback=None)

# The generate function uses as an input the sub-dictionary for one data-category and one channel
annotatFiles = annotDict['roi']['nuclei']
maskDict = binaryGen.generate(annotatFiles)


# Use a loop and the update function to add the mask dictionary to the loaded annotation dictionary
print(maskDict.keys())
keys_delete = []
for k, v in annotatFiles.items():
    
    # Make sure that key in present in mask, otherwise delete
    if k in maskDict:
        v.update(maskDict[k])
        
    else:
        keys_delete.append(k) 
        
print(keys_delete)

# Bins of histogram
binsHist = np.arange(0,100,20)
width = 0.8 * (binsHist[1] - binsHist[0])
center = (binsHist[:-1] + binsHist[1:]) / 2

# Other parameters for calculation
dist_membr_RNA = np.array([])
dist_membr_pix = np.array([])
idx = 0

# Loop over all z-slices
hist_slice ={}
print(' == Loop over slices')
N_annot = len(annotatFiles)


for idx_file, (k_annot, v_annot) in enumerate(annotatFiles.items()):

    print(f'Slice: {k_annot}')
        
    # Get Z coordinate
    m = re.search('.*_Z([0-9]*)\.tif',k_annot)
    Zmask = int(m.group(1))

    # Check if outside of specified z range
    if Zrange is not None:
        if (Zmask <= Zrange[0]) or (Zmask >= Zrange[1]):
            print('Slice outside of range')
            continue

    # Get z-range for loop
    Zloop = np.logical_and(Zrna <= Zmask + dZ,Zrna >= Zmask - dZ).flatten()
    spots_loop = spots_all[Zloop,:]
    spots_loop_XY = spots_loop[:,[16, 17]].astype(int)

    # Distance transform
    dist_nuc_outside= ndimage.distance_transform_edt(~v_annot['mask_fill'])  
    dist_nuc_inside = ndimage.distance_transform_edt(v_annot['mask_fill'])  # Negate mask
    dist_nuc = dist_nuc_outside - dist_nuc_inside

    # Indices have to be inversed to access array
    dist_nuc_RNA_loop = dist_nuc[spots_loop_XY[:,0],spots_loop_XY[:,1]]

    # Get distance from membrane for all pixel in the cell
    mask_all = v_annot['mask_fill'] + v_annot['mask_edge']
    dist_membr_pix_loop = dist_nuc[mask_all]

    # Save values
    if idx == 0:
        dist_nuc_RNA = np.copy(dist_nuc_RNA_loop)
        dist_membr_pix = np.copy(dist_membr_pix_loop)
    else:
        dist_nuc_RNA = np.append(dist_nuc_RNA,dist_nuc_RNA_loop,axis=0)
        dist_membr_pix = np.append(dist_membr_pix,dist_membr_pix_loop,axis=0)
    idx+=1



#%%
binsHist = np.arange(-40,100,10)
width = 0.8 * (binsHist[1] - binsHist[0])
center = (binsHist[:-1] + binsHist[1:]) / 2

histRNA_all, bins = np.histogram(dist_nuc_RNA,binsHist ,density=False)


# Distance from nucleus" positive when outside, negative when inside
dist_nuc_outside= ndimage.distance_transform_edt(~v_annot['mask_fill'])  # Negate mask
dist_nuc_inside = ndimage.distance_transform_edt(v_annot['mask_fill'])  # Negate mask
dist_nuc = dist_nuc_outside - dist_nuc_inside


dist_membr_RNA_loop = dist_nuc[spots_loop_XY[:,0],spots_loop_XY[:,1]]


import matplotlib.pyplot as plt
fig1, ax = plt.subplots(1,3)
ax[0].imshow(v_annot['mask_fill'],cmap="hot")
ax[0].scatter(spots_loop_XY[:,1],spots_loop_XY[:,0],color='g',s=4)

ax[1].imshow(dist_nuc,cmap="hot")
ax[1].scatter(spots_loop_XY[:,1],spots_loop_XY[:,0],color='g',s=4)

ax[2].bar(center, histRNA_all, align='center', width=width)
ax[2].set_xlabel('Distance cell cortex [pix]')
ax[2].set_ylabel('# RNAs')
ax[2].set_xticks(center)
ax[2].set_xticklabels(center.astype(int))