# Update system path to find pyfishquant
import sys
import importlib
sys.path.append('/Volumes/PILON_HD2/fmueller/Documents/code/ImJoy_dev/rna-loc/')
from rnaloc import LOCtoolbox

#%% Test function with entire analysis workflow 
importlib.reload(LOCtoolbox)
file_load =  '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/rna-loc/NucMembEnrichment/C1-N2_imb-2-670_NG-610_04_R3D__spots.txt'

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
annote_name_full = '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/rna-loc/NucMembEnrichment/zstack_segmentation/C4-N2_imb-2-670_NG-610_04_R3D_Z64__RoiSet.zip'
roi_dict_complete = read_roi_zip(annote_name_full)




annote_name_full = '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/rna-loc/NucMembEnrichment/zstack_segmentation/0250-0274.roi'
roi_dict_complete1 = read_roi_file(annote_name_full)

#%% Test workflow
importlib.reload(annotationImporter)
importlib.reload(maskGenerator)


## Open FQ results
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
                                              progress_callback=progress_callback)

# The generate function uses as an input the sub-dictionary for one data-category and one channel
annotatFiles = annotDict['roi']['cells']
maskDict = binaryGen.generate(annotatFiles)