# Update system path to find pyfishquant
import sys
import os
import importlib
sys.path.insert(0,os.path.abspath(os.path.join('..','rna-loc')))
from rnaloc import LOCtoolbox
import numpy as np

#%% Test function with entire analysis workflow 

importlib.reload(LOCtoolbox)
file_load =  '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/rna-loc/NucMembEnrichment/C1-N2_imb-2-670_NG-610_04_R3D_res_GMM.txt'
Zrange = (0,0)
dZ=0
binsHist = np.fromstring('-100,-50,-10,10,50,100,300',dtype=int, sep=',')


results_all = LOCtoolbox.calc_nuclear_enrichment(FQ_file=file_load, 
                        binsHist=binsHist,
                        show_plots = False,
                        Zrange=Zrange,
                        dZ = dZ,
                        plot_callback=None,
                        progress_callback = None)

#%% BUILD ENTIRE WORKFLOW 
import sys
import matplotlib.pyplot as plt
from scipy import ndimage
import re
import os
from read_roi import read_roi_file
# My Stuff
from rnaloc import annotationImporter, maskGenerator, FQtoolbox


#%% Test function

def calc_nuclear_enrichment(FQ_file,binsHist):

    
    ## Get folder
    drive, path_and_file = os.path.splitdrive(FQ_file)
    path_results, file_results = os.path.split(path_and_file)
    
    
    # *************************************************************************
    # Load nuclear outline
    
    # Generate binary masks for a selected data-set
    binaryGen = maskGenerator.BinaryMaskGenerator(erose_size=5,
                                                  obj_size_rem=500,
                                                  save_indiv=True,
                                                  progress_callback=None)  
    
    ## Import annotations for nuclei
    path_annot = os.path.join(drive,path_results,'zstack_segmentation')
    folderImporter = annotationImporter.FolderImporter(channels={'nuclei':'C4-'},
                                                       data_category={'roi':''},
                                                       annot_ext='__RoiSet.zip',
                                                       progress_callback=None)
    annotDict = folderImporter.load(path_annot)
    print('average roi size:', annotDict['roi_size'])
    
    

    # The generate function uses as an input the sub-dictionary for one data-category and one channel
    annotatFiles = annotDict['roi']['nuclei']
    mask_dict_nuclei = binaryGen.generate(annotatFiles)
    
    keys_delete = []
    for k, v in annotatFiles.items():
        
        # Make sure that key in present in mask, otherwise delete
        if k in mask_dict_nuclei:
            v.update(mask_dict_nuclei[k])
            
        else:
            keys_delete.append(k) 
    
    # *************************************************************************
    #  Load embryo outline
    file_embryo = os.path.join(drive,path_results,'embryo_contour.roi')
    
    # Conver dictionary & get size of ROIS
    roi_import = read_roi_file(file_embryo)
    roi_dict = {}
    roi_dict['embryo'] = {}
    roi_dict['embryo'] ['type'] = roi_import['embryo_contour']['type']
    roi_dict['embryo'] ['pos'] = np.column_stack((roi_import['embryo_contour']['y'], roi_import['embryo_contour']['x']))
    
    # Assemble structure to call mask generator
    image_size = annotatFiles.values().__iter__().__next__()['image'].shape
    image_fake = np.zeros(image_size, dtype=np.uint8)
    
    annotat_files_embryo = {} 
    annotat_files_embryo['embryo_contour'] = {}   
    annotat_files_embryo['embryo_contour']['roi'] = roi_dict
    annotat_files_embryo['embryo_contour']['image'] = image_fake
    
    
    mask_dict_embryo = binaryGen.generate(annotat_files_embryo)
    mask_embryo = mask_dict_embryo['embryo_contour']['mask_fill']
        
    
    # *************************************************************************
    #  Load and analyze FQ results
    
    fq_dict = FQtoolbox.read_FQ_matlab(file_load)
    spots_all = FQtoolbox.get_rna(fq_dict)
    
    # Z position in pixel
    Zrna = np.divide(spots_all[:,[2]],fq_dict['settings']['microscope']['pix_z']).astype(int)
    
    # Other parameters for calculation
    dist_membr_RNA = np.array([])
    dist_membr_pix = np.array([])
    idx = 0
    
    # Loop over all z-slices
    print(' == Loop over slices')
    
    for idx_file, (k_annot, v_annot) in enumerate(annotatFiles.items()):
    
        print(f'Slice: {k_annot}')
            
        # Get Z coordinate
        m = re.search('.*_Z([0-9]*)\.tif',k_annot)
        Zmask = int(m.group(1))
    
        # Check if outside of specified z range
        if Zrange is not None:
            if (Zmask < Zrange[0]) or (Zmask > Zrange[1]):
                print(f'Z-slice outside of specified range: {Zmask}')
                continue
    
        # Get z-range for loop
        Zloop = np.logical_and(Zrna <= Zmask + dZ,Zrna >= Zmask - dZ).flatten()
        Zloop = (Zrna == Zmask).flatten()
        spots_loop = spots_all[Zloop,:]
        spots_loop_XY = np.divide(spots_loop[:,[0,1]],fq_dict['settings']['microscope']['pix_xy']).astype(int)
    
    
        # Distance transform
        dist_nuc_outside= ndimage.distance_transform_edt(~v_annot['mask_fill'])  
        dist_nuc_inside = ndimage.distance_transform_edt(v_annot['mask_fill'])  # Negate mask
        dist_nuc = dist_nuc_outside - dist_nuc_inside
    
        # Indices have to be inversed to access array
        dist_nuc_RNA_loop = dist_nuc[spots_loop_XY[:,0],spots_loop_XY[:,1]]
    
        # Get distance from membrane for all pixel in the cell
        dist_membr_pix_loop = dist_nuc[mask_embryo.astype(bool)]
    
        # Save values
        if idx == 0:
            dist_membr_RNA = np.copy(dist_nuc_RNA_loop)
            dist_membr_pix = np.copy(dist_membr_pix_loop)
        else:
            dist_membr_RNA = np.append(dist_membr_RNA,dist_nuc_RNA_loop,axis=0)
            dist_membr_pix = np.append(dist_membr_pix,dist_membr_pix_loop,axis=0)
        idx+=1

    # *************************************************************************
    #  Load and analyze FQ results
    xTicks = binsHist[:-1]
    width = 0.9*np.diff(binsHist)
    width_full = np.diff(binsHist)
    center = (binsHist[:-1] + binsHist[1:]) / 2
    
    histRNA_all, bins = np.histogram(dist_membr_RNA,binsHist ,density=False)
    histPIX_all, bins = np.histogram(dist_membr_pix,binsHist ,density=False)
    histRNA_norm = np.divide(histRNA_all,histPIX_all)
    counts_total = np.nansum(np.multiply(histRNA_norm,width_full))
    histRNA_norm = np.divide(histRNA_norm,counts_total)
    
    fig1, ax = plt.subplots(3,1)
    ax[0].bar(center, histRNA_all, align='center', width=width)
    ax[0].set_xlabel('Dist to nuc')
    ax[0].set_ylabel('# RNAs')
    ax[0].set_xticks(xTicks)
    ax[0].set_xticklabels(xTicks.astype(int))
    
    ax[1].bar(center, histPIX_all, align='center', width=width)
    ax[1].set_xlabel('Dist to nuc')
    ax[1].set_ylabel('# pixels')
    ax[1].set_xticks(xTicks)
    ax[1].set_xticklabels(xTicks.astype(int))
    
    ax[2].bar(center, histRNA_norm, align='center', width=width)
    ax[2].set_xlabel('Distance to nuc')
    ax[2].set_ylabel('RNA counts [a.u.]')
    ax[2].set_xticks(xTicks)
    ax[2].set_xticklabels(xTicks.astype(int))
    
    ax[0].title.set_text('RNAs')
    ax[1].title.set_text('All pixel')
    ax[2].title.set_text('RNA renormalized with pixels')
        
    plt.tight_layout()
    

hist_all = np.stack((center,count_RNA_normSum,count_RNA,count_pix),axis=1)    
    

hist_plot_all = {'width':width,'center':center, 'bins':bins,
                 'histRNA_all':histRNA_all,
                 'histRNA_all_norm':histRNA_all_norm,
                 'histpix_all_norm':histpix_all_norm,
                 'hist_RNA_all_normPix':hist_RNA_all_normPix}

#%% ===========================================================================
#   TEST WORKFLOW ELEMENTS
#   ===========================================================================   
    
#%% Read nuclear outlines
importlib.reload(annotationImporter)
importlib.reload(maskGenerator)


## Open FQ results
file_load =  '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/rna-loc/NucMembEnrichment/C1-N2_imb-2-670_NG-610_04_R3D_res_GMM.txt'

## Get folders
drive, path_and_file = os.path.splitdrive(file_load)
path_results, file_results = os.path.split(path_and_file)

## Import annotations
path_annot = os.path.join(drive,path_results,'zstack_segmentation')
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
        

#%% Determine embryo outline
from read_roi import read_roi_file
file_embryo = os.path.join(drive,path_results,'embryo_contour.roi')


# Conver dictionary & get size of ROIS
roi_import = read_roi_file(file_embryo)
roi_dict = {}
roi_dict['embryo'] = {}
roi_dict['embryo'] ['type'] = roi_import['embryo_contour']['type']
roi_dict['embryo'] ['pos'] = np.column_stack((roi_import['embryo_contour']['y'], roi_import['embryo_contour']['x']))

# Assemble structure to call mask generator
image_size = annotatFiles.values().__iter__().__next__()['image'].shape
image_fake = np.zeros(image_size, dtype=np.uint8)

annotat_files_embryo = {} 
annotat_files_embryo['embryo_contour'] = {}   
annotat_files_embryo['embryo_contour']['roi'] = roi_dict
annotat_files_embryo['embryo_contour']['image'] = image_fake


mask_dict_embryo = binaryGen.generate(annotat_files_embryo)
mask_embryo = mask_dict_embryo['embryo_contour']['mask_fill']

    
#%% Analyze FQ results
file_load =  '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/rna-loc/NucMembEnrichment/C2-N2_imb-2-670_NG-610_04_R3D_res_GMM.txt'

fq_dict = FQtoolbox.read_FQ_matlab(file_load)
spots_all = FQtoolbox.get_rna(fq_dict)

# Z position in pixel
Zrna = np.divide(spots_all[:,[2]],fq_dict['settings']['microscope']['pix_z']).astype(int)

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
        if (Zmask < Zrange[0]) or (Zmask > Zrange[1]):
            print(f'Z-slice outside of specified range: {Zmask}')
            continue

    # Get z-range for loop
    Zloop = np.logical_and(Zrna <= Zmask + dZ,Zrna >= Zmask - dZ).flatten()
    Zloop = (Zrna == Zmask).flatten()
    spots_loop = spots_all[Zloop,:]
    #spots_loop_XY = spots_loop[:,[16, 17]].astype(int)
    spots_loop_XY = np.divide(spots_loop[:,[0,1]],fq_dict['settings']['microscope']['pix_xy']).astype(int)


    # Distance transform
    dist_nuc_outside= ndimage.distance_transform_edt(~v_annot['mask_fill'])  
    dist_nuc_inside = ndimage.distance_transform_edt(v_annot['mask_fill'])  # Negate mask
    dist_nuc = dist_nuc_outside - dist_nuc_inside

    # Indices have to be inversed to access array
    dist_nuc_RNA_loop = dist_nuc[spots_loop_XY[:,0],spots_loop_XY[:,1]]

    # Get distance from membrane for all pixel in the cell
    mask_all = v_annot['mask_fill'] + v_annot['mask_edge']
    dist_membr_pix_loop = dist_nuc[mask_embryo.astype(bool)]

    # Save values
    if idx == 0:
        dist_membr_RNA = np.copy(dist_nuc_RNA_loop)
        dist_membr_pix = np.copy(dist_membr_pix_loop)
    else:
        dist_membr_RNA = np.append(dist_membr_RNA,dist_nuc_RNA_loop,axis=0)
        dist_membr_pix = np.append(dist_membr_pix,dist_membr_pix_loop,axis=0)
    idx+=1


#%% Calculate histograms and plot
    
binsHist = np.array([-100,-50,-10,10,50,100,300])
xTicks = binsHist[:-1]
width = 0.9*np.diff(binsHist)
width_full = np.diff(binsHist)
#width = 0.8 * (binsHist[1] - binsHist[0])
center = (binsHist[:-1] + binsHist[1:]) / 2

histRNA_all, bins = np.histogram(dist_membr_RNA,binsHist ,density=False)
histPIX_all, bins = np.histogram(dist_membr_pix,binsHist ,density=False)
histRNA_norm = np.divide(histRNA_all,histPIX_all)
counts_total = np.nansum(np.multiply(histRNA_norm,width_full))
histRNA_norm = np.divide(histRNA_norm,counts_total)

fig1, ax = plt.subplots(3,1)
ax[0].bar(center, histRNA_all, align='center', width=width)
ax[0].set_xlabel('Dist to nuc')
ax[0].set_ylabel('# RNAs')
ax[0].set_xticks(xTicks)
ax[0].set_xticklabels(xTicks.astype(int))

ax[1].bar(center, histPIX_all, align='center', width=width)
ax[1].set_xlabel('Dist to nuc')
ax[1].set_ylabel('# pixels')
ax[1].set_xticks(xTicks)
ax[1].set_xticklabels(xTicks.astype(int))

ax[2].bar(center, histRNA_norm, align='center', width=width)
ax[2].set_xlabel('Distance to nuc')
ax[2].set_ylabel('RNA counts [a.u.]')
ax[2].set_xticks(xTicks)
ax[2].set_xticklabels(xTicks.astype(int))

ax[0].title.set_text('RNAs')
ax[1].title.set_text('All pixel')
ax[2].title.set_text('RNA renormalized with pixels')
    
plt.tight_layout()



#%% [ALTERNATIVE to embryo outline] Complex hull around points
#  Alternative to an outline of the embryo

from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt


# Get all points in z-range
Zind = np.logical_and(Zrna <=  Zrange[1],Zrna >=  Zrange[0]).flatten()
spots_sel = spots_all[Zind,:]
spots_XY = np.divide(spots_sel[:,[0,1]],fq_dict['settings']['microscope']['pix_xy']).astype(int)

hull = ConvexHull(spots_XY)


# Draw polygon
image_size=(541,533)
from skimage import draw as drawSK
rr, cc = drawSK.polygon(spots_XY[hull.vertices,0], spots_XY[hull.vertices,1])

# Make sure it's not outside
rr[rr < 0] = 0
rr[rr > image_size[0] - 1] = image_size[0] - 1

cc[cc < 0] = 0
cc[cc > image_size[0] - 1] = image_size[0] - 1

mask_fill_hull = np.zeros(image_size, dtype=np.uint8)
mask_fill_hull[rr, cc] = 1


fig1, ax = plt.subplots(1,2)
ax[0].plot(spots_XY[:,0], spots_XY[:,1], 'o')
ax[0].plot(spots_XY[hull.vertices,0], spots_XY[hull.vertices,1], 'r--', lw=2)
ax[0].plot(spots_XY[hull.vertices[0],0], spots_XY[hull.vertices[0],1], 'ro')

ax[1].imshow(mask_fill_hull,cmap="hot")  
ax[1].plot(spots_XY[hull.vertices,1],spots_XY[hull.vertices,0],  'r--', lw=2)
ax[1].plot(spots_XY[hull.vertices[0],1],spots_XY[hull.vertices[0],0],  'ro')
