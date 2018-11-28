# Update system path to find pyfishquant
import sys
sys.path.append('/Volumes/PILON_HD2/fmueller/Documents/code/ImJoy_dev/rna-loc/')
from rnaloc import LOCtoolbox

#%% Test function with entire analysis workflow 

file_load = '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/Segmentation__Dylan/2_ImJoy_membrane/img1/C1-img1__spots.txt'
file_load =  '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/Segmentation__Dylan/2_ImJoy_membrane/bug_img1/C1-img1__spots.txt'

results_all = LOCtoolbox.process_file(FQ_file=file_load, 
                        img_size=(960,960),
                        bin_prop=(0,90,20),
                        channels={'cells':'C3-'},
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
import annotationImporter, maskGenerator
import FQtoolbox

if  None:
    print('a')

#%% READ ANNOTATION DATA
importlib.reload(annotationImporter)
importlib.reload(maskGenerator)

# Load data with Folder importer
folderImporter = annotationImporter.FolderImporter(channels={'cells':'C3-'}, data_category={'roi':''},annot_ext ='__RoiSet.zip')

# Open folder
path_open  = '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/Segmentation__Dylan/2_ImJoy_membrane/img1/zstack_segmentation'
annotDict = folderImporter.load(path_open)
print('average roi size:', annotDict['roi_size'])

# Generate binary masks for a selected data-set
binaryGen = maskGenerator.BinaryMaskGenerator(erose_size=5, obj_size_rem=500, save_indiv=True)

# The generate function uses as an input the sub-dictionary for one data-category and one channel
annotatFiles = annotDict['roi']['cells']
maskDict     = binaryGen.generate(annotatFiles)

# Use a loop and the update function to add the mask dictionary to the loaded annotation dictonary\n",
for k, v in annotatFiles.items():
    v.update(maskDict[k])

## Print edge mask corresponding to first iterator item
data_key = next(iter(annotatFiles))
print(data_key)

fig, (ax1 ,ax2) = plt.subplots(1,2)
ax1.imshow(annotatFiles[data_key]['image'])
ax2.imshow(annotatFiles[data_key]['mask_edge'])

#%% What channel to analyze
channel = 1
channel = 2

#%% Perform analysis

## Open FQ results file
img_size = (960,960)

if channel == 1:
    file_open = '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/Segmentation__Dylan/2_ImJoy_membrane/img1/C1-erm-1_GFP_GFP-525_oEO158-erm-1-610_oEO160-set-3-670_11_R3D_D3D__spots.txt'

elif channel == 2:
    file_open = '/Volumes/PILON_HD2/fmueller/Documents/Data/ImJoy/Segmentation__Dylan/2_ImJoy_membrane/img1/C2-erm-1_GFP_GFP-525_oEO158-erm-1-610_oEO160-set-3-670_11_R3D_D3D__spots.txt'


## Prepare folder to save results
drive, path_and_file = os.path.splitdrive(file_open)
path, file = os.path.split(path_and_file)
file_base, ext = os.path.splitext(file)
    
path_save = os.path.join(path, '#analysis_MembraneDistance')
if not os.path.isdir(path_save):
    os.makedirs(path_save)


path_save_ch = os.path.join(path_save,'channel{}'.format(channel))
if not os.path.isdir(path_save_ch):
    os.makedirs(path_save_ch)

fq_dict = FQtoolbox.read_FQ_matlab(file_open)
spots_all = FQtoolbox.get_rna(fq_dict)

## Loop over all annotations

# bins of histogram
binsHist = np.arange(0,90,20)
width = 0.8 * (binsHist[1] - binsHist[0])
center = (binsHist[:-1] + binsHist[1:]) / 2

# Other parameters for calculation
Zrna = spots_all[:,[18]]
dist_membr_RNA = np.array([])
dist_membr_pix = np.array([])
idx = 0
dZ = 2

# Show new dictonar
for k, v in annotatFiles.items():
    
    # Get Z coordinate
    m = re.search('.*__Z([0-9]*)\.tif',k)
    Zmask = int(m.group(1))
    print(Zmask)
    Zloop = np.logical_and(Zrna <= Zmask + dZ,Zrna >= Zmask - dZ).flatten()
    spots_loop = spots_all[Zloop,:]
    spots_loop_XY = spots_loop[:,[16, 17]].astype(int)
        
    # Distance transform
    dist_membr = ndimage.distance_transform_edt(~v['mask_edge'])  # Negate mask   
    
    # Indices have to be inversed to access array
    dist_membr_RNA_loop = dist_membr[spots_loop_XY[:,0],spots_loop_XY[:,1]]
   
    # Get distance from membrane for all pixel in the cell
    mask_all = v['mask_fill'] + v['mask_edge']      
    dist_membr_pix_loop = dist_membr[mask_all]
    
    # Find min and max values for plotting
    pad = 10
    indMaskAx0 = np.argwhere(mask_all.sum(axis=0))
    minAx0 = indMaskAx0[0]-pad
    maxAx0 = indMaskAx0[-1]+pad
    
    indMaskAx1 = np.argwhere(mask_all.sum(axis=1))
    minAx1 = indMaskAx1[0]-pad
    maxAx1 = indMaskAx1[-1]+pad
    
    # Save values
    if idx == 0:
        dist_membr_RNA = np.copy(dist_membr_RNA_loop)   
        dist_membr_pix = np.copy(dist_membr_pix_loop)   
    else:
        dist_membr_RNA = np.append(dist_membr_RNA,dist_membr_RNA_loop,axis=0)
        dist_membr_pix = np.append(dist_membr_pix,dist_membr_pix_loop,axis=0)
    idx+=1
    
    
    #### Plot results
    
    # Set distance outside of cell to 0 for better plotting     
    dist_membr_plot = np.copy(dist_membr)
    dist_membr_plot[np.logical_not(mask_all)] = 0
    
    # Calculate histograms
    histRNA, bins = np.histogram(dist_membr_RNA_loop,binsHist ,density=False)
    histpix, bins = np.histogram(dist_membr_pix_loop,binsHist ,density=False)
       
    histRNAnorm = histRNA/histRNA.sum()
    histpixnorm = histpix/histpix.sum()
    
    histRNAnormPix = np.divide(histRNAnorm,histpixnorm)
    histRNAnormPix = np.nan_to_num(histRNAnormPix)
    
    # Generate plot
    fig1, ax = plt.subplots(2,3,num='C{}-Cell cortex analysis. Z={}'.format(channel,Zmask))
    fig1.set_size_inches((13,6))
    
    ax[0][0].imshow(v['image'],cmap="hot") 
    ax[0][0].get_xaxis().set_visible(False)
    ax[0][0].get_yaxis().set_visible(False)
    ax[0][0].set_xlim(minAx0, maxAx0)
    ax[0][0].set_ylim(minAx1, maxAx1)
    
    ax[0][1].imshow(mask_all,cmap="hot") 
    ax[0][1].get_xaxis().set_visible(False)
    ax[0][1].get_yaxis().set_visible(False)
    ax[0][1].set_xlim(minAx0, maxAx0)
    ax[0][1].set_ylim(minAx1, maxAx1)
    
    imgdum = ax[0][2].imshow(dist_membr_plot,cmap="hot") 
    ax[0][2].set_xlim(minAx0, maxAx0)
    ax[0][2].set_ylim(minAx1, maxAx1)
    
    ax[0][2].get_xaxis().set_visible(False)
    ax[0][2].get_yaxis().set_visible(False)   
    FQtoolbox.colorbar(imgdum)

    for kROI, vROI in v['roi'].items():
        roi_pos = vROI['pos']
        ax[0][2].plot(roi_pos[:,1],roi_pos[:,0],'b-')
        
    ax[0][2].scatter(spots_loop_XY[:,1],spots_loop_XY[:,0],color='g',s=4)


    ax[1][0].bar(center, histRNA, align='center', width=width)
    ax[1][0].set_xticks(center)  
    ax[1][0].set_xticklabels(center.astype(int))             
    ax[1][0].set_xlabel('Distance from cell cortex [pixel]')
    ax[1][0].set_ylabel('# RNAs')
      
    ax[1][1].bar(center, histpix, align='center', width=width)   
    ax[1][1].set_xticks(center)  
    ax[1][1].set_xticklabels(center.astype(int))   
    ax[1][1].set_xlabel('Distance from cell cortex [pixel]')
    ax[1][1].set_ylabel('# pixels')  

    ax[1][2].bar(center, histRNAnormPix, align='center', width=width)   
    ax[1][2].set_xticks(center)  
    ax[1][2].set_xticklabels(center.astype(int))   
    ax[1][2].set_xlabel('Distance from cell cortex [pixel]')
    ax[1][2].set_ylabel('Renormalized frequency')  

    # Set titles
    ax[0][0].title.set_text('Cell cortex')     
    ax[0][1].title.set_text('Cell mask')       
    ax[0][2].title.set_text('Distance transform')     
    
    ax[1][0].title.set_text('RNAs')     
    ax[1][1].title.set_text('All pixel')       
    ax[1][2].title.set_text('Renormalized RNA distance')   
    
    plt.tight_layout()
    
    # Save
    path_save_ch = os.path.join(path_save,'channel{}'.format(channel))
    file_save ='C{}_Z{}.png'.format(channel,Zmask)
    plt.savefig(os.path.join(path_save_ch, file_save),dpi=200)
    plt.close()
  
#### Save results
    
if channel == 1:    
    C1_RNAcounts = np.copy(dist_membr_RNA)   
    C1_pixcounts = np.copy(dist_membr_pix)   
elif channel == 2: 
    C2_RNAcounts = np.copy(dist_membr_RNA)   
    C2_pixcounts = np.copy(dist_membr_pix)   
      
    
#%% Make histograms for both channels

histC1RNA, bins = np.histogram(C1_RNAcounts,binsHist ,density=False)
histC1RNAnorm = histC1RNA/histC1RNA.sum()

histC2RNA, bins = np.histogram(C2_RNAcounts,binsHist ,density=False)
histC2RNAnorm = histC2RNA/histC2RNA.sum()

# Renormalize with pixel counts
histpix, bins = np.histogram(C1_pixcounts,binsHist ,density=False)
histpixnorm = histpix/histpix.sum()

histC1RNAnormPix = np.divide(histC1RNAnorm,histpixnorm)
histC1RNAnormPix = np.nan_to_num(histC1RNAnormPix)
    
histC2RNAnormPix = np.divide(histC2RNAnorm,histpixnorm)
histC2RNAnormPix = np.nan_to_num(histC2RNAnormPix)


# Plot results

fig1, ax = plt.subplots(2,3,num='Distance comparison')
fig1.set_size_inches((11,7))

ax[0][0].bar(center, histC1RNA, align='center', width=width)
ax[0][0].set_xlabel('Distance cell cortex [pix]')
ax[0][0].set_ylabel('# RNAs')
ax[0][0].set_xticks(center)  
ax[0][0].set_xticklabels(center.astype(int))  

ax[0][1].bar(center, histC2RNA, align='center', width=width)
ax[0][1].set_xlabel('Distance  cell cortex [pix]')
ax[0][1].set_ylabel('# RNAs')
ax[0][1].set_xticks(center)  
ax[0][1].set_xticklabels(center.astype(int))  
  
p1=ax[0][2].bar(center, histC1RNAnorm, align='center', width=width/2)  
p2=ax[0][2].bar(center+width/2, histC2RNAnorm, align='center', width=width/2)  
ax[0][2].set_xlabel('Distance cell cortex [pix]')
ax[0][2].set_ylabel('Frequency')
ax[0][2].set_xticks(center+width/4)
ax[0][2].set_xticklabels(center.astype(int))  
ax[0][2].legend((p1[0], p2[0]), ('C1', 'C2'),loc=3)

ax[1][0].bar(center, histC1RNAnormPix, align='center', width=width)
ax[1][0].set_xlabel('Distance cell cortex [pix]')
ax[1][0].set_ylabel('Renormalized counts')
ax[1][0].set_xticks(center)  
ax[1][0].set_xticklabels(center.astype(int))  

ax[1][1].bar(center, histC2RNAnormPix, align='center', width=width)
ax[1][1].set_xlabel('Distance cell cortex [pix]')
ax[1][1].set_ylabel('Renormalized counts')
ax[1][1].set_xticks(center)  
ax[1][1].set_xticklabels(center.astype(int))  

p1=ax[1][2].bar(center, histC1RNAnormPix, align='center', width=width/2)  
p2=ax[1][2].bar(center+width/2, histC2RNAnormPix, align='center', width=width/2)  
ax[1][2].set_xlabel('Distance cell cortex [pix]')
ax[1][2].set_ylabel('Renormalized counts')
ax[1][2].set_xticks(center+width/4)
ax[1][2].set_xticklabels(center.astype(int))  
ax[1][2].legend((p1[0], p2[0]), ('C1', 'C2'),loc=3)


ax[0][0].title.set_text('C1-absolute counts')     
ax[0][1].title.set_text('C2-absolute counts')       
ax[0][2].title.set_text('COMPARE [frequency]')     
    
ax[1][0].title.set_text('C1-renormalized')     
ax[1][1].title.set_text('C2-renormalized')       
ax[1][2].title.set_text('COMPARE [renorm]')  

plt.tight_layout()
    
plt.savefig(os.path.join(path_save, 'CellCortexDist.png'),dpi=200)
#plt.close()