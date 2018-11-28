#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:43:47 2018

@author: fmueller
"""

# Imports
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

import os
from rnaloc import annotationImporter
from rnaloc import maskGenerator
from rnaloc import FQtoolbox
import numpy as np
import re
from scipy import ndimage
import json
import time

# JSON encoder for numpy
# From https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


# Process specified FQ file

def process_file(FQ_file, img_size = (960,960), bin_prop = (0,90,20), channels={'cells':'C3-'},data_category={'roi':''},annotation_extension ='__RoiSet.zip',img_extension='.tif',show_plot=False,Zrange=None,dZ=2):
    '''
    Function uses annotations generated in FIJI and creates mask based
    on the specified parameters. The resulting files are zipped and be
    used for training of a neural network with ImJoy.

    Args:


        Zrange [tuple, 2 elements]. [Optional] Tuple specifying minimum and maximum z-value that is considered
        in analysis.

        bin_prop [Tuple, 3 elements]. Specifies the bins for the histograms (min, max,delta).

    '''
    # Get input args. Has to be FIRST call!
    input_args = locals()

    # Make sure input args are correct - assignments with 0 can come from ImJoy
    if Zrange[0] ==0 or Zrange[1] ==0:
        Zrange = None

    if bin_prop[1] == 0 or bin_prop[1] == 0:
        bin_prop = (0,90,20)

    ## Prepare folder to save results
    drive, path_and_file = os.path.splitdrive(FQ_file)
    path_results, file_results = os.path.split(path_and_file)
    file_base, ext = os.path.splitext(file_results)

    path_save = os.path.join(path_results, file_base, 'MembDist_{}'.format(time.strftime("%y%m%d-%H%M", time.localtime())))
    if not os.path.isdir(path_save):
        os.makedirs(path_save)

    ## Open FQ results
    fq_dict = FQtoolbox.read_FQ_matlab(FQ_file)
    spots_all = FQtoolbox.get_rna(fq_dict)
    Zrna = spots_all[:,[18]]

    # Open annotations
    print(' == Open annotations')
    if 'RoiSet.zip' in annotation_extension:

        path_annot = os.path.join(path_results,'zstack_segmentation')
        folderImporter = annotationImporter.FolderImporter(channels=channels,
                                                           data_category=data_category,
                                                           annot_ext=annotation_extension)
        annotDict = folderImporter.load(path_annot)
        print('average roi size:', annotDict['roi_size'])

        # Generate binary masks for a selected data-set
        binaryGen = maskGenerator.BinaryMaskGenerator(erose_size=5,
                                                      obj_size_rem=500,
                                                      save_indiv=True)

        # The generate function uses as an input the sub-dictionary for one data-category and one channel
        annotatFiles = annotDict['roi']['cells']
        maskDict = binaryGen.generate(annotatFiles)

        #Use a loop and the update function to add the mask dictionary to the loaded annotation dictonary\n",
        for k, v in annotatFiles.items():
            v.update(maskDict[k])

    # Bins of histogram
    binsHist = np.arange(bin_prop[0],bin_prop[1],bin_prop[2])
    width = 0.8 * (binsHist[1] - binsHist[0])
    center = (binsHist[:-1] + binsHist[1:]) / 2

    # Other parameters for calculation
    dist_membr_RNA = np.array([])
    dist_membr_pix = np.array([])
    idx = 0

    # Loop over all z-slices
    hist_slice ={}
    print(' == Loop over slices')
    for k, v in annotatFiles.items():

        print(f'Slice: {k}')
        # Get Z coordinate
        m = re.search('.*_Z([0-9]*)\.tif',k)
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
        dist_membr = ndimage.distance_transform_edt(~v['mask_edge'])  # Negate mask

        # Indices have to be inversed to access array
        dist_membr_RNA_loop = dist_membr[spots_loop_XY[:,0],spots_loop_XY[:,1]]

        # Get distance from membrane for all pixel in the cell
        mask_all = v['mask_fill'] + v['mask_edge']
        dist_membr_pix_loop = dist_membr[mask_all]

        # Save values
        if idx == 0:
            dist_membr_RNA = np.copy(dist_membr_RNA_loop)
            dist_membr_pix = np.copy(dist_membr_pix_loop)
        else:
            dist_membr_RNA = np.append(dist_membr_RNA,dist_membr_RNA_loop,axis=0)
            dist_membr_pix = np.append(dist_membr_pix,dist_membr_pix_loop,axis=0)
        idx+=1

        # Calculate histograms
        histRNA, bins = np.histogram(dist_membr_RNA_loop,binsHist ,density=False)
        histpix, bins = np.histogram(dist_membr_pix_loop,binsHist ,density=False)

        histRNAnorm = histRNA/histRNA.sum()
        histpixnorm = histpix/histpix.sum()

        histRNAnormPix = np.divide(histRNAnorm,histpixnorm)
        histRNAnormPix = np.nan_to_num(histRNAnormPix)

        hist_plot = {'width':width,'center':center,'bins':bins,
                     'histRNA':histRNA,'histpix':histpix,
                     'histRNAnormPix':histRNAnormPix}

        # Plot results
        name_save = os.path.join(path_save,f'Z-{Zmask}.png')
        plot_results_slice(Zmask,v,mask_all,spots_loop_XY,dist_membr,hist_plot,name_save,show_plot)

        hist_slice[f'Z{Zmask}'] = hist_plot

    # Analyze all slices
    histRNA_all, bins = np.histogram(dist_membr_RNA,binsHist ,density=False)
    histRNA_all_norm = histRNA_all/histRNA_all.sum()

    # Renormalize with pixel counts
    histpix_all, bins = np.histogram(dist_membr_pix,binsHist ,density=False)
    histpix_all_norm = histpix_all/histpix_all.sum()

    hist_RNA_all_normPix = np.divide(histRNA_all_norm,histpix_all_norm)
    hist_RNA_all_normPix = np.nan_to_num(hist_RNA_all_normPix)

    hist_plot_all = {'width':width,'center':center, 'bins':bins,
                     'histRNA_all':histRNA_all,
                     'histRNA_all_norm':histRNA_all_norm,
                     'histpix_all_norm':histpix_all_norm,
                     'hist_RNA_all_normPix':hist_RNA_all_normPix}
    name_save = os.path.join(path_save,'_DistanceEnrichmentSummary.png')
    plot_results_all(hist_plot_all,name_save)
    
    if show_plot:
        show_plot(name_save)

    # Save entire analysis results as json
    input_args.pop('show_plot', None)
    analysis_results = {'args': input_args,
                        'hist_all': hist_plot_all,
                        'hist_slice': hist_slice}

    name_json = os.path.join(path_save, 'DataAll.json')

    with open(name_json, 'w') as fp:
        json.dump(analysis_results, fp,sort_keys=True, indent=4, cls=NumpyEncoder)

    # Save histogram of pooled data as csv
    name_csv = os.path.join(path_save, '_HistogramPooled.csv')
    hist_plot_all.pop('bins', None)
    hist_plot_all.pop('width', None)
    csv_header  = ';'.join(hist_plot_all.keys())
    hist_values = np.array( list(hist_plot_all.values())).transpose()
    np.savetxt(name_csv, hist_values, delimiter=";",fmt='%f',header=csv_header,comments='')

    return analysis_results


    #return analysis_results
    #plt.savefig(os.path.join(path_save, 'CellCortexDist.png'),dpi=200)
    #plt.close()

def plot_results_all(hist_plot,name_save = None,show_plot=False):

    if not show_plot:
        plt.ioff()
    
    # Get parameters to plot histogram
    center = hist_plot['center']
    width = hist_plot['width']
    histRNA_all = hist_plot['histRNA_all']
    histRNA_all_norm = hist_plot['histRNA_all_norm']
    histpix_all_norm = hist_plot['histpix_all_norm']
    hist_RNA_all_normPix = hist_plot['hist_RNA_all_normPix']

    # Plot results
    fig1, ax = plt.subplots(2,2,num='Distance comparison')
    fig1.set_size_inches((8,8))

    ax[0][0].bar(center, histRNA_all, align='center', width=width)
    ax[0][0].set_xlabel('Distance cell cortex [pix]')
    ax[0][0].set_ylabel('# RNAs')
    ax[0][0].set_xticks(center)
    ax[0][0].set_xticklabels(center.astype(int))

    ax[0][1].bar(center, histRNA_all_norm, align='center', width=width/2)
    ax[0][1].set_xlabel('Distance cell cortex [pix]')
    ax[0][1].set_ylabel('Frequency')
    ax[0][1].set_xticks(center)
    ax[0][1].set_xticklabels(center.astype(int))

    ax[1][0].bar(center, histpix_all_norm, align='center', width=width)
    ax[1][0].set_xlabel('Distance cell cortex [pix]')
    ax[1][0].set_ylabel('Renormalized counts')
    ax[1][0].set_xticks(center)
    ax[1][0].set_xticklabels(center.astype(int))

    ax[1][1].bar(center, hist_RNA_all_normPix, align='center', width=width)
    ax[1][1].set_xlabel('Distance cell cortex [pix]')
    ax[1][1].set_ylabel('Renormalized counts')
    ax[1][1].set_xticks(center)
    ax[1][1].set_xticklabels(center.astype(int))


    ax[0][0].title.set_text('RNA-absolute counts')
    ax[0][1].title.set_text('RNA-normalized counts')
    ax[1][0].title.set_text('All pixel-renormalized')
    ax[1][1].title.set_text('RNA renormalized with pixels')

    plt.tight_layout()

    if name_save:
        plt.savefig(name_save,dpi=200)
        
    if not show_plot:    
        plt.close()

def plot_results_slice(Zmask,mask,mask_all,spots_loop_XY,dist_membr,hist_plot,name_save=None,show_plot=False):

    if not show_plot:
         plt.ioff()
    
    # Find min and max values for plotting
    pad = 10
    indMaskAx0 = np.argwhere(mask_all.sum(axis=0))
    minAx0 = indMaskAx0[0]-pad
    maxAx0 = indMaskAx0[-1]+pad

    indMaskAx1 = np.argwhere(mask_all.sum(axis=1))
    minAx1 = indMaskAx1[0]-pad
    maxAx1 = indMaskAx1[-1]+pad

    # Set distance outside of cell to 0 for better plotting
    dist_membr_plot = np.copy(dist_membr)
    dist_membr_plot[np.logical_not(mask_all)] = 0

    # Get parameters to plot histogram
    center = hist_plot['center']
    width = hist_plot['width']
    histRNA = hist_plot['histRNA']
    histpix = hist_plot['histpix']
    histRNAnormPix = hist_plot['histRNAnormPix']

    # Generate plot
    fig1, ax = plt.subplots(2,3,num='Distance to cell membrane analysis. Z={}'.format(Zmask))
    fig1.set_size_inches((13,6))

    ax[0][0].imshow(mask['image'],cmap="hot")
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

    for kROI, vROI in mask['roi'].items():
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

    if name_save:
        plt.savefig(name_save,dpi=200)
        
    if not show_plot:
        plt.close()
