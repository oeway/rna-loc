 # -*- coding: utf-8 -*-

"""
Introduction
============

Module containing different functions to work with FQ result files.

Usage
=====


"""
# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import math
import sys
import json
import numpy as np
import skimage
from nested_lookup import nested_lookup # pip install nested-lookuppython
from scipy import ndimage


# ---------------------------------------------------------------------------
# Globals
# ---------------------------------------------------------------------------

# Info about the module
__author__ = "Florian MUELLER"
__email__ = "muellerf.research@gmail.com"


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

# Encoder class that allows to have numpy arrays in dictionary
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def colorbar(mappable):
    """
    Function to place colorbars next to images and guarantee that they have the
    same size.

    From: https://joseph-long.com/writing/colorbars/
    More info: https://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
    """

    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)

    https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
    more info: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s\r' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()




def read_image(file_name):
    """
    Load image from file-name. Returns numpy array and dimensions
    """
    img = skimage.io.imread(file_name)

    # Move z axis last
    img1 = np.moveaxis(img, 0, 2)
    return img, img.shape


def save_image(file_name, img):
    """
    Save different images.
    """

    # Save renormalized image
    skimage.io.imsave(file_name, img)
    print("Filtering image saved as:{}".format(file_name))



def read_FQ_matlab(file_open):
    """ Opens FISH-quant result files generated with Matlab (tab-delimited text file).

    Args:
        file_open (string): string containing the full file name.

    Returns:
        dictionary containing outlines of cells, and if present the detected spots.
    """

    # Open file
    with open(file_open, "r") as fh:
         data = fh.readlines()

    # Strip white space characters
    data = [x.strip() for x in data]

     # Loop over read-in data
    fq_dict = {'cells':{},'file_names':{},'settings':{}}
    iLine = 0

    while iLine < len(data):

        line = data[iLine]


        # READ FILE NAMES
        if 'IMG_Raw' in line:
            img_name = line.split('\t')
            if len(img_name) == 2:
                fq_dict['file_names'].update({'smFISH':img_name[1]})

        if 'IMG_Filtered' in line:
            img_name = line.split('\t')
            if len(img_name) == 2:
                fq_dict['file_names'].update({'smFISH_filt':img_name[1]})

        if 'IMG_DAPI' in line:
            img_name = line.split('\t')
            if len(img_name) == 2:
                fq_dict['file_names'].update({'DAPI':img_name[1]})

        if 'FILE_settings' in line:
            img_name = line.split('\t')
            if len(img_name) == 2:
                fq_dict['file_names'].update({'settings':img_name[1]})

        # READ IMAGE PARAMETERS
        if 'PARAMETERS' in line:
            iLine += 2
            par_microscope = data[iLine].split('\t')
            fq_dict['settings'].update({'microscope':{'pix_xy':float(par_microscope[0]),
                                                      'pix_z':float(par_microscope[1]),
                                                      'RI':float(par_microscope[2]),
                                                      'EX':float(par_microscope[3]),
                                                      'EM':float(par_microscope[4]),
                                                      'NA':float(par_microscope[5]),
                                                      'type':par_microscope[6]}})

        # New cell
        if 'CELL_START' in line:

            # Get name of cell
            cell_id = line.split('\t')[1]

            ### POSITION OF CELL

            # Read X-POS
            iLine += 1
            pos_list = (data[iLine].replace('X_POS\t','')).split('\t')
            x_pos    = [int(s) for s in pos_list]

            # Read Y-POS
            iLine += 1
            pos_list = (data[iLine].replace('Y_POS\t','')).split('\t')
            y_pos    = [int(s) for s in pos_list]

            # Read Z-POS
            iLine += 1
            pos_list = (data[iLine].replace('Z_POS\t','')).split('\t')
            if len(pos_list) > 1:
                z_pos = [int(s) for s in pos_list]
            else:
                z_pos = ['']

            fq_dict['cells'].update({cell_id:{'cell_pos':{'x': x_pos,'y': y_pos,'z': z_pos}}})


        # New nucleus
        if 'Nucleus_START' in line:

            # Get name of cell
            nuc_id = line.split('\t')[1]

            ### POSITION OF CELL

            # Read X-POS
            iLine += 1
            pos_list = (data[iLine].replace('X_POS\t','')).split('\t')
            x_pos    = [int(s) for s in pos_list]

            # Read Y-POS
            iLine += 1
            pos_list = (data[iLine].replace('Y_POS\t','')).split('\t')
            y_pos    = [int(s) for s in pos_list]

            # Read Z-POS
            iLine += 1
            pos_list = (data[iLine].replace('Z_POS\t','')).split('\t')
            if len(pos_list) > 1:
                z_pos = [int(s) for s in pos_list]
            else:
                z_pos = ['']

            fq_dict['cells'][cell_id].update({nuc_id:{'nuc_pos':{'x': x_pos,'y': y_pos,'z': z_pos}}})



        # Position of detected RNAS
        if 'SPOTS_START' in line:
            iLine += 2 # Move over header

            RNA_prop = []
            while not('SPOTS_END' in data[iLine]):
                RNA_prop.append([float(s) for s in data[iLine].split('\t')])
                iLine += 1

            # Assign to dictionary
            fq_dict['cells'][cell_id].update({'spots': np.array(RNA_prop)})


        # Up date line counter
        iLine += 1

    return fq_dict



def get_rna(fq_dict):
    """
    Obtain a numpy array with all detected spots in the image. Detection results
    are saved in a dictionary (see read_FQ_results_matlab for more details).
    """
    RNAall = nested_lookup('spots', fq_dict)  # returns list of numpy arrays

    for idx,val in enumerate(RNAall):
        if idx == 0:
            spots_all = np.copy(val)
        else:
            spots_all = np.append(spots_all,val,axis=0)

    return spots_all


def calc_expression_density_plot(fq_dict,img_size,outline_int = 'max',flag_plot = False):
    """ Calculate expression density image.
    RNA detection results are used to calculate a 2D image where each cell
    is displayed with pixel values corresponding to the number of RNAs in the cell.

    Args:
        imageprop ('dict'): dictionary containing information about outlines of cells
        and nuclei as well as (if present) positions of RNA molecules

        img_size (tuple): specifying the size of the image.

        outline_int (string) specifying how pixel values of cell outlines in the
        density plot.'max' means that the maximum number of RNAs per cell is used.
        '*nt' with int being the integer value that should be used.

        flag_plot ('bool'): flag to indicate if results should be plotted.

    Returns:
        2D numpy arrays (i) image with outlines, (ii) image with
        expression density, (iii) image wiht expression density and outlines.

    """

    img_density = np.zeros(img_size, dtype=np.uint16)
    img_outline = np.zeros(img_size, dtype=np.uint8)

    # Generate image of outline and density
    iCell = 1
    print_progress(iCell, len(fq_dict['cells']))

    for key, value in fq_dict['cells'].items():
        print_progress(iCell, len(fq_dict['cells']))

        cell_pos = []
        cell_pos.append(value['cell_pos']['x'])
        cell_pos.append(value['cell_pos']['y'])
        cell_pos = np.array(cell_pos)

        # How many RNAs
        if 'spots' in value.keys():
            Nrna = value['spots'].shape[0]
        else:
            Nrna = 0

        # Create contour image
        [rr_cont, cc_cont] = skimage.drawSK.polygon_perimeter(cell_pos[1,:], cell_pos[0,:], shape=img_outline.shape, clip=True)
        img_outline[rr_cont, cc_cont] = 1

        # Generate coordinates of pixels within polygon.
        rr, cc = skimage.drawSK.polygon(cell_pos[1,:], cell_pos[0,:])
        img_density[rr, cc] = Nrna

        # Update cell counter
        iCell += 1


    ## Add outline mask to density plot

    # Decide which intensity the outline should have
    if outline_int == 'max':
        line_int = np.amax(img_density)
    else:
        line_int = float(outline_int)

    img_density_outline  = np.copy(img_density)
    img_density_outline[img_outline==1] = line_int


    # Plot results
    if flag_plot:
        fig, (ax1, ax2,ax3) = plt.subplots(3,1,num='density_plt')
        img1 = ax1.imshow(img_density,cmap="hot")
        ax1.get_xaxis().set_visible(False)
        ax1.get_yaxis().set_visible(False)
        colorbar(img1)

        img2 = ax2.imshow(img_density_outline,cmap="hot")
        ax2.get_xaxis().set_visible(False)
        ax2.get_yaxis().set_visible(False)
        colorbar(img2)

        ax3.imshow(img_outline,cmap="hot")
        ax3.get_xaxis().set_visible(False)
        ax3.get_yaxis().set_visible(False)

        plt.tight_layout(h_pad=1)
        plt.draw()

    # Return results
    return img_density,img_density_outline,img_outline



def calc_dist_enrichment(ref_pos,spots_pos,img_size,delta_dist = 100, img_density=[],flag_plot=False):
    """ Calculates the expression level as a function of the distance from a
    reference point.

    Args:
        ref_pos (tuple): position of reference point.

        spots_pos (np array): RNA positions

        img_size (tuple): size of image

        delta_dist (int): width of histogram to calculate to calculate spatial enrichment.
        Expressed in pixel.

        img_density (np array): image to be displayed when results are plotted.

        flag_plot (bool): show results be plotted.

    Returns:

        np array with histogram. 1st col: center of bins (in pixel), 2nd col
        raw counts, 3rd colum counts normalized with area of concentric circles,
        4th column cound normalized wiht number of pixels in image falling in distance range
        of each bin.
    """

    # Get distance transform image [for display purposes only]
    com_image = np.ones((img_size), dtype=np.uint8)
    com_image[int(ref_pos[0]),int(ref_pos[1])] = 0
    dist_tr_image = ndimage.distance_transform_edt(com_image)

    # Distance of all spots to reference point
    Nspots  = spots_pos.shape[0]
    RNAdist = np.sqrt(np.sum(np.square(spots_pos - np.matlib.repmat(ref_pos, Nspots,1 )),axis=1))

    RNAdist_max = np.round(np.amax(RNAdist))

    # Histogram calculation and center for display
    hist, bins = np.histogram(RNAdist, bins=np.arange(0,RNAdist_max,delta_dist),density=False)
    width = 0.8 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2


    # Renormalize counts with area [considering simple circles]
    area_bins = np.diff((list(map(lambda r: (math.pi * (r**2)),bins))))
    hist_norm_area = hist/area_bins

    # Renormalize considering how many pixels are really in the actual image
    pixels_bins = np.diff(list(map(lambda threshold: np.sum(dist_tr_image <= threshold),bins)))

    hist_norm_pixel= hist/pixels_bins

    # Summarize all histograms
    hist_all = np.stack((center,hist,hist_norm_area,hist_norm_pixel),axis=1)

    if flag_plot:

        # PLOT ROI and center of mass
        fig1, ax = plt.subplots(3,2,num='dist_enrich')
        fig1.set_size_inches((15,12))

        # Plot image with region of interest and reference point
        img1 = ax[0][0].imshow(img_density,cmap="hot")
        ax[0][0].get_xaxis().set_visible(False)
        ax[0][0].get_yaxis().set_visible(False)
        colorbar(img1)

        ax[0][0].scatter(ref_pos[1],ref_pos[0],color='g')

        # Plot distance map
        img2 = ax[1][0].imshow(dist_tr_image,cmap="hot")
        ax[1][0].get_xaxis().set_visible(False)
        ax[1][0].get_yaxis().set_visible(False)
        colorbar(img2)
        ax[1][0].scatter(ref_pos[1],ref_pos[0],color='g')

        # plot histogram of distances with various normalizations
        ax[2][0].hist(RNAdist, bins='auto')  # arguments are passed to np.histogram
        ax[2][0].set_xlabel('Distance [pixel]')

        ax[0][1].bar(center, hist, align='center', width=width)
        ax[0][1].set_xlabel('Distance [pixel]')
        ax[0][1].set_ylabel('# RNAs')

        ax[1][1].bar(center, hist_norm_area, align='center', width=width)
        ax[1][1].set_xlabel('Distance [pixel]')
        ax[1][1].set_ylabel('# RNAs/ring area')

        ax[2][1].bar(center, hist_norm_pixel, align='center', width=width)
        ax[2][1].set_xlabel('Distance [pixel]')
        ax[2][1].set_ylabel('# RNAs/area in image')

        # Set titles
        ax[0][0].title.set_text('Expression density map with reference point')
        ax[1][0].title.set_text('Distance from reference point [um]')
        ax[2][0].title.set_text('Expression as a function of distance')

        ax[0][1].title.set_text('Raw histogram with user defined range')
        ax[1][1].title.set_text('Renormalized with area')
        ax[2][1].title.set_text('Renormalized with number of pixels')

        fig1.tight_layout(h_pad=1)
        plt.draw()

    return hist_all
