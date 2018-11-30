# Measure RNA distance distribution to different membranes

Goal is to quantify the distance distribution of RNAs to different membranes in
the cell. This analysis requires

* [**FISH-quant**](https://bitbucket.org/muellerflorian/fish_quant/) to detect RNA positions.
* **ImJoy** with adequate adequate plugin.

## Installing the ImJoy plugin

To install the analysis plugin for a specific analysis, select the
corresponding links below. For some tasks, you might have to install several plugins.

* **Cell membrane enrichment**. Two plugins are required.
    - <a href="https://imjoy.io/#/app?&plugin=https://raw.githubusercontent.com/muellerflorian/rna-loc/master/imjoy-plugins/MembraneDistance.imjoy.html&tag=CellMemb&w=CellMemb"  target="_blank">Processing plugin</a>
    - <a href="https://imjoy.io/#/app?&plugin=https://raw.githubusercontent.com/muellerflorian/rna-loc/master/imjoy-plugins/MembraneDistProgress.imjoy.html&w=CellMemb"  target="_blank">Progress window plugin</a>

This will open a dialog box where you can install the plugin by pressing the
`install` button.

<img src="https://raw.githubusercontent.com/muellerflorian/rna-loc/master/docs/img/imjoy-install-membdist.png" width="600px"></img>

For more information for how to install and use the pluging engine, please
consult the [ImJoy documentation](https://imjoy.io/docs/#/user-manual?id=python-engine).

## Test date

You can download already processed test data for the Cell membrane enrichment plugin, from [**Dropbox**](https://www.dropbox.com/s/0sbsmbg5xlccamp/img1.zip?dl=0). The zip archive contains data following the naming conventions of the examples below.

## Analysis overview

For each RNA, we determine the closest distance of an RNA to a membrane. One thing to keep in mind more pixel close to the membrane than far away, e.g. in the center of
the cell. A simple example is a circle. The maximum distance that you can be away from the “membrane” is the radius of the circle. However, there is only one possibility to be that far away (in the center). However, there are many more “close” positions.  Plotting a histogram of the distance to the membrane for all possible positions in the circle,  will yield a distribtion that strongly enriched for small distances.

To normalized for this effect, we calculate all possible distance from the membrane
for a given cell with a **distance transformation**. This tranformation results in
an image, where the pixel values are not fluorescence intensities but distance values. An example is shown below. The blue lines are the cell outlines. The green dots are the detected RNAs. The image is the distance transform. The intensity values are distance to the membrane in pixels.

<img src="https://raw.githubusercontent.com/muellerflorian/rna-loc/master/docs/img/dist_transform.png" width="400px"></img>

We measure for all RNAs the distance to the membrane and calculate a
histgram. We report the

-   Raw RNA distance histogram.
-   Normalized RNA histogram (values add to 1).
-   Normalized distance histogram of all pixels in the cell.
-   Normalized RNA histogram with the pixel histogram.

<img src="https://raw.githubusercontent.com/muellerflorian/rna-loc/master/docs/img/memb_summaryPlot.png" width="600px"></img>
