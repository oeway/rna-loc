# Measure RNA distance distribution to cell membrane
Workflow to quantify the distance distribution of RNAs to the cell membrane.

* [**FISH-quant**](https://bitbucket.org/muellerflorian/fish_quant/) to detect RNA positions.
* **ImJoy**: you can install the plugin from  <a href="http://imjoy.io/#/app?w=rna-loc&plugin=muellerflorian/rna-loc:CellMembraneDistance@Stable"  target="_blank">**here**</a>

ImJoy plugins will be available in the  workspace: **`ExpGradient`**

You also need the, please consult the [ImJoy documentation](https://imjoy.io/docs/#/user-manual?id=python-engine).

## Test date
You can download already processed test data for the Cell membrane enrichment plugin,
from [**Dropbox**](https://www.dropbox.com/s/0sbsmbg5xlccamp/img1.zip?dl=0).
The zip archive contains data following the naming conventions of the examples below.

## Analysis overview

For each RNA, we determine the closest distance of an RNA to a membrane. One thing to keep in mind more pixel close to the membrane than far away, e.g. in the centre of
the cell. A simple example is a circle. The maximum distance that you can be away from the “membrane” is the radius of the circle. However, there is only one possibility to be that far away (in the centre). However, there are many more “close” positions.  Plotting a histogram of the distance to the membrane for all possible positions in the circle,  will yield a distribution that strongly enriched for small distances.

To normalised for this effect, we calculate all possible distance from the membrane
for a given cell with a **distance transformation**. This transformation results in
an image, where the pixel values are not fluorescence intensities but distance values. An example is shown below. The blue lines are the cell outlines. The green dots are the detected RNAs. The image is the distance transform. The intensity values are distance to the membrane in pixels.

<img src="https://raw.githubusercontent.com/muellerflorian/rna-loc/master/docs/img/dist_transform.png" width="400px"></img>

We measure for all RNAs the distance to the membrane and calculate a
histogram. We report the

-   Raw RNA distance histogram.
-   Normalised RNA histogram (values add to 1).
-   Normalised distance histogram of all pixels in the cell.
-   Normalised RNA histogram with the pixel histogram.

<img src="https://raw.githubusercontent.com/muellerflorian/rna-loc/master/docs/img/memb_summaryPlot.png" width="600px"></img>
