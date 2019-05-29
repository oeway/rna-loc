# RNA distance distribution to nuclear envelope
Workflow to quantify the distance distribution of RNAs to the cell membrane.

* [**FISH-quant**](https://bitbucket.org/muellerflorian/fish_quant/) to detect RNA positions.
* **ImJoy**: you can install the plugin from  <a href="http://imjoy.io/#/app?w=MembDist&plugin=muellerflorian/rna-loc:NuclearEnvelopeDistance@Stable"  target="_blank">**here**</a>

ImJoy plugins will be available in the  workspace: **`MembDist`**

You also need the Python plugin, please consult the [ImJoy documentation](https://imjoy.io/docs/#/user-manual?id=python-engine).


## Analysis overview

For each RNA, we determine the closest distance of an RNA to the nuclear envelope.
This distance is negative for an RNA inside a nucleus, and positive for an RNA
outside of a nucleus.

Distances for all RNAs will be summarised in histogram. These counts are then normalised
as follows
1. Normalisation for complete spatial randomness. There is "more space" away from
  a nucleus than close to it. This means that for a randomly distributed RNA
  it is less likely to be close to a nucleus, then being further away. We consider
  this by calculating a histogram with the distance of all pixels in the embryo
  to the nuclei. The RNA counts are then normalised with this pixel histogram.
2. The normalised histogram is then further normalised such that it sums up to 1.

<img src="https://raw.githubusercontent.com/muellerflorian/rna-loc/master/docs/img/nucEnvDist/dist_histogram.png" width="400px"></img>
