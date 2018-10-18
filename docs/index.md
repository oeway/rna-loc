
# RNA localization in smFISH images
This is a collection of analysis tools to study RNA localization from smFISH images.
We use a number of different open-source software packages, and detail their usage in this documentation.

[**ImJoy**](https://imjoy.io) is image processing platform with an easy to use interface powered by a Python engine running in the background.

<img src="https://raw.githubusercontent.com/muellerflorian/rna_loc/master/docs/img/imjoy-screenshot.png" width="600px"></img>

[**FISH-quant**](https://bitbucket.org/muellerflorian/fish_quant/) is a Matlab toolbox to localize RNAs in 3D from smFISH images.

<img src="https://raw.githubusercontent.com/muellerflorian/rna_loc/master/docs/img/fq-screenshot.png" width="600px"></img>

## Getting started
We provide different analysis workflows (you can select them also from the banners).
For each we specify the requirements and the required installations.

* [Analysis RNA distance distribution to cellular membranes](membraneDistance/memb-overview.md)

### Working with ImJoy
Most workflows depend on ImJoy plugins. ImJoy is an easy to use [image processing framework](https://imjoy.io/docs/#/overview).

We provide links to install the different ImJoy plugins in dedicated ImJoy workspaces. Once installed, ImJoy remembers the workspaces and plugins and you simply have to open the web app and select the appropriate workspace [https://imjoy.io/#/app](https://imjoy.io/#/app)

#### Installing plugins
If you press on the installation link, the ImJoy web app will open and display a dialog asking if you want to install the specified plugin. To confirm, press the `install` button.
<img
  src="https://raw.githubusercontent.com/muellerflorian/rna_loc/master/docs/img/install_plugin.png" width="600px"></img>

Most plugins require the **ImJoy Plugin Engine**, to perform computations in
Python. You will need to **install** it only once, but **launch** it each time
you work with ImJoy.

### ImJoy: installing and launching the Plugin Engine
We provide the Plugin Engine for different operating systems [**here**](https://github.com/oeway/ImJoy-Python/releases). Select the appropirate
one and follow the instructions.

When you **launch the Plugin Engine** for the first time it will show a webpage
with the connection token. This token (blacked out) is needed for the ImJoy app to be able to communicate with the Plugin Engine.

<img
  src="https://raw.githubusercontent.com/muellerflorian/rna_loc/master/docs/img/imjoy-connectionToken" width="600px"></img>

In the ImJoy interface press the 🚀 symbol to connect it to the Plugin Engine. This
will open an dialog where you can specify the token in the `Connection token` field. Please note, that the token will be saved and you will not need to enter it anymore to connect to this Plugin Engine.

<img
  src="https://raw.githubusercontent.com/muellerflorian/rna_loc/master/docs/img/imjoy-install-engine" width="600px"></img>
