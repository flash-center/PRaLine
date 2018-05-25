# PRaLine code (Proton Radiography Linear reconstruction, the "reconstruction" is silent)

Proton radiography is an experimental technique used to reveal the magnetic fields found in high-energy density plasma experiments using beams of protons. The proton beams are either produced by imploding capsules, which release a short pulse of mono-energetic protons at stagnation, or by thermally produced protons accelerated through extreme electric field gradients. The usual proton energies are between 3 and 30MeV and are driven through a magnetic field in a plasma experiment, where the path is shifted by the Lorentz force and then come to stop at the screen. The positions on the screen are recorded along with initial conditions. Under certain experimental conditions, the full structure of the perpendicular magnetic field can be reconstructed by solving a steady- state inhomogeneous two-dimensional diffusion equation. The code presented is a Python package developed at the University of Chicago to analyze proton radiography experimental radiographs in the linear (small-image contrast) regime. This code is based on work by Graziani et al.:

Graziani, C., Tzeferacos, P., Lamb, D. Q. & Li, C. Inferring morphology and strength of magnetic fields from proton radiographs. Review of Scientific Instruments 88, 123507 (2017).

The published manuscript may be found [here](https://doi.org/10.1063/1.5013029) and an earlier open-access draft may be found [here](https://arxiv.org/abs/1603.08617).

# Setup

## Dependencies
This module requires **Python 2.7** or **3.5**. Installation requires **git**.

**OS X users:** Prior to installing dependencies, ensure an adequate Python installation by following [this guide](https://matplotlib.org/faq/installing_faq.html#osx-notes). The Python that ships with OS X may not work well with some required dependencies.

The following Python packages are required:
* future (Cross-compatibility between Python2 and Python3)
* numpy (Scientific computing)
* scipy (Scientific computing)
* matplotlib (Plotting)
* pandas (Parsing Large files)
* pradreader (https://github.com/flash-center/PRadReader) (Reading various proton radiograph file formats)

On most systems (see above note for OS X), they can be installed using Python's [PIP package manager](https://packaging.python.org/tutorials/installing-packages/) as follows:

```shell
pip install future
pip install numpy scipy matplotlib pandas
pip install git+https://github.com/flash-center/PRadReader
```

Depending on how Python was installed on your system, `pip` may require *Administrative* or `sudo` privileges.

## Installation
Once all dependencies are satisfied, install the latest version of **PRaLine** by:

```shell
pip install git+https://github.com/flash-center/PRaLine
```

The module can also be installed by:

```shell
git clone https://github.com/flash-center/PRaLine
cd PRaLine
python setup.py install
```

# Usage
## Requirements
An intermediate file is created using the [PRadReader python package](https://github.com/flash-center/PRadReader).

An intermediate file that contains the variables such as:
* Distance from proton source to the interaction region(cm), s2r_cm
* Distance from proton source to the screen(cm), s2d_cm
* Proton Kinetic Energies (MeV), Ep_MeV
* flux image which is a matrix with number of protons per bin of the screen which is dependent on the inputted bin length, flux
* flux reference which is the flux image if there were no interaction region

## How it Works

The user has to run the PRadReader on there proton radiograph data which returns a intermediate file that can be used on both command line tools below.
## Command Line Tools

Supported file formats include any that PRadReader supports, including radiographs generated from FLASH simulations and MIT's CR39 proton radiography analysis
### Tool 1: "lin-reconstruct"

A command line tool for reconstructing the magnetic field of the data from a proton radiography experiment given an intermediate file with the requirements above. It outputs streamplots based on the reconstruction algorithim  
#### Usage
```shell
lin-reconstruct [options] [intermediate file] 
```
##### Options

| Option | Action |
|:-------|--------|
|--tol| The Gauss-Seidel tolerance. DEFAULT:1.0E-04 |
|--iter| The number of Gauss-Seidel iterations. DEFAULT:4000|

**The number of Gauss-Seidel iterations**: This number represents the number of iterations that in the Gauss-Seidel method. Changing this number may affect the results if it hasn't reached converegence.

**The Gauss-Seidel tolerance**: This number represents  limit  and a resudial value that is calculated every iteration. If the resudial value reaches this limit the Gauss-Seidel method should stop.

For more info check out pages 8 and 9: https://arxiv.org/abs/1603.08617

#### Example
```shell
lin-reconstuct --tol 1.0E-05 --iter 8000 input.txt
```
This command line script ensures that Gauss-Seidel Tolerance is 1.0E-05 and the number of Gauss-Seidel Iterations 8000 and parses the input.txt constructed by [PRadReader](https://github.com/flash-center/PRadReader).
#### Output
The tool outputs Log Reconstructed Perpendicular Magnetic Field Projection

### Tool 2: "lin-analyze"

A command line tool for analysis of a proton radiography experiment. Analysis is done by plotting a 2D matrices each value is considered a pixel on the graph and the number of pixels is determinded by the bin size that the user inputs. The 2D matrices are the flux and fluence(fluence distribution of protons) plot.

#### Usage
```shell
lin-analyze [intermediate file]
```
#### Example
```shell
lin-analyze input.txt
```
This command line parses input.txt that has been constructed by [PRadReader](https://github.com/flash-center/PRadReader).
#### Output
The tool outputs a flux and fluence contrast plot

## Example Problem
There is an example intermediate file, test_input.txt, in the `examples/` directory which was generated from the magnetic field configuartion in the paper using [PRadReader](https://github.com/flash-center/PRadReader). 

**Instructions**:
1. The user should follow the instructions for the installation of package above and also clone the repository 
2. Navigate into the `examples/` directory and copy the file,`test_input.txt`, out. 
3. From the same directory run this command
```shell
lin-reconstruct test_input.txt
``` 
Output (after 4000 iterations)

<img src="examples/reference_images/B_Reconstructed.png" width="425"/>

```shell
lin-analyze test_input.txt
```
Output
<p float="left">
<img src="examples/reference_images/Flux.png" width="400" height="300"/>
<img src="examples/reference_images/Fluence.png" width="400" height="300"/>
</p>

# Updating/Uninstalling
Write up how to update a current installation (and how to update dependencies as well), how to uninstall it

To update **PRaLine** at a later date
```shell
pip install --upgrade git+https://github.com/flash-center/PRadReader
pip install --upgrade git+https://github.com/flash-center/PRaLine
```

To uninstall **PRaLine**
```shell
pip uninstall praline
```
