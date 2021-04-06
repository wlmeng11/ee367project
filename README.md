# Single-Element Ultrasound Imaging with Compressed Sensing
William Meng  
wlmeng@stanford.edu  
EE 367 Final Project  
March 19, 2021

A current version of the code can be found in the [Github repository](https://github.com/wlmeng11/ee367project).
Please see the [project page on my website](https://williammeng.com/ee367.html) for additional media.

# Installation
## MATLAB
Make sure you have a recent version of [MATLAB](https://www.mathworks.com/products/matlab.html) installed.

## Field II
This project utilizes the Field II ultrasound simulation library.
Please download and install the appropriate version of Field II from [here](https://field-ii.dk/?./downloading_8_20.html).

To install the library, you will have to add the Field II subdirectory to your MATLAB path.
You can do this by right-clicking the Field II subdirectory in MATLAB's file explorer, and then choosing *Add to Path -> Selected Folders*.

# Usage
A few versions of the script are provided:

* *NoMask.m* - Simulation and reconstruction without a delay mask. Some information about the scene can still be reconstructed due to non-uniformity in the near-field of the transducer.
* *SingleRotation.m* - Simulation and reconstruction with a delay mask, sampled with only a single rotation. The reconstructed image looks accurate in terms of features, but there is a lot of background noise.
* *MultiRotation.m* - Simulation and reconstruction with a delay mask rotated multiple times. With R=4 rotations, the reconstructed image looks very good.

The key parameters (R, electronic SNR, etc.) can be adjusted at the top of the file.
When you run each script in MATLAB, it will display some figures and save them as images in your current working directory.

# Acknowledgements
The techniques used in this project are primarily based on the following [paper](https://advances.sciencemag.org/content/3/12/e1701423):
> P. Kruizinga, et al, "Compressive 3D ultrasound imaging using a single sensor," *Science Advances*, Vol. 3, No. 12, December 2017.

Additional references are listed in my report.

Code for the ultrasound simulation was adapted from RAD 235 workshop examples.
Code for PCG and ADMM were adapted from the EE 367 homework.
