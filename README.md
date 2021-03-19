# Single-Element Ultrasound Imaging with Compressed Sensing
William Meng  
EE 367 Final Project  
March 19, 2021

# Installation
## MATLAB
Make sure you have a recent version of [MATLAB](https://www.mathworks.com/products/matlab.html) installed.

## Field II
This project utilizes the Field II ultrasound simulation library.
The Linux version is bundled in this repository for convenience.
If you are using another operating system, please download and install the appropriate version of Field II from [here](https://field-ii.dk/?./downloading_8_20.html).

To install the library, you will have to add the Field II subdirectory to your MATLAB path.
You can do this by right-clicking the Field II subdirectory in MATLAB's file explorer, and then choosing *Add to Path -> Selected Folders*.

# Usage
A few versions of the script are provided:

* NoMask.m
* SingleRotation.m
* MultiRotation.m

The key parameters can be adjusted at the top of the file.
When you run each script, it will display some figures and save them as images in your current working directory.

# Acknowledgements
The techniques used in this project are primarily based on the following [paper](https://advances.sciencemag.org/content/3/12/e1701423):
> P. Kruizinga, et al, "Compressive 3D ultrasound imaging using a single sensor," *Science Advances*, Vol. 3, No. 12, December 2017.

Additional references are listed in my report.
