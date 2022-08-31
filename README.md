[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7037853.svg)](https://doi.org/10.5281/zenodo.7037853)

# 3D_spots_count

## Description

These two scripts allow for different analysis about counting spots.

* [count_3D_FISH.py](https://github.com/imcf-shareables/3D_spots_count/blob/main/count_3D_FISH.py) counts the number of spots and their density in 3D across channels in regions of interest selected by the user.
* [H_watershed_3D_nuclei.py](https://github.com/imcf-shareables/3D_spots_count/blob/main/H_watershed_3D_nuclei.py) segments nuclei in 3D and measures volume and mean intensity.

## Requirements

These scripts are made in Python and requires [Fiji](https://doi.org/10.1038/nmeth.2019) to be run. On top of this, multiple update sites need to be activated following [this guide](https://imagej.net/update-sites/#following-an-update-site): 
* 3DImageSuite
* IJPB-Plugins
* SCF-MPI-CBG

Once activated, just drag and drop the script in the main Fiji window and click on the RUN button.

As Fiji is operating system independant, this should run on any Windows, Mac and Linux. 

## Run scripts

### count_3D_FISH

#### Input

The script asks for a folder containing a list of files and a specific extension to look for. Some thresholds are fixed in the script and defined empirically based on the input datasets. 

#### Runtime 

Once initiated, the script loops through the channels of all files and uses [TrackMate](https://www.biorxiv.org/content/10.1101/2021.09.03.458852v2) and its log detector to count the spots in the selected regions. An additional step for background subtraction is applied for the last channel as the signal is less easy to identify in this one. Based on the number of spots found and the area of the ROIs, the script will measure the density of these spots and report results in a CSV.

#### Output

The script saves a CSV returning the files name that were analyzed, the ROI numbers (multiple ROIs can be applied for a single image), the area of the ROIs as well as the spots counts and densities in all channels.

### H_watershed_3D_nuclei

#### Input

The script asks for a folder containing a list of files and a specific extension to look for. A volume and an intensity thresholds for the DAPI are needed to filter the segmented objects and keep the ones of interest. A checkbox to filter objects touching the borders in Z is also there.

#### Runtime 

Once initiated, the script loops through the channels of all files. A background subtraction is applied to the DAPI channel and H-watershed with empirically selected settings are used to segment the nucleis. These are then dilated with a radius of 2 in X, Y and Z, and filtered according to the previously selected thresholds and the ones left are then saved using the [3D ROI Manager](https://academic.oup.com/bioinformatics/article/29/14/1840/231770).

#### Output

The script saves a ZIP file per image analyzed, containing the 3D nucleis which can be reopened using the 3D ROI Manager to be checked and verified.
