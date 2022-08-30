'''
Author: Laurent Guerard
Group: IMCF
Email: laurent.guerard@unibas.ch
Creation Date: Wednesday, 16th October 2019 3:35:12 pm
-----
Last Modified: Wednesday, 13th November 2019 16:36:52
Modified By: Laurent Guerard
-----
HISTORY:
Date            By      Comments
----------      --      ---------------------------------------------------------
2019-11-13		LG 		Fixed the density count.
2019-10-23		LG 		Fixed ROI numbering. Works only on ImageJ 1.52r55 for now
2019-10-22		LG 		Added subpixel resolution. Fixed C3.
2019-10-21		LG 		Spots positions will now be based on the whole image. However, 
                        Z and C are not working at the moment.
2019-10-21      LG      Fixed some broken imports and only makedirs if not existing.
2019-10-18      LG      Finished 1st version of the script
'''

# ─── SCRIPT PARAMETERS ──────────────────────────────────────────────────────────

#@ File(label="Select the directory with your cropped images", style="directory") src_dir
#@ String(label="Extension for the images to look for") filename_filter

# ─── IMPORTS ────────────────────────────────────────────────────────────────────

from fiji.plugin.trackmate.detection import LogDetector
from net.imglib2.img.display.imagej import ImageJFunctions

import os
import csv
import glob
from itertools import izip

from ij import IJ, ImagePlus, ImageStack, WindowManager as wm
from ij.plugin.frame import RoiManager
from ij.gui import PointRoi, WaitForUserDialog
from ij.measure import ResultsTable
from ij.process import ImageConverter
from ij.plugin import Duplicator, ImageCalculator


# Bioformats imports
from loci.plugins import BF
from loci.plugins.in import ImporterOptions


# ─── VARIABLES ──────────────────────────────────────────────────────────────────

# ############################# #
# TRACKMATE DETECTION VARIABLES #
# ############################# #

# Radius for the spots /!\ NOT DIAMETER
radius_C2     = 0.425
# Threshold
threshold_C2  = 100
# Radius for the spots /!\ NOT DIAMETER
radius_C3     = 0.45
# Threshold
threshold_C3  = 90
# Radius for the spots /!\ NOT DIAMETER
radius_C4     = 0.425
# Threshold
threshold_C4  = 25
# Do Subpixel detection
doSubpixel = False
# Apply median before detection
doMedian   = False


# ─── FUNCTIONS ──────────────────────────────────────────────────────────────────


def checkForFiles(filepath):
    """Check if files are there no matter the extension

    Arguments:
        filepath {string} -- Path and name to check if exists

    Returns:
        bool -- Returns true if exists otherwise returns false
    """
    for filepath_object in glob.glob(filepath):
        if os.path.isfile(filepath_object):
            return True

    return False

def getFileList(directory, filteringString):
    """
    Returns a list containing the file paths in the specified directory
    path. The list is recursive (includes subdirectories) and will only
    include files whose filename contains the specified string.
    """
    files = []
    for (dirpath, dirnames, filenames) in os.walk(directory):
        # if out_dir in dirnames: # Ignore destination directory
            # dirnames.remove(OUT_SUBDIR)
        for f in filenames:
            if filteringString in f:
                files.append(os.path.join(dirpath, f))
    return files

def BFImport(indivFile):
    """
    Import the files using BioFormats.
    """
    options = ImporterOptions()
    options.setId(str(indivFile))
    options.setColorMode(ImporterOptions.COLOR_MODE_COMPOSITE)
    imps = BF.openImagePlus(options)
    return imps


def count_cellDetection3D(implus, current_channel, rad, thresh, subpix, med, bbox, save_file):
    """Function to detect the cells in 3D using TrackMate

    Arguments:
        implus {imagePlus}    -- ImagePlus of the image to use for detection
        current_channel {int} -- Current channel for creating the ROI
        rad    {int}          -- Radius of the cell to detect, half the diameter
        thresh {int}          -- Intensity threshold for the detection
        subpix {bool}         -- Option for subpixel detection
        med {bool}            -- Option for median filter before detection
        bbox {Rectangle}      -- Rectangle corresponding to the bounding box
        save_file {str}       -- Path to the output file containing the ROIs

    Returns:
        cellCount {int} -- Number of cells found
    """
    dim = implus.getDimensions()
    cal = implus.getCalibration()

    implus2 = implus.duplicate()
    implus2.setCalibration(cal)

    # Set the parameters for LogDetector
    img           = ImageJFunctions.wrap(implus2)
    interval      = img
    # cal           = implus.getCalibration()
    # calibration = [round(cal.pixelWidth,3), round(cal.pixelHeight,3), round(cal.pixelDepth,3)]
    calibration   = [cal.pixelWidth, cal.pixelHeight, cal.pixelDepth]
    # print calibration

    radius     = rad  # the radius is half the diameter
    threshold  = thresh
    doSubpixel = subpix
    doMedian   = med

    # print cal.pixelDepth

    # Setup spot detector (see http://javadoc.imagej.net/Fiji/fiji/plugin/trackmate/detection/LogDetector.html)
    #
    # public LogDetector(RandomAccessible<T> img,
    #            Interval interval,
    #            double[] calibration,
    #            double radius,
    #            double threshold,
    #            boolean doSubPixelLocalization,
    #            boolean doMedianFilter)

    detector = LogDetector(img, interval, calibration, radius,
                           threshold, doSubpixel, doMedian)

    # Start processing and display the results
    if detector.process():
        # Get the list of peaks found
        peaks = detector.getResult()
        # print str(len(peaks)), "peaks were found."

        # Add points to ROI manager
        rm = RoiManager(False)
        # rm.reset()

        # Loop through all the peak that were found
        for peak in peaks:
            # Print the current coordinates
            # print peak.getDoublePosition(0), peak.getDoublePosition(1), peak.getDoublePosition(2)
            # Add the current peak to the Roi manager
            # print ((peak.getDoublePosition(1) + bbox.y)/ cal.pixelHeight)
            # print cal.pixelWidth
            # print cal.pixelHeight
            # Adding 0.5 to have subpixel resolution
            roi = PointRoi((peak.getDoublePosition(0) / cal.pixelWidth) + bbox.x + 0.5,
                           (peak.getDoublePosition(1) / cal.pixelHeight) + bbox.y + 0.5)
            # roi = PointRoi((peak.getDoublePosition(0) + bbox.x),
            #                (peak.getDoublePosition(1) + bbox.y))
            # roi = PointRoi((peak.getDoublePosition(0) / cal.pixelWidth) + 0.5,
            #                (peak.getDoublePosition(1) / cal.pixelHeight) + 0.5)
            # print peak.getDoublePosition(2)/cal.pixelDepth
            # roi.setPosition(
            #     int(round(peak.getDoublePosition(2) / cal.pixelDepth)) + 1)
            roi.setPosition(
                current_channel,
                int(round(peak.getDoublePosition(2) / cal.pixelDepth)) + 1,
                1)
            # print roi
            # print roi.subPixelResolution()
            # roi.enableSubPixelResolution() 
            rm.addRoi(roi)
            # print roi
            # print roi.getZPosition()
        
        # Show all ROIs on the image
        # rm.runCommand(imp2, "Show All")
        cellCount = rm.getCount()
        # cellCount = len(peaks)
        # Close the duplicate
        implus2.changes = False
        implus2.close()
        # rm.close()
        # print rm.getCount()
        if rm.getCount() != 0:
            rm.runCommand("Save", save_file)
        rm.close()
    else:
        print "The detector could not process the data."
    return cellCount

    # print(type(markerimage))

# ─── MAIN CODE ──────────────────────────────────────────────────────────────────

# Retrieve list of files
src_dir = str(src_dir)
files = getFileList(src_dir, filename_filter)
IJ.setBackgroundColor(0, 0, 0)

# If the list of files is not empty
if files:
    # Lists for the result file
    name_list = []
    ch2_count = []
    ch3_count = []
    ch4_count = []
    ch2_density = []
    ch3_density = []
    ch4_density = []
    
    # For each file finishing with the filtered string
    for file in files:
        # Get info for the files
        IJ.log("\\Clear")
        folder   = os.path.dirname(file)
        basename = os.path.basename(file)
        basename = os.path.splitext(basename)[0]
        roi_zip  = os.path.join(folder, basename + ".zip")
        # print roi_zip

        if not os.path.exists(roi_zip):
            IJ.log("Couldn't find the ROIs for image " + basename + ", will skip it.")
            continue

        # Import the file with BioFormats
        
        IJ.log("Currently opening " + basename + "...")
        imps = BFImport(str(file))

        for imp in imps:

            # imp.show()
            # Add points to ROI manager
                # Add points to ROI manager
            rm_image = RoiManager(False)
            # rm_image.reset()

            rm_image.runCommand("Open", roi_zip)
            rois_image = rm_image.getRoisAsArray()
            if rm_image.getCount() == 0:
                IJ.log("Couldn't load the ROIs for this image. Check what happened.")
                continue

            zip_folder = os.path.join(folder,"zip_folder")

            roi_area_list = []
            # rm.close()

            for roi_index in range(rm_image.getCount()):
                out_ROI_folder = os.path.join(folder,basename,"ROI" + str(roi_index+1))
                if not os.path.exists(out_ROI_folder):
                    os.makedirs(out_ROI_folder)
                IJ.log("Working on ROI " + str(roi_index))
                # imp.setRoi(roi)
                rm_image.select(imp, roi_index)
                roi_area = imp.getStatistics().area
                # print roi_area
                roi_area_list.append(roi_area)
                current_roi = rm_image.getRoi(roi_index)
                bounding_box = current_roi.getBounds()
                # print bounding_box
                # sys.exit(0)

                # Calculate the number of cells in channel 2
                IJ.log("Looking into Channel 2")
                channel_of_interest = 2
                imp_for_tm2 = Duplicator().run(imp, channel_of_interest,
                                            channel_of_interest, 1, imp.getNSlices(), 1, 1)
                imp_for_tm2.setCalibration(imp.getCalibration())
                # Clear outside the ROI 
                rm_image.select(imp_for_tm2, roi_index)
                
                # imp_for_tm2.setRoi(roi)
                # imp_for_tm2.show()

                IJ.run(imp_for_tm2, "Clear Outside", "stack")
                # imp_for_tm2.show()

                roi_C2_zip  = os.path.join(out_ROI_folder, basename + "_ROI_" + str(roi_index+1) + "_dots_C2.zip")


                # Get the marker image with the peaks of cells using TrackMate
                cell_count_ch2 = count_cellDetection3D(
                    imp_for_tm2, channel_of_interest ,radius_C2, threshold_C2, doSubpixel, doMedian, bounding_box, roi_C2_zip)
                # print rm_image.getCount()

                # rois_C2 = rm.getRoisAsArray()
                # rm_image.runCommand("Save", roi_C2_zip)

                

                # Calculate the number of cells in channel 3
                IJ.log("Looking into Channel 3")
                channel_of_interest = 3
                imp_for_tm3 = Duplicator().run(imp, channel_of_interest,
                                            channel_of_interest, 1, imp.getNSlices(), 1, 1)
                
                imp_for_tm3.setCalibration(imp.getCalibration())

                # Clear outside the ROI 
                rm_image.select(imp_for_tm3, roi_index)
                # imp_for_tm3.setRoi(roi)
                IJ.run(imp_for_tm3, "Clear Outside", "stack")

                roi_C3_zip  = os.path.join(out_ROI_folder, basename + "_ROI_" + str(roi_index+1) + "_dots_C3.zip")

                # Get the marker image with the peaks of cells using TrackMate
                cell_count_ch3 = count_cellDetection3D(
                    imp_for_tm3, channel_of_interest ,radius_C3, threshold_C3, doSubpixel, doMedian, bounding_box, roi_C3_zip)

                # rois_C2 = rm.getRoisAsArray()
                # rm.runCommand("Save", roi_C3_zip)

                # Calculate the number of cells in channel 4
                IJ.log("Looking into Channel 4")
                channel_of_interest = 4
                imp_for_tm4 = Duplicator().run(imp, channel_of_interest,
                                            channel_of_interest, 1, imp.getNSlices(), 1, 1)

                imp_for_tm4.setCalibration(imp.getCalibration())

                imp_for_bgd = Duplicator().run(imp, channel_of_interest,
                                            channel_of_interest, 1, imp.getNSlices(), 1, 1)
                imp_for_bgd.setCalibration(imp.getCalibration())

                # Clear outside the ROI
                # rm.select(imp_for_tm4, roi_index)
                # IJ.run(imp_for_tm4, "Clear Outside", "")
                # rm.select(imp_for_bgd, roi_index)
                # IJ.run(imp_for_bgd, "Clear Outside", "")

                # Background subtraction
                ic = ImageCalculator()
                IJ.run(imp_for_bgd, "Gaussian Blur...", "sigma=20 stack")
                imp_minus_bgd = ic.run("Subtract create stack", imp_for_tm4, imp_for_bgd)
                imp_minus_bgd.setCalibration(imp.getCalibration())

                IJ.run(imp_minus_bgd, "Median 3D...", "x=2 y=2 z=2")
                # Clear outside the ROI
                rm_image.select(imp_minus_bgd, roi_index)
                # imp_minus_bgd.setRoi(roi)
                IJ.run(imp_minus_bgd, "Clear Outside", "stack")                    

                roi_C4_zip  = os.path.join(out_ROI_folder, basename + "_ROI_" + str(roi_index+1) + "_dots_C4.zip")

                # imp_minus_bgd.show()

                # Get the marker image with the peaks of cells using TrackMate
                cell_count_ch4 = count_cellDetection3D(
                    imp_minus_bgd, channel_of_interest, radius_C4, threshold_C4, doSubpixel, doMedian, bounding_box, roi_C4_zip)
                # rois_C2 = rm.getRoisAsArray()
                # rm.runCommand("Save", roi_C4_zip)

                imp.close()
                imp_for_tm2.close()
                imp_for_tm3.close()
                imp_for_tm4.close()
                imp_for_bgd.close()
                imp_minus_bgd.close()                

                name_list.append(basename)
                ch2_count.append(cell_count_ch2)
                ch3_count.append(cell_count_ch3)
                ch4_count.append(cell_count_ch4)
            
        
        outCSV = os.path.join(folder,basename,basename + "_Results.csv")

        ch2_density = [x / y for x,y in zip(ch2_count, roi_area_list)]
        ch3_density = [x / y for x,y in zip(ch3_count, roi_area_list)]
        ch4_density = [x / y for x,y in zip(ch4_count, roi_area_list)]

        # print(volList)
        with open(outCSV, 'wb') as f:
            writer = csv.writer(f)
            writer.writerow(
                ["Filename, ROI_index, ROI_area, Channel 2 count, Channel 2 density, Channel 3 count, Channel 3 density, Channel 4 count, Channel 4 density"])
            writer.writerows(izip(name_list, range(1,rm_image.getCount() + 1), roi_area_list, ch2_count, ch2_density, ch3_count, ch3_density, ch4_count, ch4_density))

        rm_image.close()

            
IJ.log('###########################')
IJ.log('Script done')
IJ.log('###########################')
