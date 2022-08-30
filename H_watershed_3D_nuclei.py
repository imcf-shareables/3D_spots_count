'''
Author: Laurent Guerard
Group: IMCF
Email: laurent.guerard@unibas.ch
Email bis: l.guerard42@gmail.com
Creation Date: Monday, 21st October 2019 4:39:24 pm
-----
Last Modified: Wednesday, 20th November 2019 18:50:35
Modified By: Laurent Guerard
-----
HISTORY:
Date              By        Comments
----------        --        ---------------------------------------------------------
2019-11-7        LG         Added batch and changed the log in Fiji
2019-10-30        LG         Added the dilation
2019-10-24        LG         Finished having the H-watershed results in 3D ROI Manager format
2019-10-23        LG         Also displaying the other channels
2019-10-21        LG         Added the display of the channel of interest for the H watershed
2019-10-18        LG         V0.1 starting the interactive H watershed to check the settings
'''


# ─── SCRIPT PARAMETERS ──────────────────────────────────────────────────────────

#@ OpService ops

#@ File(label="Select the directory with your cropped images", style="directory") src_dir
#@ String(label="Extension for the images to look for", value="tif") filename_filter
#@ Integer(label="Volume threshold", description="Discard objects with volume BELOW that threshold", value=0) min_volume
#@ Integer(label="DAPI intensity threshold", description="Discard objects with intensity value in DAPI channel BELOW that threshold", value=0) min_intensity_DAPI
#@ Boolean(label="Filter objects touching in Z", description="Discard objects touching in the first and last slice", value=False) filter_objects_touching_z

# ─── IMPORTS ────────────────────────────────────────────────────────────────────

import os

from ij import IJ, ImagePlus
from ij.plugin import Duplicator, ImageCalculator

# 3DSuite imports
from mcib3d.geom import Objects3DPopulation
from mcib3d.image3d import ImageInt, ImageHandler, Segment3DImage

# MorpholibJ imports
from inra.ijpb.morphology import Strel3D
from inra.ijpb.morphology import Morphology

# Bioformats imports
from loci.plugins import BF
from loci.plugins.in import ImporterOptions

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


# ─── MAIN CODE ──────────────────────────────────────────────────────────────────

IJ.log("\\Clear")
IJ.log("STARTING")

# ─── MAIN CODE ──────────────────────────────────────────────────────────────────

# Retrieve list of files
src_dir = str(src_dir)
files = getFileList(src_dir, filename_filter)

# If the list of files is not empty
if files:

    # For each file finishing with the filtered string
    for file in files:
        # Get info for the files
        folder   = os.path.dirname(file)
        basename = os.path.basename(file)
        basename = os.path.splitext(basename)[0]

        # Import the file with BioFormats
        IJ.log("Currently opening " + basename + "...")
        imps = BFImport(str(file))

        for imp in imps:

            # Get info about the image
            input_dir = imp.getOriginalFileInfo().directory
            filename  = os.path.splitext(imp.getOriginalFileInfo().fileName)[0]
            filename  = filename.replace(" ", "_")

            out_folder = os.path.join(input_dir,filename)
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)

            # Still open the image for testing H-watershed
            # TODO: Remove it when settings are decided
            channel_of_interest = 1
            imp_for_tm1 = Duplicator().run(imp, channel_of_interest,
                                        channel_of_interest, 1, imp.getNSlices(), 1, 1)
            channel_of_interest = 2
            imp_for_tm2 = Duplicator().run(imp, channel_of_interest,
                                        channel_of_interest, 1, imp.getNSlices(), 1, 1)


            IJ.log("    Looking into Channel 3")
            channel_of_interest = 3
            imp_for_tm3 = Duplicator().run(imp, channel_of_interest,
                                        channel_of_interest, 1, imp.getNSlices(), 1, 1)
            imp_for_bgd = Duplicator().run(imp, channel_of_interest,
                                        channel_of_interest, 1, imp.getNSlices(), 1, 1)


            # Background subtraction
            IJ.log("    Pre processing")
            ic = ImageCalculator()
            IJ.log("        Gaussian")
            IJ.run(imp_for_bgd, "Gaussian Blur...", "sigma=20 stack")
            IJ.log("        Background subtraction")
            imp_minus_bgd = ic.run("Subtract create stack", imp_for_tm3, imp_for_bgd)
            IJ.log("        Median filter")
            IJ.run(imp_minus_bgd, "Median 3D...", "x=6 y=6 z=2")
            # imp_for_tm1.show()
            # imp_for_tm2.show()
            # imp_for_tm3.show()
            # imp_minus_bgd.show()
            # IJ.run("Interactive H_Watershed")
            # IJ.selectWindow("interactive watershed-Z");

            # ─── H WATERSHED ────────────────────────────────────────────────────────────────
            IJ.log("    H-watershed")
            h_value                = 50
            segmentation_threshold = 4
            peakFlooding           = 86
            outputMask             = True # necessary so it can be used by ops connected component to create an ImgLabeling
            allowSplit             = True


            all_nuclei_mask = ops.run("H_Watershed", imp_minus_bgd, h_value, segmentation_threshold, peakFlooding, outputMask, allowSplit)

            all_nuclei_mask.setTitle("All nuclei mask")
            # all_nuclei_mask.show()

            # ─── 3D ROI MANAGER ─────────────────────────────────────────────────────────────
            IJ.log("    Dilation")
            # Segment the image for 3D Manager
            segment_3D  = Segment3DImage(all_nuclei_mask, 1, 255)
            segment_3D.segment()
            stack_label = segment_3D.getLabelledObjectsStack()
            imp_label   = ImagePlus("3D Labelled", stack_label)
            imp_label.setCalibration(imp.getCalibration())

            # Filter with Morphological opening to get rif of the rings
            # create structuring element (ball of x,y,z-radius in px)
            strel = Strel3D.Shape.BALL.fromRadiusList(2, 2, 2)

            # apply morphological opening filter to input image

            imStWTH = Morphology.dilation(imp_label.getImageStack(), strel)

            impWTH  = ImagePlus("WTH results", imStWTH)
            # assign correct calibration
            impWTH.setCalibration(imp.getCalibration())
            # impWTH.show()



            # wrap ImagePlus into 3D suite image format
            img           = ImageInt.wrap(impWTH)
            # create a population of 3D objects
            pop           = Objects3DPopulation(img)
            unit          = imp.getCalibration().getUnits()
            # print(nb)
            IH_imp_C1     = ImageHandler.wrap(imp_for_tm1)
            IH_imp_C2     = ImageHandler.wrap(imp_for_tm2)
            IH_imp_C3     = ImageHandler.wrap(imp_for_tm3)
            volList       = []
            meanIntList   = []
            feretList     = []

            obj_to_remove = []

            IJ.log("    Saving raw objects")
            raw_objects_path = os.path.join(out_folder, filename + "_raw_objects.zip")
            pop.removeObjectsTouchingBorders(img, filter_objects_touching_z)
            nb            = pop.getNbObjects()
            pop.saveObjects(raw_objects_path)

            IJ.log("    Filtering")
            # loop over the objects
            for i in range(0, nb):
                obj = pop.getObject(i)
                # if(obj.touchBorders(img, filter_objects_touching_z)):
                #     obj_to_remove.append(obj)
                #     continue
                if(obj.getVolumeUnit() < min_volume):
                    obj_to_remove.append(obj)
                    continue
                if(obj.getPixMeanValue(IH_imp_C1) < min_intensity_DAPI):
                    obj_to_remove.append(obj)
                    continue


                volList.append(obj.getVolumeUnit())
                # print(obj.getIntegratedDensity())

                # Measure volume unit
                # print(obj.getMeasure(2))

                # Measure mean intensity
                meanIntList.append(obj.getPixMeanValue(IH_imp_C2))


            IJ.log("    Saving removed objects")
            removed_objects_path = os.path.join(out_folder, filename + "_removed_objects.zip")
            pop_3D_removed = Objects3DPopulation(obj_to_remove)
            pop_3D_removed.saveObjects(removed_objects_path)

            for obj in obj_to_remove:
                pop.removeObject(obj)

            IJ.log("    Saving filtered objects")
            filtered_objects_path = os.path.join(out_folder, filename + "_filtered_objects.zip")
            pop.saveObjects(filtered_objects_path)

        IJ.log("DONE")
            # outCSV = outFullPath + ".csv"
