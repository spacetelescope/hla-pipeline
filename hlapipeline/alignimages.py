#!/usr/bin/env python

"""
This script is a modernized replacement of tweakreg.
"""

import argparse
from astropy.io import fits
from astropy.table import Table
import pdb
from stwcs.wcsutil import HSTWCS
import sys
from utils import astrometric_utils as amutils

# Module-level dictionary contains instrument/detector-specific parameters used later on in the script.
detector_specific_params = {"acs":
                                {"hrc":
                                     {"fwhmpsf":0.073, # TODO: Verify value
                                      "platescale":0.025}, # TODO: hla DGSL pipleine value: 0.025. Instrument handbook value: 0.027, and later ~0.028 × 0.025"/pixel. Verify correct value.
                                 "sbc":
                                     {"fwhmpsf":0.065, # TODO: Verify value
                                      "platescale":0.03}, # TODO: hla DGSL pipleine value: 0.030. Instrument handbook value: 0.032, and later ~0.034 × 0.030"/pixel Verify correct value.
                                 "wfc":{"fwhmpsf":0.076, # TODO: Verify value
                                      "platescale":0.05}}, # TODO: nothing. pipeline value in agreement with instrument handbook value.
                            "wfc3":
                                {"ir":{"fwhmpsf":0.14, # TODO: nothing. pipeline value in agreement with instrument handbook value.
                                      "platescale":0.09}, # TODO: hla DGSL pipleine value: 0.09. Instrument handbook value: 0.13, also later 0.135×0.121. Verify correct value.
                                 "uvis":{"fwhmpsf":0.076, # TODO: nothing. pipeline value in agreement with instrument handbook value.
                                      "platescale":0.04}}} # TODO: hla DGSL pipleine value: 0.04. Instrument handbook value: 0.04, but also later 0.0395×0.0395. Verify correct value.

#-----------------------------------------------------------------------------------------------------------------------
def main(imgList):
    """
    Main calling function.

    :param imgList: list of images
    :type imgList: list
    :return: nothing for now.
    """

    # 1: Interpret input data and optional parameters

    # 2: Apply filter to input observations to insure that they meet minimum criteria for being able to be aligned

    # 3: Build WCS for full set of input observations

    # 4: Retrieve list of astrometric sources from database

    # 5: Extract catalog of observable sources from each input image

    # 6: Cross-match source catalog with astrometric reference source catalog

    # 7: Perform fit between source catalog and reference catalog

def generate_source_catalogs(imgList):
    """
    Genreates a dictionary of source catalogs keyed by image name.

    :param imgList:
    :return:
    """
    sourceCatalogDict = {}

    #Generate composite WCS based on input images
    refwcs = amutils.build_reference_wcs(imgList)

    for imgName in imgList:
        print("Image name:           ", imgName)

        sourceCatalogDict[imgName] = {}

        # open image
        imgHDU = fits.open(imgName)
        imgPrimaryHeader = imgHDU[0].header
        instrument = imgPrimaryHeader['INSTRUME'].lower()
        detector = imgPrimaryHeader['DETECTOR'].lower()

        # get instrument/detector-specific image alignment parameters
        sourceCatalogDict[imgName]["params"] = detector_specific_params[instrument][detector]

        # Identify sources in image, convert coords from chip x, y form to reference WCS sky RA, Dec form.
        fwhmpsf_pix = sourceCatalogDict[imgName]["params"]['fwhmpsf']/sourceCatalogDict[imgName]["params"]['platescale']
        sourceCatalogDict[imgName]["catalog_table"] = amutils.generate_source_catalog(imgHDU,refwcs,threshold = 1000,fwhm = fwhmpsf_pix)

        # write out coord lists to files for diagnostic purposes. Protip: To display the sources in these files in DS9,
        # set the "Coordinate System" option to "Physical" when loading the region file.
        regfilename = imgName[0:9]+".reg"
        out_table = Table(sourceCatalogDict[imgName]["catalog_table"])
        out_table.write(regfilename, include_names=["xcentroid", "ycentroid"], format="ascii.basic")
        print("Wrote region file ",regfilename)

        print()
        imgHDU.close()
    return()
#=======================================================================================================================
if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description='Align images')
    PARSER.add_argument('inputList', nargs='+', help='A space-separated list of fits files to align, or a simple text '
                                                     'file containing a list of fits files to align, one per line')
    ARGS = PARSER.parse_args()

    # Build list of input images
    imgList=[]
    for item in ARGS.inputList:
        if item.endswith(".fits"):
            if item.endswith("asn.fits"):
                sys.exit("ADD SUPPORT FOR ASN FILES!") #TODO: Add support for asn.fits files
            else:
                imgList.append(item)
        else:
            with open(item,'r') as infile:
                fileLines = infile.readlines()
            for fileLine in fileLines:
                imgList.append(fileLine.strip())

    # Get to it!
    main(imgList)
