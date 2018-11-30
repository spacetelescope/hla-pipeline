#!/usr/bin/env python

"""
This script is a modernized replacement of tweakreg.
"""

import argparse
from astropy.io import fits
from utils import astrometric_utils as amutils
import pdb
from stwcs.wcsutil import HSTWCS
import sys

#-----------------------------------------------------------------------------------------------------------------------
def return_hardware_specific_parameters(instrument,detector):
    """
    Returns a dictionary containing parameters used in the image alignment process that are specific to each
    instrument/detector.

    .. note::
        Summary of parameter units:
        - fwhmpsf: arcseconds
        - platescale: arcseconds per pixel

    :param instrument: instrument name.
    :param detector: detector name
    :type instrument: string
    :type detector: string
    :return: dictionary of parameters.
    """
    instrument = instrument.lower()
    detector = detector.lower()
    print("INSTRUMENT: ",instrument)
    print("DETECTOR:   ",detector)
    paramDict = {}
    if instrument == 'acs':
        if detector == 'hrc':
            paramDict["fwhmpsf"] = 0.073 # TODO: Verify value
            paramDict['platescale'] = 0.025 # TODO: hla DGSL pipleine value: 0.025. Instrument handbook value: 0.027, and later ~0.028 × 0.025"/pixel. Verify correct value.
        elif detector == 'sbc':
            paramDict["fwhmpsf"] = 0.065 # TODO: Verify value
            paramDict['platescale'] = 0.03 # TODO: hla DGSL pipleine value: 0.030. Instrument handbook value: 0.032, and later ~0.034 × 0.030"/pixel Verify correct value.
        elif detector == 'wfc':
            paramDict["fwhmpsf"] = 0.076 # TODO: Verify value
            paramDict['platescale'] = 0.05 # TODO: nothing. pipeline value in agreement with instrument handbook value.
        else:
            sys.exit("ERROR! '%s/%s' is an unrecognized detector!"%(instrument,detector))
    elif instrument == 'wfc3':
        if detector == 'ir':
            paramDict["fwhmpsf"] = 0.14 # TODO: nothing. pipeline value in agreement with instrument handbook value.
            paramDict['platescale'] = 0.09 # TODO: hla DGSL pipleine value: 0.09. Instrument handbook value: 0.13, also later 0.135×0.121. Verify correct value.
        elif detector == 'uvis':
            paramDict["fwhmpsf"] = 0.076 # TODO: nothing. pipeline value in agreement with instrument handbook value.
            paramDict['platescale'] = 0.04  # TODO: hla DGSL pipleine value: 0.04. Instrument handbook value: 0.04, but also later 0.0395×0.0395. Verify correct value.
        else:
            sys.exit("ERROR! '%s/%s' is an unrecognized detector!" % (instrument, detector))
    else:
        sys.exit("ERROR! '%s' is an unrecognized instrument!" % (instrument))

    return(paramDict)
#-----------------------------------------------------------------------------------------------------------------------
def main(imgList, refImage):
    """
    Main calling function.

    :param imgList: list of images
    :param refImage: name of reference image.
    :type imgList: list
    :type refImage: string
    :return: nothing for now.
    """
    sourceCatalogDict = {}

    #get reference image WCS information
    refImgHDU = fits.open(refImage)
    refwcs = HSTWCS(refImgHDU, 1)

    for imgName in imgList:
        imgHDU = fits.open(imgName)
        imgPrimaryHeader = imgHDU[0].header
        # get instrument/detector-specific image alignment parameters
        detectorSpecificParams = return_hardware_specific_parameters(imgPrimaryHeader['INSTRUME'],imgPrimaryHeader['DETECTOR'])
        # Identify sources in image, convert coords from chip x, y form to reference WCS sky RA, Dec form.
        fwhmpsf_pix = detectorSpecificParams['fwhmpsf']/detectorSpecificParams['platescale']
        sourceCatalogDict[imgName] = amutils.generate_source_catalog(imgHDU,refwcs,threshold = 1000,fwhm = fwhmpsf_pix)

    pdb.set_trace()
    return()
#=======================================================================================================================
if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description='Align images')
    PARSER.add_argument('inputList', nargs='+', help='A space-separated list of fits files to align, or a simple text '
                                                     'file containing a list of fits files to align, one per line')
    PARSER.add_argument('-r', '--refImage', required=False,default=None,help='Name of the image that will be used as the '
                                                                             'reference image for the alignment. If not '
                                                                             'explicitly specified, the default is the  '
                                                                             'first image in the input list.')
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

    # Set reference image
    if ARGS.refImage:
        refImage = ARGS.refImage
    else:
        refImage = imgList[0]

    # Get to it!
    main(imgList, refImage)
