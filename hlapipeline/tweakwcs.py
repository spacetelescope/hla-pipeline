#! /bin/env python
"""
This script is a modernized replacement of tweakreg.
"""
import argparse
from astropy.io import fits
import pdb

def generate_source_catalogs(img1):
    """
    Generate source catalog for a specified input image
    :param img1: image name

    :return:
    """

def return_hardware_specific_parameters(instrument,detector):
    """
    Returns a dictionary containing parameters used in the image alignment process that are specific to each
    instrument/detector.

    :param instrument: instrument name.
    :param detector: detector name
    :type instrument: string
    :type detector: string
    :return: dictionary of parameters.
    """
    instrument = instrument.lower()
    detector = detector.lower()
    paramDict = {}
    if instrument == 'acs':
        if detector == 'hrc':
            print("ACS/HRC")
        elif detector == 'sbc':
            print("ACS/SBC")
        elif detector == 'wfc':
            print("ACS/WFC")
        else:
            sys.exit("ERROR! '%s/%s' is an unrecognized detector!"%(instrument,detector))
    if instrument == 'wfc3':
        if detector == 'ir':
            print("WFC3/IR")
        elif detector == 'uvis':
            print("WFC3/UVIS")
        else:
            sys.exit("ERROR! '%s/%s' is an unrecognized detector!" % (instrument, detector))
    else:
        sys.exit("ERROR! '%s' is an unrecognized instrument!" % (instrument))

    return(paramDict)

def main(imgList, refImage):
    """
    Main calling function.

    :param imgList: list of images
    :param refImage: name of reference image.
    :type imgList: list
    :type refImage: string
    :return: nothing for now.
    """

    # get instrument/detector-specific image alignment parameters
    # ids_paramDict = return_hardware_specific_parameters(inst,det) #TODO: flesh this out.

    for imgName in imgList:
        imgHDU = fits.open(imgName)
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