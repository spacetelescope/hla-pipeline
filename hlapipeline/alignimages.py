#!/usr/bin/env python

"""This script is a modernized replacement of tweakreg.

"""

from astropy.io import fits
from astropy.table import Table
import glob
import os
import pdb
from stwcs.wcsutil import HSTWCS
import sys
from utils import astrometric_utils as amutils
from utils import astroquery_utils as aqutils
from utils import filter


# Module-level dictionary contains instrument/detector-specific parameters used later on in the script.
detector_specific_params = {"acs":
                                {"hrc":
                                     {"fwhmpsf": 0.073,
                                      "classify": True,
                                      "threshold": None},
                                 "sbc":
                                     {"fwhmpsf": 0.065,
                                      "classify": False,
                                      "threshold": 10},
                                 "wfc":
                                     {"fwhmpsf": 0.076,
                                      "classify": True,
                                      "threshold": None}},  # TODO: Verify ACS fwhmpsf values
                            "wfc3":
                                {"ir":
                                     {"fwhmpsf": 0.14,
                                      "classify": False,
                                      "threshold": None},
                                 "uvis":
                                     {"fwhmpsf": 0.076,
                                      "classify": True,
                                      "threshold": None}}}
# ----------------------------------------------------------------------------------------------------------------------


def check_and_get_data(input_list):
    """Verify that all specified files are present. If not, retrieve them from MAST.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    input_file_list : list
        list of full filenames

    """
    totalInputList=[]
    for input_item in input_list:
        if input_item.endswith("0"): #asn table
            totalInputList.append(aqutils.retrieve_observation(input_item))

        else: #single file rootname.
            fitsfilename = glob.glob("{}_flc.fits".format(input_item))
            if not fitsfilename:
                fitsfilename = glob.glob("{}_flt.fits".format(input_item))
            fitsfilename = fitsfilename[0]

            if not os.path.exists(fitsfilename):
                imghdu = fits.open(fitsfilename)
                imgprimaryheader = imghdu[0].header
                totalInputList.append(aqutils.retrieve_observation(imgprimaryheader['ASN_ID'].strip()))
            else: totalInputList.append(fitsfilename)
    print("TOTAL INPUT LIST: ",totalInputList)
    # TODO: add trap to deal with non-existent (incorrect) rootnames
    # TODO: add "Clobber" and "Archive" options to aqutils.retrieve_observation
    return(totalInputList)




def perform_align(input_list):
    """Main calling function.

    Parameters
    ----------
    input_list : list
        List of one or more IPPSSOOTs (rootnames) to align.

    Returns
    -------
    Nothing for now.

    """

    # Define astrometric catalog list in priority order
    catalogList = ['GAIADR2', 'GSC241']
    numCatalogs = len(catalogList)

    # 1: Interpret input data and optional parameters
    imglist = check_and_get_data(input_list)

    # 2: Apply filter to input observations to insure that they meet minimum criteria for being able to be aligned
    filteredTable = filter.analyze_data(imglist)
    pdb.set_trace()

    # 3: Build WCS for full set of input observations
    refwcs = amutils.build_reference_wcs(imglist)

    # 4: Retrieve list of astrometric sources from database
    reference_catalog = generate_astrometric_catalog(imglist, catalog='GAIADR2', existing_wcs=refwcs)

    # 5: Extract catalog of observable sources from each input image
    extracted_sources = generate_source_catalogs(imglist, refwcs, threshold=1000)

    # 6: Cross-match source catalog with astrometric reference source catalog

    # 7: Perform fit between source catalog and reference catalog
# ----------------------------------------------------------------------------------------------------------------------


def generate_astrometric_catalog(imglist, **pars):
    """Generates a catalog of all sources from an existing astrometric catalog are in or near the FOVs of the images in
        the input list.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    ref_table : object
        Astropy Table object of the catalog

    """
    # generate catalog
    out_catalog = amutils.create_astrometric_catalog(imglist,**pars)

    # write catalog to ascii text file
    catalog_fileanme = "refcatalog.cat"
    out_catalog.write(catalog_fileanme, format="ascii.fast_commented_header")

    print("Wrote reference catalog {}.".format(catalog_fileanme))
    return(out_catalog)
# ----------------------------------------------------------------------------------------------------------------------


def generate_source_catalogs(imglist, refwcs, **pars):
    """Generates a dictionary of source catalogs keyed by image name.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for source detection.
    refwcs : astropy WCS object
        Composite WCS information generated by build_reference_wcs() for the full set of input images.

    Returns
    -------
    sourcecatalogdict : dictionary
        a dictionary (keyed by image name) of two element dictionaries which in tern contain 1) a dictionary of the
        detector-specific processing parameters and 2) an astropy table of position and photometry information of all
        detected sources
    """
    sourcecatalogdict = {}
    for imgname in imglist:
        print("Image name: ", imgname)

        sourcecatalogdict[imgname] = {}

        # open image
        imghdu = fits.open(imgname)
        imgprimaryheader = imghdu[0].header
        instrument = imgprimaryheader['INSTRUME'].lower()
        detector = imgprimaryheader['DETECTOR'].lower()

        # get instrument/detector-specific image alignment parameters
        if instrument in detector_specific_params.keys():
            if detector in detector_specific_params[instrument].keys():
                sourcecatalogdict[imgname]["params"] = detector_specific_params[instrument][detector]
            else:
                sys.exit("ERROR! Unrecognized detector '{}'. Exiting...".format(detector))
        else:
            sys.exit("ERROR! Unrecognized instrument '{}'. Exiting...".format(instrument))

        # Identify sources in image, convert coords from chip x, y form to reference WCS sky RA, Dec form.
        imgwcs = HSTWCS(imghdu, 1)
        fwhmpsf_pix = sourcecatalogdict[imgname]["params"]['fwhmpsf']/imgwcs.pscale
        sourcecatalogdict[imgname]["catalog_table"] = amutils.generate_source_catalog(imghdu, refwcs, fwhm=fwhmpsf_pix, **pars)

        # write out coord lists to files for diagnostic purposes. Protip: To display the sources in these files in DS9,
        # set the "Coordinate System" option to "Physical" when loading the region file.
        regfilename = imgname[0:9]+".reg"
        out_table = Table(sourcecatalogdict[imgname]["catalog_table"])
        out_table.write(regfilename, include_names=["xcentroid", "ycentroid"], format="ascii.basic")
        print("Wrote region file {}\n".format(regfilename))
        imghdu.close()
    return(sourcecatalogdict)
# ======================================================================================================================


if __name__ == '__main__':
    import argparse
    PARSER = argparse.ArgumentParser(description='Align images')
    PARSER.add_argument('raw_input_list', nargs='+', help='A space-separated list of fits files to align, or a simple text '
                                                     'file containing a list of fits files to align, one per line')
    ARGS = PARSER.parse_args()

    # Build list of input images
    input_list = []
    for item in ARGS.raw_input_list:
        if os.path.exists(item):
            with open(item, 'r') as infile:
                fileLines = infile.readlines()
            for fileLine in fileLines:
                input_list.append(fileLine.strip())
        else:
            input_list.append(item)


    # Get to it!
    perform_align(input_list)

