""" Utility to filter datasets which cannot be aligned 

The function analyze_data opens the input FITS file to access the primary 
header data.  Based upon the values of certain FITS keywords, the function
determines whether or not this dataset can be reconciled against an astrometric 
catalog and, for multiple images, used to create a mosaic.

"""
from astropy.io import fits
from astropy.table import Table
import math

__all__ = ['analyze_data']

def analyze_data(inputFileName):
    """
    Determine if the dataset can be aligned
   
    Parameters
    ==========
    inputFileName: str
        Input FLT or FLC filename

    Returns
    =======
    outputTable: object
        Astropy Table object containing data pertaining to the dataset, including the doProcess bool

    The keyword/value pairs below define the "cannot process categories".
    OBSTYPE: not IMAGING
    MTFLAG: T
    SCAN-TYP: C or D (or !N)
    FILTER: G*, *POL*
    APERTURE: *GRISM*, G*-REF, RAMP, *POL*
    TARGNAME: DARK, TUNGSTEN, BIAS, FLAT, EARTH-CALIB, DEUTERIUM
    EXPOTIME: 0 
    """
    OBSKEY = 'OBSTYPE'
    MTKEY  = 'MTFLAG'
    SCNKEY = 'SCAN_TYP'
    FILKEY = 'FILTER'
    APKEY  = 'APERTURE'
    TARKEY = 'TARGNAME'
    EXPKEY = 'EXPTIME'

    header_hdu  = 0
    header_data = fits.getheader(inputFileName, header_hdu)

    obstype  = (header_data[OBSKEY]).upper()
    mtflag   = (header_data[MTKEY]).upper()
    scan_typ = (header_data[SCNKEY]).upper()

    sfilter  = (header_data[FILKEY]).upper()
    aperture = (header_data[APKEY]).upper()
    targname = (header_data[TARKEY]).upper()
    expotime = header_data[EXPKEY]

    # Keywords to use potentially for analysis
    instrume = (header_data['INSTRUME']).upper()
    detector = (header_data['DETECTOR']).upper()
    subarray = header_data['SUBARRAY']

    # Determine if the image has one of these conditions.  The routine
    # will exit processing upon the first satisfied condition.

    noProcKey   = None
    noProcValue = None
    doProcess = True
    # Imaging vs spectroscopic or coronagraphic
    if obstype != 'IMAGING':
        noProcKey   = OBSKEY
        noProcValue = obstype 

    # Moving target
    elif mtflag == 'T':
        noProcKey   = MTKEY
        noProcValue = mtflag 

    # Bostrophidon without or with dwell
    elif any ([scan_typ == 'C', scan_typ == 'D']):
        noProcKey   = SCNKEY
        noProcValue = scan_typ

    # Filter which begins with !F (e.g., 'POL*' or 'G*')
    elif sfilter[0] != 'F': 
        noProcKey   = FILKEY
        noProcValue = sfilter

    # Ramp, polarizer, or grism 
    elif any (['RAMP' in aperture, 'POL' in aperture, 'GRISM' in aperture, 
              '-REF' in aperture]):
        noProcKey   = APKEY
        noProcValue = aperture 

    # Calibration target
    elif any (['DARK' in targname, 'TUNG' in targname, 'BIAS' in targname, 
              'FLAT' in targname, 'DEUT' in targname, 'EARTH-CAL' in targname]):
        noProcKey   = TARKEY
        noProcValue = targname

    # Exposure time of zero
    elif math.isclose(expotime, 0.0, abs_tol=1e-5)
        noProcKey   = EXPKEY
        noProcValue = expotime 

    if (noProcKey is not None):
        doProcess = False

        # Issue message to log file 
        issue_msg(inputFileName, noProcKey, noProcValue)

    # Create an astropy table and start to populate it 
    catalog = None
    foundSources = None
    rms_x = None
    rms_y = None
    matchSources = None
    isSuccess = False
    namesArray=('instrument', 'detector', 'filter', 'aperture', 'obstype', 
            'subarray', 'doProcess', 'catalog', '# found sources', 'rms_x', 
            'rms_y', '# match sources', 'isSuccess')
    outputTable = (namesArray, [instrume, detector, sfilter, aperture, obstype, subarray, doProcess, catalog, 
                   foundSources, rms_x, rms_y, matchSources, isSuccess])

    return(outputTable)


def issue_msg(filename, key, value):
    """ Generate a message for the output log indicating the file/association will not
        be processed as the characteristics of the data are known to be inconsistent
        with alignment.
    """

    print('Dataset ' + filename + ' has keyword = value of ' + key + ' = ' + str(value) + '.\n')
    print('Dataset cannot be aligned.\n')

