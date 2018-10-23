import os
import pytest
from astropy.table import Table

from stwcs import updatewcs

from base_test import BaseHLATest
from hlapipeline import align_to_gaia

@pytest.mark.bigdata
class TestAlignMosaic(BaseHLATest):
    """ Tests which validate whether mosaics can be aligned to an astrometric standard.

        Characeteristics of these tests:
          * A single astrometric catalog was obtained with both GAIA and non-GAIA
            (PanSTARRS?) sources for the entire combined field-of-view using the GSSS
            server.
              * These tests assume combined input exposures for each test have enough
                astrometric sources from the external catalog to perform a high quality
                fit.
          * This test only determines the fit between sources
            extracted from the images by Tweakreg and the source positions included in
            the astrometric catalog.
          * The WCS information for the input exposures do not get updated in this test.
          * No mosaic gets generated.

        Success Criteria:
          * Success criteria hard-coded for this test represents 10mas RMS for the
            WFC3 images based on the fit of source positions to the astrometric catalog
            source positions.
              * RMS values are extracted from optional shiftfile output from `tweakreg`
              * Number of stars used for the fit and other information is not available
                with the current version of `tweakreg`.

    """

    ref_loc = ['truth']

    def test_align_ngc188(self):
        """ Verify whether NGC188 exposures can be aligned to an astrometric standard.

        Characeteristics of this test:
          * Input exposures include both ACS and WFC3 images of the same general field-of-view
            of NGC188 suitable for creating a combined mosaic using both instruments.
        """
        self.input_loc = 'mosaic_ngc188'
        input_filenames = ['iaal01hxq_flc.fits', 'iaala3btq_flc.fits',
                            'iaal01hyq_flc.fits', 'iaala3bsq_flc.fits',
                            'j8boa1m8q_flc.fits', 'j8boa1m4q_flc.fits',
                            'j8boa1maq_flc.fits', 'j8boa1m6q_flc.fits']
        self.output_shift_file = 'test_mosaic_ngc188_shifts.txt'
        shift_file = self.run_align(input_filenames)

        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)

    def test_align_47tuc(self):
        """ Verify whether 47Tuc exposures can be aligned to an astrometric standard.

        Characeteristics of this test:
          * Input exposures include both ACS and WFC3 images of the same general field-of-view
            of 47Tuc suitable for creating a combined mosaic using both instruments.
        """
        self.input_loc = 'mosaic_47tuc'
        input_filenames = ['ib6v06c4q_flc.fits','ib6v06c7q_flc.fits',
                                'ib6v25aqq_flc.fits','ib6v25atq_flc.fits',
                                'jddh02gjq_flc.fits','jddh02glq_flc.fits',
                                'jddh02goq_flc.fits']
        self.output_shift_file = 'test_mosaic_47tuc_shifts.txt'
        shift_file = self.run_align(input_filenames)

        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)

    def test_astroquery(self):
        """Verify that new astroquery interface will work"""
        self.curdir = os.getcwd()
        self.input_loc = ''

        shift_file = self.run_align('ib6v06060')
        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)

    def test_align_randomFields(self):
        """ Process a large number of randomly selected fields which are defined in
            an input ascii file (CSV).

            Each field is used as input to determine if it can be aligned to an 
            astrometric standard.  The success or fail status for each test is retained
            as the overall success or fail statistic is the necessary output from
            this test.
        """

        # Need to accommodate a pre-defined list existing in local store or Artifactory
        # via the setting of the environment variable TEST_BIGDATA
        #input_list_file= ['ACSList50.csv']
        input_list_file = ['ACSList5.csv']

        self.input_loc = 'master_lists'
        self.curdir = os.getcwd()

        numSuccess = 0
        numAllDatasets = 0

        for input_file in input_list_file:
            # Obtain the full path to the file containing the dataset names
            input_file_path = self.get_data(input_file)

            # Read the file and extract a list of each dataset name in IPPSSOOT format
            # which is either an association ID or an individual filename
            dataset_list = get_dataset_list(input_file_path[0])

            numAllDatasets += len(dataset_list)

            # Reset the input location as the actual data will be obtained
            # from the MAST archive via astroquery
            self.input_loc = ''

            # Process the dataset names in the list
            #
            # If the dataset name represents an association ID, the multiplicity 
            # of images within the association need to be processed.  Otherwise,
            # the dataset is a single image.
            for dataset in dataset_list:
                filenames = self.get_input_file(dataset, docopy=True)
                print(filenames)
                for infile in filenames:
                    updatewcs.updatewcs(infile)
  
                try:
                    output_shift_file = 'test_mosaic_shifts.txt'
                    align_to_gaia.align(filenames, shift_name=output_shift_file)

                    shift_file = Table.read(output_shift_file, format='ascii')
                    rms_x = max(shift_file['col6'])
                    rms_y = max(shift_file['col7'])

                    if ((rms_x <= 0.25) and (rms_y <= 0.25)):
                        numSuccess += 1

                    # Clean up
                    os.remove("test_mosaic_shifts.txt")

                except ValueError:
                    pass
             
        # Determine the percent success over all datasets processed
        percentSuccess = numSuccess/numAllDatasets
        print('Number of successful tests: ', numSuccess, ' Total number of tests: ', numAllDatasets, ' Percent success: ', percentSuccess)
 
        assert(percentSuccess >= 0.70)

def get_dataset_list(filename):
    """ Standalone function to read the master file list and get the dataset names"""

    dataFromTable = Table.read(filename, format='ascii')
    datasetIDs = dataFromTable['observationID'][:10]
    asnIDs     = dataFromTable['asnID'][-10:]

    datasetNames = []

    # Determine if the data is part of an association or is an individual image
    for imgid,asnid in zip(datasetIDs,asnIDs):

        # If the asnID is the string NONE, this is an individual image,
        # and it is necessary to get the individual image dataset name.  
        # Otherwise, this is an association dataset, so just add the asnID.
        if (asnid.upper() == "NONE"):
            datasetNames.append(imgid)
        else:
            datasetNames.append(asnid)

    return(datasetNames)
