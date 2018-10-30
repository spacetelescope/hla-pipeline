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
        self.curdir = os.getcwd()

        for infile in input_filenames:
            self.get_input_file(infile, docopy=True)
            updatewcs.updatewcs(infile)

        output_shift_file = 'test_mosaic_ngc188_shifts.txt'
        align_to_gaia.align(input_filenames, shift_name=output_shift_file)

        shift_file = Table.read(output_shift_file, format='ascii')
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
        self.curdir = os.getcwd()

        for infile in input_filenames:
            self.get_input_file(infile, docopy=True)
            updatewcs.updatewcs(infile)

        output_shift_file = 'test_mosaic_47tuc_shifts.txt'
        align_to_gaia.align(input_filenames, shift_name=output_shift_file)

        shift_file = Table.read(output_shift_file, format='ascii')
        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)

    @pytest.mark.parametrize("input_filenames", [['j92c01b4q_flc.fits','j92c01b5q_flc.fits',
                                                  'j92c01b7q_flc.fits','j92c01b9q_flc.fits'],
                                                 ['jbqf02gzq_flc.fits', 'jbqf02h5q_flc.fits',
                                                  'jbqf02h7q_flc.fits', 'jbqf02hdq_flc.fits',
                                                  'jbqf02hjq_flc.fits', 'jbqf02hoq_flc.fits',
                                                  'jbqf02hqq_flc.fits', 'jbqf02hxq_flc.fits',
                                                  'jbqf02i3q_flc.fits', 'jbqf02i8q_flc.fits',
                                                  'jbqf02iaq_flc.fits'],
                                                 ['ib2u12kaq_flt.fits', 'ib2u12keq_flt.fits',
                                                  'ib2u12kiq_flt.fits', 'ib2u12klq_flt.fits'],
                                                 ['ibjt01a1q_flc.fits', 'ibjt01a8q_flc.fits',
                                                  'ibjt01aiq_flt.fits', 'ibjt01amq_flt.fits',
                                                  'ibjt01aqq_flt.fits', 'ibjt01auq_flt.fits',
                                                  'ibjt01yqq_flc.fits', 'ibjt01z0q_flc.fits',
                                                  'ibjt01zwq_flc.fits', 'ibjt01a4q_flc.fits',
                                                  'ibjt01acq_flc.fits', 'ibjt01akq_flt.fits',
                                                  'ibjt01aoq_flt.fits', 'ibjt01asq_flt.fits',
                                                  'ibjt01avq_flt.fits', 'ibjt01yuq_flc.fits',
                                                  'ibjt01ztq_flc.fits'],
                                                 ['ibnh02coq_flc.fits','ibnh02cmq_flc.fits',
                                                  'ibnh02c7q_flc.fits','ibnh02c5q_flc.fits',
                                                  'ibnh02cpq_flc.fits','ibnh02c9q_flc.fits',
                                                  'ibnh02bfq_flc.fits','ibnh02beq_flc.fits'],
                                                 ])
    def test_align_single_visits(self,input_filenames):
        """ Verify whether single-visit exposures can be aligned to an astrometric standard.

        Characteristics of these tests:
          * Input exposures include exposures from a number of single visit datasets to explore what impact differing
            observing modes (differing instruments, detectors, filters, subarray size, etc.) have on astrometry.

        The following datasets are used in these tests:

            * ACS dataset 10265_01: 4x F606W full-frame ACS/WFC images (assert err)
            * ACS dataset 12580_02: 5x F475W & 6x F814W ACS/WFC images (OK!)
            * WFC3 dataset 11663_12: 4x F160W full-frame WFC3/IR images (other err)
            * WFC3 dataset 12219_01: 8x F160W full-frame WFC3/IR images, 9x F336W full-frame WFC3/UVIS images (OK!)
            * WFC3 dataset 12379_02: 4X F606W, 4x F502N full-frame WFC3/UVIS images (assert err)

        """
        self.input_loc = 'base_tests'

        self.curdir = os.getcwd()
        for infile in input_filenames:
            self.get_input_file(infile, docopy=True)
            updatewcs.updatewcs(infile)

        output_shift_file = 'test_single_visits_shifts.txt'
        align_to_gaia.align(input_filenames, shift_name=output_shift_file)

        shift_file = Table.read(output_shift_file, format='ascii')
        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)
