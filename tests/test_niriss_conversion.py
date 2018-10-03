#!/usr/bin/env python
"""Tests for the NIRISS magnitude conversion.

Authors
-------

    Johannes Sahlmann

"""

import os
import sys

from astropy.table import Table
from configobj import ConfigObj
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../'))
import jwst_magnitude_converter


# import importlib
# importlib.reload( jwst_magnitude_converter )
# import jwst_magnitude_converter

test_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data')


def test_niriss():
    """Test magnitude conversion to several NIRISS filters."""
    instrument = 'NIRISS'

    out_dir = os.path.join(test_data_dir, 'temporary_data')
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # filter_name = '{} F200W'.format(instrument)
    # niriss_filters = np.array([filter_name])
    niriss_filters = np.array(['{} {}'.format(instrument, filter) for filter in
                               ['F200W', 'F277W', 'F380M', 'F480M']])


    catalog_file = os.path.join(test_data_dir, instrument.lower(), 'catalog.txt')
    catalog = Table.read(catalog_file, format='ascii.basic', delimiter=',')


    temp_cat_file = os.path.join(out_dir, 'temp_catalog.txt')
    catalog.write(temp_cat_file, format='ascii.no_header', overwrite=True)

    # filter names in jwst_magnitude_conversion convention
    filter1_name = 'HST_ACS_F606W'
    filter2_name = '2MASS_Ks'

    # column names in input catalog
    column1_name = 'vcal'
    column2_name = 'kuse'

    # write configuration file for converter
    mag_conv_config = ConfigObj()
    mag_conv_config['Input_Magnitude_Parameters'] = {}
    mag_conv_config['Input_Magnitude_Parameters']['filter1'] = filter1_name
    mag_conv_config['Input_Magnitude_Parameters']['filter2'] = filter2_name
    mag_conv_config['Input_Magnitude_Parameters']['column1'] = catalog.colnames.index(
        column1_name) + 1
    mag_conv_config['Input_Magnitude_Parameters']['column2'] = catalog.colnames.index(
        column2_name) + 1
    mag_conv_config['Input_Magnitude_Parameters']['column1type'] = 'magnitude'
    mag_conv_config['Input_Magnitude_Parameters']['column2type'] = 'magnitude'
    mag_conv_config['Input_Magnitude_Parameters']['yvalue'] = 1
    mag_conv_config['Input_Magnitude_Parameters']['datafile'] = temp_cat_file
    mag_conv_config['Input_Magnitude_Parameters']['racolumn'] = -1
    mag_conv_config['Input_Magnitude_Parameters']['deccolumn'] = -1

    mag_conv_config['Output_Filter_Values'] = {}
    mag_conv_config['Output_Filter_Values']['modelset'] = 'Kurucz'
    mag_conv_config['Output_Filter_Values']['fitorder'] = 4
    mag_conv_config_filename = os.path.join(out_dir, 'mag_conv.cfg')
    mag_conv_config.filename = mag_conv_config_filename

    # loop through filters
    for niriss_filter in niriss_filters:
        filename = os.path.join(out_dir, 'converted_magnitudes_{}.txt'
                                .format(niriss_filter).replace(' ', '_'))
        mag_conv_config['Output_Filter_Values']['outfilename'] = filename
        # have to set niriss1 and niriss2 fields to same value to return single filter output
        # mag_conv_config['Output_Filter_Values']['niriss1'] = niriss_filter
        # mag_conv_config['Output_Filter_Values']['niriss2'] = niriss_filter
        mag_conv_config['Output_Filter_Values']['jwst1'] = niriss_filter
        mag_conv_config['Output_Filter_Values']['jwst2'] = niriss_filter
        mag_conv_config.write()

        jwst_magnitude_converter.main(['jwst_magnitude_converter.py', mag_conv_config_filename])

    for niriss_filter in niriss_filters:
        filename = os.path.join(out_dir, 'converted_magnitudes_{}.txt'
                                .format(niriss_filter).replace(' ', '_'))
        niriss_mags = Table.read(filename, format='ascii.no_header',
                                 names=('tmp1', 'tmp2', niriss_filter))
        # add column to catalog containing desired NIRISS filter magnitude
        catalog[niriss_filter] = niriss_mags[niriss_filter]

        assert len(niriss_mags[niriss_filter]) == len(catalog)


    catalog.pprint()
