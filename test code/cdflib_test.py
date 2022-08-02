# First party imports

import os
import platform

# Third party imports

# dotenv and environment variable stuff
from dotenv import load_dotenv
load_dotenv()

if (platform.system() == "Windows"):
    TEST_DIR = os.environ["WIN_TEST_DIR"]
elif (platform.system() == "Linux"):
    TEST_DIR = os.environ["LINUX_TEST_DIR"]
elif (platform.system() == "Darwin"):
    TEST_DIR = os.environ["MAC_TEST_DIR"]

import cdflib
from datetime import datetime as dt

import pprint
pp = pprint.PrettyPrinter(indent=4)

start_time = dt.strptime("2014-06-20 19:36:30", "%Y-%m-%d %H:%M:%S")
end_time = dt.strptime("2014-06-20 19:55:48", "%Y-%m-%d %H:%M:%S")


test_file_1 = (f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/pitchangle/rel04/2014/"
               "rbspa_rel04_ect-hope-pa-l3_20140826_v7.1.0.cdf")

test_file_2 = (f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/pitchangle/rel04/2014/"
               "rbspa_rel04_ect-hope-pa-l3_20140827_v7.1.0.cdf")

test_file_3 = (f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/pitchangle/rel04/2014/"
               "rbspa_rel04_ect-hope-pa-l3_20140828_v7.1.0.cdf")

test_cdf_1 = cdflib.cdf_to_xarray(test_file_1, fillval_to_nan=True, 
                                  to_unixtime=True)

test_cdf_2 = cdflib.cdf_to_xarray(test_file_2, fillval_to_nan=True, 
                                  to_unixtime=True)

test_cdf_3 = cdflib.cdf_to_xarray(test_file_3, fillval_to_nan=True, 
                                  to_unixtime=True)

# xr.concat(objs=[test_cdf_1, test_cdf_2, test_cdf_3])

low = [-1, 1]
high = [6, 7]

test_data = test_cdf_1['HOPE_ENERGY_Ion'][0]
print(test_data)
print()
print(test_data.where(low[0] < test_data.HOPE_ENERGY_Ion_dim).where(
    test_data.HOPE_ENERGY_Ion_dim < high[0]))
