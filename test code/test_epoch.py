import os
import platform

from dotenv import load_dotenv
from datetime import datetime

load_dotenv()

# CDF Operateles
from spacepy import pycdf

import numpy

# Find OS-specific CDF_LIB path in .env file. Make sure you
# set whichever is patform-relevant to you
if (platform.system() == "Windows"):
    TEST_DIR = os.environ["WIN_TEST_DIR"]
    CDF_LIB = os.environ["WIN_CDF_LIB"]
    os.environ["CDF_LIB"] = CDF_LIB
elif (platform.system() == "Linux"):
    TEST_DIR = os.environ["LINUX_TEST_DIR"]
    CDF_LIB = os.environ["LINUX_CDF_LIB"]
    os.environ["CDF_LIB"] = CDF_LIB
elif (platform.system() == "Mac"):
    TEST_DIR = os.environ["MAC_TEST_DIR"]
    CDF_LIB = os.environ["MAC_CDF_LIB"]
    os.environ["CDF_LIB"] = CDF_LIB


cdf_path = (f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/pitchangle/rel04/2014/"
            "rbspa_rel04_ect-hope-pa-l3_20140620_v7.1.0.cdf")
test_cdf = pycdf.CDF(cdf_path)
test_dt = datetime.fromtimestamp(1403236800)


def transpose(given_list):
    """Transposes a list (swaps rows/columns)

        Parameters:
            given_list (list): A 2D list

        Returns:
            transposed_list (list): given_list transposed
    """
    arr = numpy.array(given_list)
    transposed = arr.T
    transposed_list = transposed.tolist()
    return transposed_list


def reform(given_list):
    """Acts similar to how the IDL reform functele works in the original code

        Parameters:
            given_list (list): A list, potentially w/ a dimensele of size 1

        Returns:
            reformed_list (list): Original list without dimensele of size 1
    """
    arr = numpy.array(given_list)
    reformed = arr.squeeze()
    reformed_list = reformed.tolist()
    return reformed_list


t_ion = list(test_cdf["Epoch_Ion"])
print(type(t_ion[0]))
t_ele = list(test_cdf["Epoch_Ele"])
e_data_ion = list(test_cdf["HOPE_ENERGY_Ion"])
e_data_ele = test_cdf["HOPE_ENERGY_Ele"]
mode_ion = list(test_cdf["Mode_Ion"])
mode_ele = list(test_cdf["Mode_Ele"])
fpdu_data = list(test_cdf["FPDU"])
fhedu_data = list(test_cdf["FHEDU"])
fodu_data = list(test_cdf["FODU"])
fedu_data = list(test_cdf["FEDU"])
pitchangle_data = list(test_cdf["PITCH_ANGLE"])
test_cdf.close()

t_ion_arr = numpy.array(t_ion)
t_ele_arr = numpy.array(t_ele)
e_data_ion_arr = numpy.array(e_data_ion)
e_data_ele_arr = numpy.array(e_data_ele)
mode_ion_arr = numpy.array(mode_ion)
mode_ele_arr = numpy.array(mode_ele)
fpdu_data_arr = numpy.array(fpdu_data)
fhedu_data_arr = numpy.array(fhedu_data)
fodu_data_arr = numpy.array(fodu_data)
fedu_data_arr = numpy.array(fedu_data)
pitchangle_data_arr = numpy.array(pitchangle_data)
