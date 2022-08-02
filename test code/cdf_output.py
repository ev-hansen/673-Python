import os
import platform

from dotenv import load_dotenv
load_dotenv()

# CDF Operations
from spacepy import pycdf

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

test_cdf = pycdf.CDF(f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/pitchangle/"
                     "rel04/2014/rbspa_rel04_ect-hope-pa-l3_20140620_v7.1.0"
                     ".cdf")
print(test_cdf)
