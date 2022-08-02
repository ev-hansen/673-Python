import os
import platform

from dotenv import load_dotenv
load_dotenv()

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

# CDF Operations
from spacepy import pycdf
import numpy

test_cdf = pycdf.CDF(f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/pitchangle/"
                     "rel04/2014/rbspa_rel04_ect-hope-pa-l3_20140620_v7.1.0"
                     ".cdf")
fedu_list = list(test_cdf["FEDU"])
fpdu_list = list(test_cdf["FPDU"])
fhedu_list = list(test_cdf["FHEDU"])
fodu_list = list(test_cdf["FODU"])
arr = numpy.array(fodu_list)

fill_value = -1.0e30

print("       val             -- [i,j,k]")
for i in range(len(arr)):
    print(f"i:{i}\n")
    for j in range(len(arr[i])):
        print()        
        for k in range(len(arr[i][j])):
            data = str(arr[i][j][k]) if (arr[i][j][k] > fill_value) else str(
                float("NaN"))
            if ("e" in data): 
                print(f"{data} -- arr[{i},{j},{k}]\n{fodu_list[i][j][k]} -- "
                      f"lst[{i},{j},{k}]")
            else: 
                print(data)
