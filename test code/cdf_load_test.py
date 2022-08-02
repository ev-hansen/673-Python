import os
import platform

# dotenv and environment variable stuff
from dotenv import load_dotenv
load_dotenv()
# set up .env file to specify the home directory
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

from spacepy import pycdf

test_cdf_path = (f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/moments/rel04/"
                 "2015/rbspa_rel04_ect-hope-mom-l3_20151126_v7.1.0.cdf")
print(test_cdf_path)


def print_cdf_file(given_cdf_path=test_cdf_path):
    cdf_file = pycdf.CDF(given_cdf_path)
    print(type(cdf_file))
    print(cdf_file)
    cdf_file.close()


def main():
    print_cdf_file()


if (__name__ == "__main__"):
    main()
