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

import pyspedas.rbsp as rbsp
import cdflib

from pprint import PrettyPrinter

pp = PrettyPrinter(indent=4)

yr_str_1 = "2014"
m_str_1 = "08"
d_str_1 = "26"
date_str_1 = f"{yr_str_1}-{m_str_1}-{d_str_1}"

yr_str_2 = "2014"
m_str_2 = "08"
d_str_2 = "28"
date_str_2 = f"{yr_str_2}-{m_str_2}-{d_str_2}"

probe = 'a'

rbsp.efw(trange=[f"{date_str_1}", f"{date_str_2}"], probe=probe, 
         downloadonly=True)


test_efw_cdf_path = (f"{TEST_DIR}/rbsp_data/rbspa/l3/efw/{yr_str_1}"
                     f"/rbspa_efw-l3_{yr_str_1}{m_str_1}{d_str_1}_v04.cdf")

test_hope_cdf_path = (f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/pitchangle/"
                      f"rel04/{yr_str_1}/rbspa_rel04_ect-hope-pa-l3_"
                      f"{yr_str_1}{m_str_1}{d_str_1}_v7.1.0.cdf")

test_efw_cdf = cdflib.CDF(test_efw_cdf_path)
test_hope_cdf = cdflib.CDF(test_hope_cdf_path)

print(len(test_efw_cdf['epoch']))
print(test_efw_cdf['epoch'][0])
print(len(test_efw_cdf['spacecraft_potential']))
print(test_efw_cdf['spacecraft_potential'][0])
