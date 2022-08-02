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

import pandas as pd
import random

from pyspedas.rbsp import hope

possible_datetimes = (pd.date_range(start="2012-10-25", end="2019-07-16")
                      .to_pydatetime().tolist())
random_datetime = random.choice(possible_datetimes)

probe = 'a'
yr_str = random_datetime.year  # year
m_str = random_datetime.month   # month
d_str = random_datetime.day   # day
level = "l3"
instrument = "ect"
component = "hope"
rel = "rel04"

hope(trange=[f"{yr_str}-{m_str}-{d_str}", f"{yr_str}-{m_str}-{d_str}"],
     probe=probe, downloadonly=True)  


# mirror some formatting from the IDL version
ymd_str = yr_str + m_str + d_str


# example: rbspa_rel04_ect-hope-mom-l3_20140317_v7.1.0
nm_hd = f"rbsp{probe}_{rel}_{instrument}-{component}-mom-{level}_"
fln_tmp = f"{nm_hd}{ymd_str}_"

cdf_file_parent_dir = (f"{TEST_DIR}/rbsp_data/rbsp{probe}/{level}/"
                       f"{instrument}/{component}/moments/{rel}/{yr_str}")

# Search for file in dir
file_list = os.listdir(cdf_file_parent_dir)
str_match = list(filter(lambda x: fln_tmp in x, file_list))
cdf_file_path = f"{cdf_file_parent_dir}/{str_match[0]}"

cdfObject = pycdf.CDF(cdf_file_path)
print(cdfObject)
