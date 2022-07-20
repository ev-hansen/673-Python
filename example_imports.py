# Standard Imports #

import datetime as dt
import argparse
import os
import platform
from typing import List, Dict, Tuple, Any


# Third-Party Inports #

from dotenv import load_dotenv
load_dotenv()

# CDF Operations
from spacepy import pycdf
from pyspedas.rbsp import hope

# Calculations and array/list ops
import numpy

from tqdm import tqdm


# Find OS-specific CDF_LIB path in .env file. Make sure you
# set whichever is patform-relevant to you
if (platform.system() == "Windows"):
    HOME_DIR = os.environ["WIN_HOME_DIR"]
    CDF_LIB = os.environ["WIN_CDF_LIB"]
    os.environ["CDF_LIB"] = CDF_LIB
elif (platform.system() == "Linux"):
    HOME_DIR = os.environ["LINUX_HOME_DIR"]
    CDF_LIB = os.environ["LINUX_CDF_LIB"]
    os.environ["CDF_LIB"] = CDF_LIB
elif (platform.system() == "Darwin"):
    HOME_DIR = os.environ["MAC_HOME_DIR"]
    CDF_LIB = os.environ["MAC_CDF_LIB"]
    os.environ["CDF_LIB"] = CDF_LIB
