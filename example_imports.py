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
import cdflib
from spacepy import pycdf
import xarray as xr
from xarray.core.dataset import Dataset as XarrDataset
from pyspedas.rbsp import hope, efw

import matplotlib as mpl
import matplotlib.pyplot as plt

from tqdm import tqdm

# Find OS-specific paths from .env file. Make sure you
# set whichever is patform-relevant to you
if (platform.system() == "Windows"):
    HOME_DIR = os.environ["WIN_HOME_DIR"]
elif (platform.system() == "Linux"):
    HOME_DIR = os.environ["LINUX_HOME_DIR"]
elif (platform.system() == "Darwin"):
    HOME_DIR = os.environ["MAC_HOME_DIR"]
