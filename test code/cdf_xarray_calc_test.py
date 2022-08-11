# First party imports

import os
import platform

# Third party imports
import numpy
from tqdm import tqdm

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

test_cdf = cdflib.cdf_to_xarray(f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/"
                                "pitchangle/rel04/2014/rbspa_rel04_ect-"
                                "hope-pa-l3_20140620_v7.1.0.cdf")


pa = test_cdf['PITCH_ANGLE']
fodu = test_cdf['FODU']

pi = 3.141592653589793
pa_arr = pa * (pi / 180)
npa = len(pa_arr)

del_pa = [pa_arr[1] - pa_arr[0] if ii == 0 else 
          pa_arr[npa - 1] - pa_arr[npa - 2] if ii == npa - 1 else
          (pa_arr[ii + 1] - pa_arr[ii - 1]) / 2.0 
          for ii in tqdm(range(0, npa), desc='del_pa')]

print(f"npa {npa}\ndel_pa len {len(del_pa)}")

ang_perp = 0.5 * (numpy.sin(pa_arr) ** 3)
ang_para = numpy.sin(pa_arr) * (numpy.cos(pa_arr) ** 2)

# species dims - time (3k+), PA (11), Energy (72)

test_time = fodu[0]

test_perp = test_time
test_para = test_time

for i in range(len(test_time)):
    test_perp = test_perp[i] * (ang_perp * del_pa)
    test_para = test_para[i] * (ang_para * del_pa)

print(pa)
print(test_time['PITCH_ANGLE'])
