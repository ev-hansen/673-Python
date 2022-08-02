import os
import platform

from datetime import datetime
from datetime import timedelta

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

from spacepy import pycdf
import pyspedas


fmax = 6.
fmin = 2.
df = fmax - fmin
iener = 1
ipa = 5
date1 = '2015-3-18/00:00'
date2 = '2015-3-18/06:00'

iyear = 2015   # year
idoy = 77      # day of the year
dt = 5         # dt in minute
ihourf = 24    # end of simulation from t=0 in hour

# read date 
date0 = datetime.strptime(str(iyear) + ' ' + str(idoy), '%Y %j')
datef = date0 + timedelta(hours=ihourf)
print(date0)
print(datef)
date1 = date0.strftime('%Y-%m-%d/%H:%M')  # Original was '%Y-%-m-%d/%H:%M' 
date2 = datef.strftime('%Y-%m-%d/%H:%M')
print(date1)
print(date2)
ndata = int(ihourf * 60 / dt) + 1
date_I = [date0 + timedelta(minutes=dt * x) for x in range(ndata)]

pyspedas.rbsp.hope(trange=[date1, date2])

cdf_path = (f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/moments/rel04/2015/"
            "rbspa_rel04_ect-hope-mom-l3_20150318_v7.1.0.cdf")
cdf_file = pycdf.CDF(cdf_path)
print(cdf_file)
