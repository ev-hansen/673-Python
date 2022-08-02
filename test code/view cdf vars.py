import os
import platform
from spacepy import pycdf


# dotenv and environment variable stuff
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


TEST_CDF = (f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/moments/rel04/2015/"
            "rbspa_rel04_ect-hope-mom-l3_20151126_v7.1.0.cdf")

cdf_file = pycdf.CDF(TEST_CDF)
t_ion = cdf_file["Epoch_Ion"]
t_ele = cdf_file["Epoch_Ele"]

print(f"\nepoch_ion: {t_ion} len: {len(t_ion)}\n")
print(f"\nepoch_ele: {t_ele} len: {len(t_ele)}\n")

item = str(t_ion[0])
print(item.split())

cdf_file.close()
