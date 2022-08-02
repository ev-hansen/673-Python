import os

# dotenv and environment variable stuff
from dotenv import load_dotenv
load_dotenv()
# set up .env file to specify the home directory
TEST_DIR = os.environ["TEST_DIR"]

from spacepy import pycdf

import pprint
pp = pprint.PrettyPrinter(indent=4)

with pycdf.CDF(f"{TEST_DIR}/rbsp_data/rbspa/l3/ect/hope/moments/rel04/2015/"
               "rbspa_rel04_ect-hope-mom-l3_20151126_v7.1.0.cdf") as cdf_file:
    pp.pprint(cdf_file)
