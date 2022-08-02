import os
from dotenv import load_dotenv
load_dotenv()

print(os.environ["CDF_LIB"])
print(os.environ["TEST_DIR"])
print(os.listdir(os.environ["TEST_DIR"]))
