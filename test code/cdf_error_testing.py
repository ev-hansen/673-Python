from dotenv import load_dotenv
load_dotenv()

from spacepy import pycdf


try:
    cdf = pycdf.CDF("doesnt_exist")
except pycdf.CDFError as e:
    print(e)

print("nice")
    