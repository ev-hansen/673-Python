from bs4 import BeautifulSoup
import requests

url = ('https://cdaweb.gsfc.nasa.gov/pub/data/rbsp/rbspa/l3/ect/hope/'
       'pitchangle/rel04/2013/?')
ext = 'cdf'
ymd_string = "20131030"


def listFD(url, ext=''):
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') if 
            node.get('href').endswith(ext)]


for file in listFD(url, ext):
    if (ymd_string in file):
        print(file)
