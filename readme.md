# Description
This repo exists to host code that is being moved from IDL to Python. For the most part, I try to follow flake8 formatting standards. My internship ends 2022/08/12 so code may not be maintained after then.


# Requirements
Python is required via a conda installation. My current conda environment uses python 3.9.12.

# Setting up Conda and the Programming .env
1) Install conda (I suggest miniconda, but anaconda might be better depending on the system)
2) Create a new conda environment, or install everything into the ``base`` environment if you're only using it for this. Regardless, install numpy in a conda environment with ``conda install numpy``. 
3) Run ``pip install pyspedas python-dotenv tqdm cdflib``. The Conda site suggest adding ``--upgrade-strategy only-if-needed`` to pip install commands, I'm not sure how essential that is.
4) After that installs, create a file called ``.env`` in the same directory as the files from this repo.
5) In the ``.env`` file,  add a line that is your OS type in all caps (e.g. MAC, WINDOWS, or LINUX), followed by ``_HOME_DIR = `` then the path to the folder of the files from this repo.
6) Some of the test code requires spacepy's pycdf. You need to follow instructions on https://spacepy.github.io/pycdf.html to install on your machine, there are some platform-dependent there.
7) You can now run the Python code from this repo. Email ev.hansen@umbc.edu for assistance.

# TODO (``calc_hope_pressure.py``)
- Correct fluxes (skip for now)
- *Calculate Pressure*
- Pitchangle avg flux
- *Calculate mean energy*
- Smooth data
- Replace 0s w/ NANs
- *Average Pressures*
- *Create plot data*
- Change CDF dict data object to an Xarray for easier
  operations
- Restructure class to be compatable with other 
  data operations, maybe multiple files?
- Rename class