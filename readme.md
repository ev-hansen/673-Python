# Description
This repo exists to host code that is being moved from IDL to Python. For the most part, I try to follow flake8 formatting standards.


# Requirements
Python is required, I suggest using a Conda installation.

## Linux
- C and fortran compilers are needed (e.x. ``gcc`` and ``gfortran``)
- Conda packages needed are ``numpy``, ``scipy``, ``matplotlib``, and ``h5py``
- An ncurses library should be installed
- The CDF library can be obtained from https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_1/linux/ or https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf37_1/linux/. After installing, note the path of the ``lib/`` directory.
- See https://spacepy.github.io/install_linux.html for more installation guidance, as you may have to compile the CDF library yourself.

## Mac
- A C compiler is needed (run ``xcode-select --install``)
- conda packages needed are ``gfortran_osx-64``, ``numpy``, ``scipy``, ``matplotlib``, ``h5py``
- The CDF library can be obtained from https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_1/macosx/ or https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf37_1/macosx/. After installing, note the path of the ``lib/`` directory (you may have to view the contents of the app in your ``Applications`` folder).

## Windows
- Conda packages needed are ``m2w64-gcc-fortran``, ``libpython``, ``numpy``, ``scipy``, ``matplotlib``, ``h5py``
- The CDF library can be obtained from https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_1/windows/ or https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf37_1/windows/. After installing, note the path of the ``lib/`` directory.

# Setting up Conda and the Programming .env
1) In your conda environment, install the packages for your system listed above. After that, run ``pip install spacepy pyspedas python-dotenv tqdm``. The Conda site suggest adding ``--upgrade-strategy only-if-needed`` to pip install commands, I'm not sure how essential that is.
2) After that installs, create a file called ``.env`` in the same directory as the files from this repo.
3) In the ``.env``file, add a line that is your OS type in all caps (e.g. MAC, WINDOWS, or LINUX), followed by ``_CDF_LIB = `` then the path to the ``lib/`` directory you noted earlier.
4) In the ``.env`` file,  add a line that is your OS type in all caps (e.g. MAC, WINDOWS, or LINUX), followed by ``_HOME_DIR = `` then the path to the folder of the files from this repo.
5) You can now run the Python code from this repo. Email ev.hansen@umbc.edu for assistance.

# TODO (calc_hope_pressure.py)
- Correct fluxes
- Calculate Pressure
- Smooth data
- Replace 0s w/ NANs
- Average Pressures
- Create plot data