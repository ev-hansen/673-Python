# Description
2022 Summer GSFC work

This repo exists to host code that is being moved from IDL to Python. For the most part, I try to follow flake8 formatting standards (Warnings/Errors I ignore are W291, W504, E402, and E741).

[Flake8 etc setup for Visual Studio Code](https://medium.com/@jackklpan/auto-format-and-lint-by-black-isort-flake8-in-vs-visual-studio-code-a62a3f5d940e)
[Flake8 etc setup for vim](https://github.com/nvie/vim-flake8)

Cristian's IDL code (``calc_hope_pressure.pro``) uses formulas from [this paper](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019JA026695) to calculate pressure of ion species (H+, He+, O+).

My internship ended 2022/08/12 and code is not fully implemented and may not be maintained in the future. Feel free to make [issues](https://docs.github.com/en/issues/tracking-your-work-with-issues/about-issues) and [pull requests](https://docs.github.com/en/pull-requests), please keep in mind the [license](https://opensource.org/licenses/MIT). 


# Requirements
Python is required, a conda installation is reccomended. My current conda environment uses python 3.9.12.

# Setting up Conda and the .env file
1) Install conda (I suggest miniconda, but anaconda might be better depending on the system)
2) Create a new conda environment ([see this link](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)), or install everything into the ``base`` environment if you're only using it for this (not reccomended, but it's what I did). Regardless, install numpy in the conda environment you chose with ``conda install numpy``. 
3) Run ``pip install pyspedas python-dotenv tqdm cdflib`` in the conda environment. [The Conda blog](https://www.anaconda.com/blog/using-pip-in-a-conda-environment) suggests adding ``--upgrade-strategy only-if-needed`` to pip install commands. In all honesty I didn't because I forgot, but it's probably best practice.
4) After the packages install, create a file called ``.env`` in the same directory as the files from this repo.
5) In the ``.env`` file,  add a line that is your OS type in all caps (e.g. ``MAC``, ``WINDOWS``, or ``LINUX``), followed by ``_HOME_DIR = `` then the path to the folder of the files from this repo.
6) You can now run the Python code from this repo. Email ev.hansen@umbc.edu for assistance.

# Running code from the ``test code`` folder
1) Some of the test code requires spacepy's ``pycdf``. You need to follow instructions [from spacepy's website](https://spacepy.github.io/) to install on your machine, as there are some platform-specific instructions for the ``spacepy`` package. Also see [this page](https://spacepy.github.io/pycdf.html) when configuring ``pycdf``.
2) Copy the ``.env`` file from the project's parent directory, but replace the word ``HOME`` with ``TEST`` and make sure all the paths are for the ``test code`` folder instead of the parent.
3) You will need to include the ``CDF_LIB`` path in the ``.env`` file in the ``test code`` directory. To do this, add a line that is your OS type in all caps (e.g. ``MAC``, ``WINDOWS``, or ``LINUX``), followed by ``_CDF_LIB = `` then the path to your local CDF installation from step 1.
4) When using the test code, make sure to ``cd`` into that directory instead of running from the parent directory, as the parent directory has a different ``.env`` file
5) Some test code may have additional packages installed. I suggest checking the imports and installing missing packages via pip before running code in that directory.

# TODO (``calc_hope_pressure.py``)
- Correct fluxes (skip for now)
- ***Calculate Pressure*** (``calc_pressure`` in ``HopeCalculations``)
- Pitchangle avg flux
- ***Calculate mean energy***
- Smooth data
- Replace 0s w/ NANs
- ***Average Pressures***
- ***Create plot data***
- Change CDF dict data object to an Xarray for easier operations
- Restructure class to be compatable with other data operations, maybe multiple files?
- Rename class

# TODO (other)
- Look at other plasma parameters (other calculations, density, mean energy)
- Recreate loop (handler/wrapper)
- Recreate binning
- Plot data
- Recreate dial plot (may be possible w/ matplotlib)
- Plasma density plots?
- Solar wind + curvature of magnetic field (?)
