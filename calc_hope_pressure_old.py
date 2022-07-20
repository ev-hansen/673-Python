#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Python version of Cristian's IDL code of the same name

This code is mostly ported code from the IDL calc_hope_pressure.pro code. 
There are some modifications, such as the use of an object (currently named
HopeCalculations) to store parameter data and perform calculations so
functionality can be imported and extended, as well as commandline argument
handling, so the file does not need ot be manually edited for parameters
when directly running the file. 
Use the --help flag for info on the parameters.

This code was written to work in a conda environment, it is suggested to set
one up. It is suggested to add “--upgrade-strategy only-if-needed” when using 
pip to install packages in your conda environment.

This code requires cdflib to be installed, see 
https://spacepy.github.io/pycdf.html#module-spacepy.pycdf. 
Evan suggests using 3.7.1, as 3.8.1 seems to have issues importing in python.
Make a '.env' file in the same directory as this file and specify the CDF_LIB
location.  See above link for what to set CDF_LIB to in your .env file. 

"""

__authors__ = ["Evan Hansen", "Cristian Ferradas"]
# only including Cristian's email, as Evan is an intern.
__contact__ = "cristian.ferradasalva@nasa.gov" 

# add your name here if you modify the code for a new version, 
# perhaps w/ Cristian's approval first
__credits__ = [["Evan Hansen", "Python code"], 
               ["Cristian Ferradas", "IDL code"]] 

# perhaps change these if/when code is changed / is production-ready
__date__ = "2022/06/08 - 2022/08/12"
__status__ = "Development"
__version__ = "0.0.1"


############################################################
# TODO: - Correct fluxes
# TODO: - Calculate Pressure
# TODO: - Smooth data
# TODO: - Replace 0s w/ NANs
# TODO: - Average Pressures
# TODO: - Create plot data
############################################################


# Standard Imports #

import datetime as dt
import argparse
import os
import platform
from typing import List, Dict, Tuple, Any


# Third-Party Inports #

from dotenv import load_dotenv
load_dotenv()

# CDF Operations
from spacepy import pycdf
from pyspedas.rbsp import hope

# Calculations and array/list ops
import numpy

from tqdm import tqdm


# Find OS-specific CDF_LIB path in .env file. Make sure you
# set whichever is patform-relevant to you
if (platform.system() == "Windows"):
    HOME_DIR = os.environ["WIN_HOME_DIR"]
    CDF_LIB = os.environ["WIN_CDF_LIB"]
    os.environ["CDF_LIB"] = CDF_LIB
elif (platform.system() == "Linux"):
    HOME_DIR = os.environ["LINUX_HOME_DIR"]
    CDF_LIB = os.environ["LINUX_CDF_LIB"]
    os.environ["CDF_LIB"] = CDF_LIB
elif (platform.system() == "Darwin"):
    HOME_DIR = os.environ["MAC_HOME_DIR"]
    CDF_LIB = os.environ["MAC_CDF_LIB"]
    os.environ["CDF_LIB"] = CDF_LIB


# Main File #


class HopeCalculations:
    """
    Class to process HOPE data and calculations
    """
    def __init__(self, given_time_s_str: str, given_time_e_str: str, 
                 given_probe: str, given_level: str, given_factor, 
                 given_low_energy: List[float], given_up_energy: List[float], 
                 given_pot_corr: int, given_relat: int, given_swindow: int):
        """
        Constructor for the HopeCalculations class. 

            Parameters: 
                given_time_s_str (str): start time, format yyyy-mm-dd/hh:mm:ss
                given_time_e_str (str): end time, format yyyy-mm-dd/hh:mm:ss
                given_probe (char or str): probe 'a' or 'b'
                given_level (str): level of data (l1, l2, l3, or l4)
                given_factor: 
                given_low_energy (float array): [i+, e-] in eV 
                given_up_energy (float array): [i+, e-] in eV
                given_pot_corr (int): Use sc potential corrected fluxes? 
                                                            0:No, 1:Yes
                given_relat (int): Use relativistic energy equation? 
                                                            0:No, 1:Yes 
                given_swindow (int): Smoothing window in sec; 
                                                swindow=0: No smoothing
        """
        self.__time_s_str = given_time_s_str
        self.__time_e_str = given_time_e_str
        self.__time_range = [given_time_s_str, given_time_e_str]
        self.__probe = given_probe
        self.__datatype = "pitchangle"
        self.__instrument = "ect" 
        self.__component = "hope"
        self.__level = given_level
        self.__relat = given_relat
        datetimes = self.init_datetime_objects(
            self.__time_s_str,
            self.__time_e_str)
        self.__time_s_datetime = datetimes[0]
        self.__time_e_datetime = datetimes[1]
        self.__day_diff = self.calc_day_diff(
            self.__time_s_datetime,
            self.__time_e_datetime)
        self.__factor = given_factor
        self.__low_energy = given_low_energy
        self.__up_energy = given_up_energy
        self.__pot_corr = given_pot_corr
        self.__swindow = given_swindow

    def init_datetime_objects(self, given_time_s_str: str, 
                              given_time_e_str: str) -> Tuple[dt.datetime]:
        """
        Return datetime objects for the start and end times.

            Parameters:
                given_time_s_str (str): start time, format yyyy-mm-dd/hh:mm:ss
                given_time_e_str (str): end time, format yyyy-mm-dd/hh:mm:ss

            Returns:
                time_s_datetime (datetime object): start time
                time_e_datetime (datetime object): end time
        """
        time_s_datetime = dt.datetime.strptime(given_time_s_str, 
                                               '%Y-%m-%d/%H:%M:%S')
        time_e_datetime = dt.datetime.strptime(given_time_e_str, 
                                               '%Y-%m-%d/%H:%M:%S')
        return time_s_datetime, time_e_datetime

    def calc_day_diff(self, given_datetime_s: dt.datetime, 
                      given_datetime_e: dt.datetime) -> int:
        """
        Find the difference of days between start and end time.

            Parameters:
                given_datetime_s (datetime object): start time
                given_datetime_e (datetime object): end time

            Returns:
                day_diff (int): difference of days
        """
        given_datetime_delta = given_datetime_e - given_datetime_s
        day_diff = given_datetime_delta.days
        return day_diff

    def datetime_range_list(self) -> list[dt.datetime]:
        """
        Generate a list of datetime objects between the start and end time.

            Returns:
                date_list (list): list of datetime objects
        """
        date_list = [self.__time_s_datetime +
                     dt.timedelta(days=i) for i in range(self.__day_diff + 1)]
        return date_list

    def get_cdf_objs(self, 
                     given_datetime_list: 
                         list[dt.datetime]) -> List[pycdf.CDF]:
        """
        Gets a lsit of CDF objects for the given time range.
            Parameters:
                given_datetime_list (list): List of datetime objects

            Returns:
                cdf_obj_list c(list): cdf file objects in a list
        """
        cdf_obj_list = []
        for date_obj_i in given_datetime_list:
            now = dt.datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print(f"{current_time} getting file for {date_obj_i}")
            cdf_obj_i = self.get_cdf(date_obj_i)
            cdf_obj_list.append(cdf_obj_i)
        return cdf_obj_list

    def get_cdf(self, given_datetime_obj: dt.datetime) -> pycdf.CDF:
        """
        Gets cdf file object, either from disk or the internet.

            Parameters:
                given_datetime_obj (datetime object): datetime object

            Returns:
                cdf_obj (cdf object): cdf file object
        """
        probe = self.__probe
        datatype = self.__datatype

        if (datatype == "moments"):
            datatype_abbr = "mom"
        elif (datatype == "pitchangle"):
            datatype_abbr = "pa"

        rel = "rel04"
        instrument = self.__instrument
        component = self.__component
        level = self.__level

        yr_str = given_datetime_obj.strftime("%Y")  # year
        m_str = given_datetime_obj.strftime("%m")   # month
        d_str = given_datetime_obj.strftime("%d")   # day

        # given_datetime_obj.strftime("%Y%m%d") probably could have been done 
        # but Evan wanted to mirror some formatting from the IDL version
        ymd_str = yr_str + m_str + d_str

        # example: rbspa_rel04_ect-hope-mom-l3_20140317_v7.1.0
        nm_hd = (f"rbsp{probe}_{rel}_{instrument}-{component}-"
                 f"{datatype_abbr}-{level}_")
        fln_tmp = f"{nm_hd}{ymd_str}_"

        cdf_file_dir = (f"{HOME_DIR}/rbsp_data/rbsp{probe}/{level}/"
                        f"{instrument}/{component}/pitchangle/"
                        f"{rel}/{yr_str}")

        # Check of cdf file exists
        error = 0

        try:
            # search for file directory
            file_list = os.listdir(cdf_file_dir)
            try:
                # Search for file in dir
                str_match = list(filter(lambda x: fln_tmp in x, file_list))
                cdf_file_path = f"{cdf_file_dir}/{str_match[0]}"
            except IndexError:
                print(f"Could not find {fln_tmp} in directory {cdf_file_dir}")
                error += 1
        except FileNotFoundError:
            print(f"Could not find directory {cdf_file_dir}")
            error += 2

        if (error > 0):
            print(f"Attempting to download {fln_tmp}")
            hope(trange=[f"{yr_str}-{m_str}-{d_str}", 
                         f"{yr_str}-{m_str}-{d_str}"], probe=probe,
                 datatype=datatype, downloadonly=True)  
            file_list = os.listdir(cdf_file_dir)
            str_match = list(filter(lambda x: fln_tmp in x, file_list))
            cdf_file_path = f"{cdf_file_dir}/{str_match[0]}"

        cdfObject = pycdf.CDF(cdf_file_path)
        return cdfObject

    def read_cdf_data(self, given_cdf_object: pycdf.CDF) -> Dict[str, Any]:
        """Gather data from the given CDF object

            Args:
                given_cdf_object (CDF): A pitchangle CDF object

            Returns:
                cdf_data (dict): A dict w/ CDF data
                    ion_data_array (numpy.arry): An array w/ epoch, ion data, 
                        mode, fpdu, fhedu, and fodu values from the given CDF
                        file
                    ele_data_array (numpy.array): An array w/ epoch, 
                        electron data, mode, and fedu values from the given CDF 
                        file.

        """

        start_time = self.__time_s_datetime
        end_time = self.__time_e_datetime

        # Datetime objs represening timestamps
        t_ion = self.reform(given_cdf_object["Epoch_Ion"])  # Ion
        t_ele = self.reform(given_cdf_object["Epoch_Ele"])  # Electron

        # 
        e_data_ion = given_cdf_object["HOPE_ENERGY_Ion"] 
        e_data_ele = given_cdf_object["HOPE_ENERGY_Ele"]

        # Mode of data collection
        # 0 is Apogee Mode, 1 is Perigee Mode, 2 is Burst Mode
        # (0) Apogee mode: normal operation.
        # (1) Perigee mode: minimum energy channel
        # (2) burst mode: subset of energies sampled at rapid cadence.
        # burst mode is only used for electron operations.
        mode_ion = self.reform(given_cdf_object["Mode_Ion"])
        mode_ele = self.reform(given_cdf_object["Mode_Ele"])

        # Data of the "species," function of pitchangle"
        fpdu_data = given_cdf_object["FPDU"]  # proton flux
        fhedu_data = given_cdf_object["FHEDU"]  # hydrogen flux
        fodu_data = given_cdf_object["FODU"]  # oxygen flux
        fedu_data = given_cdf_object["FEDU"]  # electron fluxe

        pitchangle_data = numpy.array(given_cdf_object["PITCH_ANGLE"])
        # given_cdf_object.close()

        fill_value = 0
        apogee_mode = 0

        # read ion related data, replace fill values in species data w/ NaN
        # data from Parigee and Burst are also replaced w/ NaNs
        # epoch type: dt.datetime
        # data type: numpy.array[numpy.array[float]]
        # mode type: numpy.array[int]
        # fpdu, fhedu, fodu type: numpy.array[numpy.array[float]]
        ion_data_array = numpy.array([
            {"epoch": t_ion[i],  # type: dt.datetime
             "data": numpy.array(list(e_data_ion[i])),
             "mode": numpy.array(mode_ion[i]), 
             "fpdu": numpy.array([[
                 k if ((k >= fill_value) and (mode_ion[i]) == apogee_mode) 
                 else float("NaN")
                 for k in fpdu_data[i][j]] for j in range(len(
                     fpdu_data[i]))]),
             "fhedu": numpy.array([[
                 k if ((k >= fill_value) and (mode_ion[i]) == apogee_mode) 
                 else float("NaN")
                 for k in fhedu_data[i][j]] for j in range(len(
                     fhedu_data[i]))]),
             "fodu": numpy.array([[
                 k if ((k >= fill_value) and (mode_ion[i]) == apogee_mode) 
                 else float("NaN")
                 for k in fodu_data[i][j]] for j in range(len(
                     fodu_data[i]))])
             } for i in tqdm(range(len(t_ion)), desc="ion") 
            if (start_time < t_ion[i] < end_time)])

        ele_data_array = numpy.array([
            {"epoch": t_ele[i],
             "data": numpy.array(list(e_data_ele[i])),
             "mode": numpy.array(mode_ele[i]),
             "fedu": numpy.array([[
                 k if ((k >= fill_value) and (mode_ele[i]) == apogee_mode) 
                 else float("NaN")
                 for k in fedu_data[i][j]] for j in range(len(
                     fedu_data[i]))])
             } for i in tqdm(range(len(t_ele)), desc="ele")
            if (start_time < t_ele[i] <= end_time)])     

        cdf_data = {"ion_data_dict": ion_data_array, 
                    "ele_data_dict": ele_data_array,
                    "pitchangle_data": pitchangle_data}

        return cdf_data

    def pressure_calculations(self, 
                              cdf_data: Dict[dt.datetime, 
                                             Dict[str, Any]]):

        print()

    def transpose(self, given_list: List[Any]) -> numpy.array:
        """Transposes a list (swaps rows/columns)

            Parameters:
                given_list (list): A 2D list

            Returns:
                transposed_list (arr): given_list transposed as an arr
        """
        arr = numpy.array(given_list)
        transposed = arr.T
        return transposed

    def reform(self, given_list: List[Any]) -> numpy.array:
        """Acts similar to how the IDL reform function works in the original code

            Parameters:
                given_list (list): A list, potentially w/ a dimension of size 1

            Returns:
                reformed (arr): Numpy arr of list w/o dimension of size 1
        """
        arr = numpy.array(given_list)
        reformed = arr.squeeze()
        return reformed

    def wrapper(self):
        datetime_list = self.datetime_range_list()
        cdf_objs = self.get_cdf_objs(datetime_list)
        cdf_objs_data = numpy.array([])

        for i in tqdm(range(len(cdf_objs)), desc="getting data"):
            cdf_obj = cdf_objs[i]
            now = dt.datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print(f"{current_time} gathering information from " 
                  f"{datetime_list[i]}")
            this_cdf_data = self.read_cdf_data(cdf_obj)
            cdf_obj.close()
            cdf_objs_data = numpy.append(cdf_objs_data, this_cdf_data)

        self.pressure_calculations(cdf_objs_data)

    # Getters, mostly for testing purposes
    @property
    def time_s_str(self) -> str:
        return self.__time_s_str  

    @property
    def time_e_str(self) -> str:
        return self.__time_e_str 

    @property
    def time_range(self) -> list[str]:
        return self.__time_range

    @property
    def probe(self) -> str:
        return self.__probe

    @property
    def instrument(self) -> str:
        return self.__instrument

    @property
    def component(self) -> str:
        return self.__component

    @property
    def relat(self) -> int:
        return self.__relat

    @property
    def time_s_datetime(self) -> dt.datetime:
        return self.__time_s_datetime  

    @property
    def time_e_datetime(self) -> dt.datetime:
        return self.__time_e_datetime  

    @property
    def day_diff(self) -> int:
        return self.__day_diff   

    @property
    def factor(self):
        return self.__factor   

    @property
    def low_energy(self) -> list[float]:
        return self.__low_energy  

    @property
    def up_energy(self) -> list[float]:
        return self.__up_energy  

    @property
    def pot_corr(self) -> int:
        return self.__pot_corr

    @property
    def swindow(self) -> int:
        return self.__swindow


###############################################################################


def init_parser() -> argparse.ArgumentParser:
    """
    Initializes console command parser

        Returns: parser: parser object
    """
    # Initialize parser
    parser = argparse.ArgumentParser(
        description='Calculate the pressure from a ECT-HOPE cdf')
    parser.add_argument('-s', '--start_time', 
                        help='Start of the time range, YYYY-MM-DD/HH:MM:SS.',
                        required=False, type=str, 
                        default='2013-03-17/00:00:00')
    parser.add_argument('-e', '--end_time',
                        help='End of the time range, YYYY-MM-DD/HH:MM:SS.',
                        required=False, type=str, 
                        default='2013-03-20/00:00:00')
    parser.add_argument('-p', '--probe',
                        help='RBSP probe: a or b.', 
                        required=False, type=str, default='a')
    parser.add_argument('-l', '--level',
                        help='Data level (l1 - l4)', 
                        required=False, type=str, default='l3')
    parser.add_argument('-r', '--relat', 
                        help='Use relativistic energy equation? 0:No, 1:Yes.', 
                        required=False, type=int, default=1)
    parser.add_argument('-f', '--factor', 
                        help='Factor to be multiplied to the particle fluxes.', 
                        required=False, type=float, default=1.0)
    parser.add_argument('-l0', '--low_energy_0', 
                        help='Lower energy limit i+ for the particle fluxes.', 
                        required=False, type=float, default=1e2) 
    parser.add_argument('-l1', '--low_energy_1', 
                        help='Lower energy limit e+ for the particle fluxes.', 
                        required=False, type=float, default=1e2) 
    parser.add_argument('-u0', '--up_energy_0', 
                        help='Upper energy limit i+ for the particle fluxes.', 
                        required=False, type=float, default=5.5e4)
    parser.add_argument('-u1', '--up_energy_1', 
                        help='Upper energy limit e+ for the particle fluxes.', 
                        required=False, type=float, default=5.5e4)
    parser.add_argument('-c', '--pot_corr', 
                        help='Use sc potential corrected fluxes? 0:No, 1:Yes.', 
                        required=False, type=int, default=1)
    parser.add_argument('-w', '--swindow', 
                        help='Smoothing window in sec; 0 is No smoothing', 
                        required=False, type=int, default=60)

    return parser


def main():
    """
    Main function to run the program via command line. 
    Functionality included to parse and validate command line arguments so 
    modifications to the program are not required.
    """    

    parser = init_parser()

    args = parser.parse_args()
    args_dict = vars(args)

    time_s_str = args_dict['start_time']
    time_e_str = args_dict['end_time']
    level = args_dict['level']
    probe = args_dict['probe']
    factor = args_dict['factor']
    low_energy = [args_dict['low_energy_0'], args_dict['low_energy_1']]
    up_energy = [args_dict['up_energy_0'], args_dict['up_energy_1']]
    pot_corr = args_dict['pot_corr']
    relat = args_dict['relat']
    swindow = args_dict['swindow']

    main_calculations = HopeCalculations(time_s_str, time_e_str, probe, level, 
                                         relat, factor, low_energy, up_energy, 
                                         pot_corr, swindow)

    main_calculations.wrapper()


# Run main() if called from command line
if (__name__ == "__main__"):
    main()
