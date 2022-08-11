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
__license__ = "MIT"


###############################################################################
# TODO: - Correct fluxes (skip for now)
# TODO: - *Calculate Pressure*
# TODO: - Pitchangle avg flux
# TODO: - *Calculate mean energy*
# TODO: - Smooth data
# TODO: - Replace 0s w/ NANs
# TODO: - *Average Pressures*
# TODO: - *Create plot data*
# TODO: - Change CDF dict data object to an Xarray for easier
#         operations
# TODO: - Restructure class to be compatable with other 
#         data operations, maybe multiple files?
# TODO: - Rename class
###############################################################################


# Standard Imports #

import datetime as dt
import argparse
import os
import platform
from typing import List, Dict, Tuple, Any
import contextlib
import math


# Third-Party Inports #

from dotenv import load_dotenv
load_dotenv()

# CDF Operations
import cdflib
import xarray as xr
from xarray.core.dataset import Dataset as XarrDataset
from pyspedas.rbsp import hope, efw

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy

from tqdm import tqdm

# matplotlib settings, feel free to change
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams.update({'font.size': 10})
mpl.rcParams["figure.figsize"] = [14.1, 9.7]
mpl.rcParams["figure.dpi"] = 226

# Find OS-specific paths from .env file. Make sure you
# set whichever is patform-relevant to you
if (platform.system() == "Windows"):
    HOME_DIR = os.environ["WIN_HOME_DIR"]
elif (platform.system() == "Linux"):
    HOME_DIR = os.environ["LINUX_HOME_DIR"]
elif (platform.system() == "Darwin"):
    HOME_DIR = os.environ["MAC_HOME_DIR"]


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
                         list[dt.datetime]) -> List[XarrDataset]:
        """
        Gets a lsit of CDF objects for the given time range.
            Parameters:
                given_datetime_list (list): List of datetime objects

            Returns:
                cdf_obj_list c(list): cdf file objects in a list
        """
        hope_cdf_obj_list = []
        efw_cdf_obj_list = []
        for date_obj_i in tqdm(given_datetime_list, desc="Getting files"):
            now = dt.datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print(f"{current_time} getting files for {date_obj_i}")
            hope_cdf_obj_i, efw_cdf_obj_i = self.get_cdfs(date_obj_i)
            hope_cdf_obj_list.append(hope_cdf_obj_i)
            efw_cdf_obj_list.append(efw_cdf_obj_i)
        return hope_cdf_obj_list, efw_cdf_obj_list

    def get_cdfs(self, given_datetime_obj: dt.datetime):
        hope_cdf = self.get_hope_cdf(given_datetime_obj)
        efw_cdf = self.get_efw_cdf(given_datetime_obj)
        return hope_cdf, efw_cdf

    def get_hope_cdf(self, given_datetime_obj: dt.datetime) -> XarrDataset:
        """
        Gets hope cdf file object, either from disk or the internet.

            Parameters:
                given_datetime_obj (datetime object): datetime object

            Returns:
                cdf_obj (cdf object): cdf file object
        """

        probe = self.__probe
        datatype = "pitchangle"
        datatype_abbr = "pa"

        rel = "rel04"
        instrument = "ect"
        component = "hope"
        level = self.__level

        yr_str = given_datetime_obj.strftime("%Y")  # year
        m_str = given_datetime_obj.strftime("%m")   # month
        d_str = given_datetime_obj.strftime("%d")   # day

        # given_datetime_obj.strftime("%Y%m%d") probably could have been done 
        # but Evan wanted to mirror some formatting from the IDL version
        ymd_str = yr_str + m_str + d_str

        # example: rbspa_rel04_ect-hope-mom-l3_20140317_v7.1.0
        hope_nm_hd = (f"rbsp{probe}_{rel}_{instrument}-{component}-"
                      f"{datatype_abbr}-{level}_")
        hope_fln_tmp = f"{hope_nm_hd}{ymd_str}_"

        hope_cdf_file_dir = (f"{HOME_DIR}/rbsp_data/rbsp{probe}/{level}/"
                             f"{instrument}/{component}/pitchangle/"
                             f"{rel}/{yr_str}")

        # Check if hope cdf file exists
        hope_error = 0
        try:
            # search for file directory
            file_list = os.listdir(hope_cdf_file_dir)
            try:
                # Search for file in dir
                str_match = list(filter(lambda x: hope_fln_tmp in x, 
                                        file_list))
                cdf_file_path = f"{hope_cdf_file_dir}/{str_match[0]}"
            except IndexError:
                print(f"Could not find {hope_fln_tmp} in directory "
                      f"{hope_cdf_file_dir}")
                hope_error += 1
        except FileNotFoundError:
            print(f"Could not find directory {hope_cdf_file_dir}")
            hope_error += 2

        if (hope_error > 0):
            print(f"Attempting to download {hope_fln_tmp}")
            hope(trange=[f"{yr_str}-{m_str}-{d_str}", 
                         f"{yr_str}-{m_str}-{d_str}"], probe=probe,
                 datatype=datatype, downloadonly=True)  
            file_list = os.listdir(hope_cdf_file_dir)
            str_match = list(filter(lambda x: hope_fln_tmp in x, file_list))
            cdf_file_path = f"{hope_cdf_file_dir}/{str_match[0]}"

        dataset = cdflib.cdf_to_xarray(cdf_file_path, to_unixtime=True, 
                                       fillval_to_nan=True)
        return dataset

    def get_efw_cdf(self, given_datetime_obj: dt.datetime):
        """
        Gets hope cdf file object, either from disk or the internet.

            Parameters:
                given_datetime_obj (datetime object): datetime object

            Returns:
                cdf_obj (cdf object): cdf file object
        """

        probe = self.__probe

        rel = "rel04"
        instrument = "efw"
        level = self.__level

        yr_str = given_datetime_obj.strftime("%Y")  # year
        m_str = given_datetime_obj.strftime("%m")   # month
        d_str = given_datetime_obj.strftime("%d")   # day

        # given_datetime_obj.strftime("%Y%m%d") probably could have been done 
        # but Evan wanted to mirror some formatting from the IDL version
        ymd_str = yr_str + m_str + d_str

        # example: rbspa_rel04_ect-hope-mom-l3_20140317_v7.1.0
        efw_fln_tmp = f"rbsp{probe}_{instrument}-{level}_{ymd_str}_"

        efw_cdf_file_dir = (f"{HOME_DIR}/rbsp_data/rbsp{probe}/{level}/"
                            f"{instrument}/{yr_str}/")

        # Check if hope cdf file exists
        efw_error = 0
        try:
            # search for file directory
            file_list = os.listdir(efw_cdf_file_dir)
            try:
                # Search for file in dir
                str_match = list(filter(lambda x: efw_fln_tmp in x, 
                                        file_list))
                efw_cdf_file_path = f"{efw_cdf_file_dir}/{str_match[0]}"
            except IndexError:
                print(f"Could not find {efw_fln_tmp} in directory "
                      f"{efw_cdf_file_dir}")
                efw_error += 1
        except FileNotFoundError:
            print(f"Could not find directory {efw_cdf_file_dir}")
            efw_error += 2

        if (efw_error > 0):
            print(f"Attempting to download {efw_fln_tmp}")
            efw(trange=[f"{yr_str}-{m_str}-{d_str}", 
                        f"{yr_str}-{m_str}-{d_str}"], probe=probe,
                downloadonly=True)  
            file_list = os.listdir(efw_cdf_file_dir)
            str_match = list(filter(lambda x: efw_fln_tmp in x, file_list))
            efw_cdf_file_path = f"{efw_cdf_file_dir}/{str_match[0]}"

        dataset = cdflib.CDF(efw_cdf_file_path)
        return dataset

    def read_cdf_data(self, 
                      given_datasets: List[cdflib.CDF]) -> Dict[str, Any]:
        """Gather data from the given CDF object

            Args:
                given_cdf_object (CDF): A pitchangle CDF object

            Returns:
                cdf_data (dict): A dict w/ CDF data
                    ion_data_array (List[Dict]): A list of Dicts w/ epoch, ion 
                        data, mode, fpdu, fhedu, and fodu values from the given
                        CDF file
                    ele_data_array (List[Dict]): A list of Dicts w/ epoch, 
                        electron data, mode, and fedu values from the given CDF 
                        file.

        """

        start_time = self.__time_s_datetime
        end_time = self.__time_e_datetime

        factor = self.__factor
        low_energy = self.__low_energy
        up_energy = self.__up_energy

        pi = 3.141592653589793

        # Datetime objs represening timestamps
        # 

        # Mode of data collection
        # 0 is Apogee Mode, 1 is Perigee Mode, 2 is Burst Mode
        # (0) Apogee mode: normal operation.
        # (1) Perigee mode: minimum energy channel
        # (2) burst mode: subset of energies sampled at rapid cadence.
        # burst mode is only used for electron operations.

        apogee_mode = 0

        # read ion related data, replace fill values in species data w/ NaN
        # data from Parigee and Burst are also replaced w/ NaNs
        # epoch type: dt.datetime
        # data type: XarrDataset[XarrDataset[float]]
        # mode type: XarrDataset[int]
        # fpdu, fhedu, fodu type: XarrDataset[XarrDataset[float]]

        ion_data_list = [
            {'epoch': 
                dt.datetime.fromtimestamp(float(
                    given_datasets[i]['Epoch_Ion'].squeeze()[j].values)),
                'energy': 
                    given_datasets[i]['HOPE_ENERGY_Ion']
                    .where(low_energy[0] <= 
                           given_datasets[i]['HOPE_ENERGY_Ion']
                           .HOPE_ENERGY_Ion_dim).where(
                               given_datasets[i]['HOPE_ENERGY_Ion']
                               .HOPE_ENERGY_Ion_dim <= up_energy[0])[j],
                'fpdu': 
                    given_datasets[i]['FPDU'][j],
                'daty_avg_int_H1': 
                    given_datasets[i]['FPDU'][j].transpose() * factor,
                'fhedu': 
                    given_datasets[i]['FHEDU'][j],
                'daty_avg_int_He1': 
                    given_datasets[i]['FHEDU'][j].transpose() * factor,
                'fodu': 
                    given_datasets[i]['FODU'][j],
                'daty_avg_int_O1': 
                    given_datasets[i]['FODU'][j].transpose() * factor
             } for i in range(len(given_datasets))
            for j in tqdm(range(given_datasets[i]['Epoch_Ion'].size), 
                          desc=f"Ion data #{i+1}") 
            if ((start_time <= dt.datetime.fromtimestamp(float(
                given_datasets[i]['Epoch_Ion'].values[j])) <= end_time) and 
                (given_datasets[i]['Mode_Ion'][j].values == apogee_mode))]

        ele_data_list = [
            {'epoch':
                dt.datetime.fromtimestamp(float(
                    given_datasets[i]['Epoch_Ele'].squeeze()[j].values)),
             'energy': 
                 given_datasets[i]['HOPE_ENERGY_Ele']
                 .where(low_energy[1] <= 
                        given_datasets[i]['HOPE_ENERGY_Ele']
                        .HOPE_ENERGY_Ele_dim).where(
                            given_datasets[i]['HOPE_ENERGY_Ele']
                            .HOPE_ENERGY_Ele_dim <= up_energy[0])[j],
             'fedu': 
                 given_datasets[i]['FEDU'][j]
             } for i in range(len(given_datasets))
            for j in tqdm(range(given_datasets[i]['Epoch_Ele'].size),
                          desc=f"Ele data #{i+1}") 
            if ((start_time <= dt.datetime.fromtimestamp(float(
                given_datasets[i]['Epoch_Ele'].values[j])) <= end_time) and
                (given_datasets[i]['Mode_Ele'][j].values == apogee_mode))]

        cdf_data = {"ion_data_list": 
                    ion_data_list, 
                    "ele_data_list": 
                        ele_data_list,
                    "pitchangle_data": 
                        given_datasets[0]['PITCH_ANGLE'] * (pi / 180)}

        return cdf_data

    def correct_fluxes(self, given_cdf_data: List[Dict[str, Any]]):
        print()

    def calc_pressure(self, cdf_data: List[Dict[str, 
                                                Any]]) -> List[Dict[str, Any]]:
        """Perform the calculations, mostly cgs but also ergs

        Args:
            cdf_data (List[Dict[str, Any]]): List storing relevant CDF file 
                                            data, possibly spanning multiple 
                                            days
        """
        ion_data = cdf_data['ion_data_list']
        ele_data = cdf_data['ele_data_list']
        c_ion_time = len(ion_data)
        c_ele_time = len(ele_data)
        pa_arr = cdf_data['pitchangle_data']
        relat = self.__relat

        c = 2.9979e10  # Speed of light, cm/s
        pi = 3.141592653589793

        # Mass, grams
        mass_ele = 9.11e-28
        mass_pro = 1.67e-24
        mass_hel = 4 * mass_pro
        mass_oxy = 16 * mass_pro

        # Energy
        ev_to_erg = 1.60219e-12  # erg/eV

        en_ion_erg = [ion_data[i]['energy'] * ev_to_erg 
                      for i in tqdm(range(len(cdf_data['ion_data_list'])),
                                    desc='ion ev to erg')]
        en_ele_erg = [ele_data[i]['energy'] * ev_to_erg 
                      for i in tqdm(range(len(cdf_data['ele_data_list'])),
                                    desc='ele ev to erg')]

        en_ion_kev = [[ion_data[i]['energy'][j] / 1e3 
                       for j in range(len(ion_data[i]['energy']))]
                      for i in tqdm(range(c_ion_time), desc='ion to keV')]
        en_ele_kev = [[ele_data[i]['energy'][j] / 1e3
                       for j in range(len(ele_data[i]['energy']))]
                      for i in tqdm(range(c_ele_time), desc='ele to keV')]
        # need to speed up this section but might have to work on that list 
        # comprehension
        # 8:29
        # why is this faster than the list compr. below? Even if by a few secs.
        del_ion_en = [[float('NaN')] * 72] * c_ion_time

        for it in tqdm(range(c_ion_time), desc='del_ion'):
            for ie in range(70):
                del_ion_en[it][ie + 1] = [(en_ion_kev[it][ie + 1] - 
                                           en_ion_kev[it][ie - 1]) / 2]
            del_ion_en[it][0] = en_ion_kev[it][1] - en_ion_kev[it][0]
            del_ion_en[it][71] = en_ion_kev[it][71] - en_ion_kev[it][70]

        # 8:45
        """del_ion_en = [(en_ion_kev[it][1] - en_ion_kev[it][0] if ie == 0 else 
                      en_ion_kev[it][71] - en_ion_kev[it][70] if ie == 71 else 
                      en_ion_kev[it][ie] - en_ion_kev[it][ie - 1] / 2)
                      for it in tqdm(range(len(ion_data)), desc='del_ion') 
                      for ie in range(72)]"""

        del_ele_en = [[float('NaN')] * 72] * c_ele_time

        for it in tqdm(range(c_ele_time), desc='del_ele'):
            for ie in range(70):
                del_ele_en[it][ie + 1] = (en_ele_kev[it][ie + 1] - 
                                          en_ele_kev[it][ie - 1]) / 2
            del_ele_en[it][0] = en_ele_kev[it][1] - en_ele_kev[it][0]
            del_ele_en[it][71] = en_ele_kev[it][71] - en_ele_kev[it][70]

        gam_pro = [1 + en_ion_erg[i] / (mass_pro * (c * c))
                   for i in tqdm(range(len(en_ion_erg)), desc='gam_pro')]
        gam_hel = [1 + en_ion_erg[i] / (mass_hel * (c * c))
                   for i in tqdm(range(len(en_ion_erg)), desc='gam_hel')]
        gam_oxy = [1 + en_ion_erg[i] / (mass_oxy * (c * c))
                   for i in tqdm(range(len(en_ion_erg)), desc='gam_oxy')]
        gam_ele = [1 + en_ele_erg[i] / (mass_ele * (c * c))
                   for i in tqdm(range(len(en_ele_erg)), desc='gam_ele')]

        npa = len(pa_arr)

        del_pa = [pa_arr[1] - pa_arr[0] if ii == 0 else 
                  pa_arr[npa - 1] - pa_arr[npa - 2] if ii == npa - 1 else
                  (pa_arr[ii + 1] - pa_arr[ii - 1]) / 2.0 
                  for ii in tqdm(range(0, npa), desc='del_pa')]

        h_perp = [
            [float('NaN')] * len(ion_data[0]['daty_avg_int_H1'])] * c_ion_time
        h_para = [
            [float('NaN')] * len(ion_data[0]['daty_avg_int_H1'])] * c_ion_time
        he_perp = [
            [float('NaN')] * len(ion_data[0]['daty_avg_int_He1'])] * c_ion_time
        he_para = [
            [float('NaN')] * len(ion_data[0]['daty_avg_int_He1'])] * c_ion_time
        o_perp = [
            [float('NaN')] * len(ion_data[0]['daty_avg_int_O1'])] * c_ion_time
        o_para = [
            [float('NaN')] * len(ion_data[0]['daty_avg_int_O1'])] * c_ion_time
        e_perp = [
            [float('NaN')] * len(ele_data[0]['fedu'])] * c_ele_time
        e_para = [
            [float('NaN')] * len(ele_data[0]['fedu'])] * c_ele_time

        p_perp_h = [float('NaN')] * c_ion_time
        p_para_h = [float('NaN')] * c_ion_time
        p_perp_he = [float('NaN')] * c_ion_time
        p_para_he = [float('NaN')] * c_ion_time
        p_perp_o = [float('NaN')] * c_ion_time
        p_para_o = [float('NaN')] * c_ion_time
        p_perp_e = [float('NaN')] * c_ele_time
        p_para_e = [float('NaN')] * c_ele_time

        ang_perp = 0.5 * (numpy.sin(pa_arr) ** 3)
        ang_para = numpy.sin(pa_arr) * (numpy.cos(pa_arr) ** 2)

        # ion perp and para
        for it in tqdm(range(c_ion_time), desc='ele - perp and para'):
            temp_h = ion_data[it]['daty_avg_int_H1']
            temp_he = ion_data[it]['daty_avg_int_He1']
            temp_o = ion_data[it]['daty_avg_int_O1']

            for ie in range(len(ion_data[it]['daty_avg_int_H1'])):
                with contextlib.redirect_stdout(None):
                    h_perp[it][ie] = ((temp_h[ie] * ang_perp * del_pa)
                                      .sum(dim='PITCH_ANGLE', skipna=True)
                                      )
                    print(h_perp[it][ie])
                    h_para[it][ie] = ((temp_h[ie] * ang_para * del_pa)
                                      .sum(dim='PITCH_ANGLE', skipna=True)
                                      )
                    he_perp[it][ie] = ((temp_he[ie] * ang_perp * del_pa)
                                       .sum(dim='PITCH_ANGLE', skipna=True)
                                       )
                    he_para[it][ie] = ((temp_he[ie] * ang_para * del_pa)
                                       .sum(dim='PITCH_ANGLE', skipna=True)
                                       )
                    o_perp[it][ie] = ((temp_o[ie] * ang_perp * del_pa)
                                      .sum(dim='PITCH_ANGLE', skipna=True)
                                      )
                    o_para[it][ie] = ((temp_o[ie] * ang_para * del_pa)
                                      .sum(dim='PITCH_ANGLE', skipna=True)
                                      )

        for it in tqdm(range(c_ion_time)):
            if relat == 0:
                print(h_perp[it].dims)
                print(del_ion_en[it].dims)
                temp_h_perp = (numpy.sqrt(
                    2 * mass_pro * en_ion_erg[it]) * (
                        h_perp[it] * del_ion_en[it]))
                print(temp_h_perp.dims)
                temp_h_para = (numpy.sqrt(
                    2 * mass_pro * en_ion_erg[it]) * (
                        h_para[it] * del_ion_en[it]))
                temp_he_perp = (numpy.sqrt(
                    2 * mass_hel * en_ion_erg[it]) * (
                        he_perp[it] * del_ion_en[it]))
                temp_he_para = (numpy.sqrt(
                    2 * mass_hel * en_ion_erg[it]) * (
                        he_para[it] * del_ion_en[it]))
                temp_o_perp = (numpy.sqrt(
                    2 * mass_oxy * en_ion_erg[it]) * (
                        o_perp[it] * del_ion_en[it]))
                temp_o_para = (numpy.sqrt(
                    2 * mass_oxy * en_ion_erg[it]) * (
                        o_para[it] * del_ion_en[it]))

                '''p_perp_h[it] = 2 * pi * temp_h_perp.sum(dim='', skipna=True)
                p_para_h[it] = 2 * pi * temp_h_para.sum()
                p_perp_he[it] = 2 * pi * temp_he_perp.sum()
                p_para_he[it] = 2 * pi * temp_he_para.sum()
                p_perp_o[it] = 2 * pi * temp_o_perp.sum()
                p_para_o[it] = 2 * pi * temp_o_para.sum()
            elif relat == 1:'''

        print()

    def smooth_data(self, given_corrected_data, given_pressure):
        print()

    def average_pressures(self, given_pressure):
        print()

    def create_plot_vars(self, given_smooth_data, given_avg_pressure):
        print()

    def plot_data(self, given_plot_vars):
        print()

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

    main_calculations = HopeCalculations(time_s_str, time_e_str, probe,
                                         level, factor, low_energy, up_energy, 
                                         pot_corr, relat, swindow)

    datetime_list = main_calculations.datetime_range_list()
    hope_cdf_objs, efw_cdf_objs = main_calculations.get_cdf_objs(datetime_list)
    cdf_objs_data = main_calculations.read_cdf_data(hope_cdf_objs)

    # corrected_data = main_calculations.correct_fluxes(cdf_objs_data)
    pressure = main_calculations.calc_pressure(cdf_objs_data)


# Run main() if called from command line
if (__name__ == "__main__"):
    main()
