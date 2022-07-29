from calc_hope_pressure import HopeCalculations

import matplotlib as mpl
import matplotlib.pyplot as plt

from tqdm import tqdm

# matplotlib settings, feel free to change
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams.update({'font.size': 10})
mpl.rcParams["figure.figsize"] = [14.1, 9.7]
mpl.rcParams["figure.dpi"] = 226


def plot_data(given_cdf_data):
    test_pa = 5
    test_bin_1 = 20
    test_bin_2 = 40
    test_bin_3 = 60
    test_bin_4 = 72

    ion_data = given_cdf_data['ion_data_list']

    test_ion_data_1 = [
        ion_data[i]['energy'][test_bin_1 - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'Ion energy bin {test_bin_1}')]
    test_ion_data_2 = [
        ion_data[i]['energy'][test_bin_2 - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'Ion energy bin {test_bin_2}')]
    test_ion_data_3 = [
        ion_data[i]['energy'][test_bin_3 - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'Ion energy bin {test_bin_3}')]
    test_ion_data_4 = [
        ion_data[i]['energy'][test_bin_4 - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'Ion energy bin {test_bin_4}')]

    test_h1_data_1 = [
        ion_data[i]['daty_avg_int_H1'][test_bin_1 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'H1 bin {test_bin_1}')]
    print(f"h1 bin {test_bin_1} energy: "
          f"{test_h1_data_1[0]['HOPE_ENERGY_Ion'].values}")
    test_h1_data_2 = [
        ion_data[i]['daty_avg_int_H1'][test_bin_2 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'H1 bin {test_bin_2}')]
    print(f"h1 bin {test_bin_2} energy: "
          f"{test_h1_data_2[0]['HOPE_ENERGY_Ion'].values}")
    test_h1_data_3 = [
        ion_data[i]['daty_avg_int_H1'][test_bin_3 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'H1 bin {test_bin_3}')]
    print(f"h1 bin {test_bin_3} energy: "
          f"{test_h1_data_3[0]['HOPE_ENERGY_Ion'].values}")
    test_h1_data_4 = [
        ion_data[i]['daty_avg_int_H1'][test_bin_4 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'H1 bin {test_bin_4}')]
    print(f"h1 bin {test_bin_4} energy: "
          f"{test_h1_data_4[0]['HOPE_ENERGY_Ion'].values}")

    test_he1_data_1 = [
        ion_data[i]['daty_avg_int_He1'][test_bin_1 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'He1 bin {test_bin_1}')]
    print(f"he1 bin {test_bin_1} energy: "
          f"{test_he1_data_1[0]['HOPE_ENERGY_Ion'].values}")
    test_he1_data_2 = [
        ion_data[i]['daty_avg_int_He1'][test_bin_2 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'He1 bin {test_bin_2}')]
    print(f"he1 bin {test_bin_2} energy: "
          f"{test_he1_data_2[0]['HOPE_ENERGY_Ion'].values}")
    test_he1_data_3 = [
        ion_data[i]['daty_avg_int_He1'][test_bin_3 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'He1 bin {test_bin_3}')]
    print(f"he1 bin {test_bin_3} energy: "
          f"{test_he1_data_3[0]['HOPE_ENERGY_Ion'].values}")
    test_he1_data_4 = [
        ion_data[i]['daty_avg_int_He1'][test_bin_4 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'He1 bin {test_bin_4}')]
    print(f"he1 bin {test_bin_4} energy: "
          f"{test_he1_data_4[0]['HOPE_ENERGY_Ion'].values}")

    test_o1_data_1 = [
        ion_data[i]['daty_avg_int_O1'][test_bin_1 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'O1 bin {test_bin_1}')]
    print(f"o1 bin {test_bin_1} energy: "
          f"{test_o1_data_1[0]['HOPE_ENERGY_Ion'].values}")
    test_o1_data_2 = [
        ion_data[i]['daty_avg_int_O1'][test_bin_2 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'O1 bin {test_bin_2}')]
    print(f"o1 bin {test_bin_2} energy: "
          f"{test_o1_data_2[0]['HOPE_ENERGY_Ion'].values}")
    test_o1_data_3 = [
        ion_data[i]['daty_avg_int_O1'][test_bin_3 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'O1 bin {test_bin_3}')]
    print(f"o1 bin {test_bin_3} energy: "
          f"{test_o1_data_3[0]['HOPE_ENERGY_Ion'].values}")
    test_o1_data_4 = [
        ion_data[i]['daty_avg_int_O1'][test_bin_4 - 1][test_pa - 1]
        for i in tqdm(range(len(ion_data)), 
                      desc=f'O1 bin {test_bin_4}')]
    print(f"o1 bin {test_bin_4} energy: "
          f"{test_o1_data_4[0]['HOPE_ENERGY_Ion'].values}")

    test_dates = [ion_data[i]['epoch'] for i in tqdm(range(len(ion_data)), 
                                                     desc='getting dates')]

    ion_fig, ion_axs = plt.subplots(2, 2)
    h1_fig, h1_axs = plt.subplots(2, 2)
    he1_fig, he1_axs = plt.subplots(2, 2)
    o1_fig, o1_axs = plt.subplots(2, 2)

    ion_axs[0, 0].set_title('Ion energy, Bin 20')
    ion_axs[0, 0].plot(test_dates, test_ion_data_1)
    ion_axs[0, 1].set_title('Ion energy, Bin 40')
    ion_axs[0, 1].plot(test_dates, test_ion_data_2)
    ion_axs[1, 0].set_title('Ion energy, Bin 60')
    ion_axs[1, 0].plot(test_dates, test_ion_data_3)
    ion_axs[1, 1].set_title('Ion energy, Bin 72')
    ion_axs[1, 1].plot(test_dates, test_ion_data_4)

    h1_axs[0, 0].set_title('H1 PA 5, Energy Bin 20')
    h1_axs[0, 0].plot(test_dates, test_h1_data_1)
    h1_axs[0, 1].set_title('H1 PA 5, Energy Bin 40')
    h1_axs[0, 1].plot(test_dates, test_h1_data_2)
    h1_axs[1, 0].set_title('H1 PA 5, Energy Bin 60')
    h1_axs[1, 0].plot(test_dates, test_h1_data_3)
    h1_axs[1, 1].set_title('H1 PA 5, Energy Bin 72')
    h1_axs[1, 1].plot(test_dates, test_h1_data_4)

    he1_axs[0, 0].set_title('He1 PA 5, Energy Bin 20')
    he1_axs[0, 0].plot(test_dates, test_he1_data_1)
    he1_axs[0, 1].set_title('He1 PA 5, Energy Bin 40')
    he1_axs[0, 1].plot(test_dates, test_he1_data_2)
    he1_axs[1, 0].set_title('He1 PA 5, Energy Bin 60')
    he1_axs[1, 0].plot(test_dates, test_he1_data_3)
    he1_axs[1, 1].set_title('He1 PA 5, Energy Bin 72')
    he1_axs[1, 1].plot(test_dates, test_he1_data_4)

    o1_axs[0, 0].set_title('O1 PA 5, Energy Bin 20')
    o1_axs[0, 0].plot(test_dates, test_o1_data_1)
    o1_axs[0, 1].set_title('O1 PA 5, Energy Bin 40')
    o1_axs[0, 1].plot(test_dates, test_o1_data_2)
    o1_axs[1, 0].set_title('O1 PA 5, Energy Bin 60')
    o1_axs[1, 0].plot(test_dates, test_o1_data_3)
    o1_axs[1, 1].set_title('O1 PA 5, Energy Bin 72')
    o1_axs[1, 1].plot(test_dates, test_o1_data_4)

    for axs in [h1_axs, he1_axs, o1_axs]:
        for ax in axs.flat:
            ax.set_yscale('log')

    for ax in ion_axs.flat:
        ax.set(xlabel='time', ylabel='HOPE_ENERGY_Ion')

    for ax in h1_axs.flat:
        ax.set(xlabel='time', ylabel='daty_avg_int_H1')

    for ax in he1_axs.flat:
        ax.set(xlabel='time', ylabel='daty_avg_int_He1')

    for ax in o1_axs.flat:
        ax.set(xlabel='time', ylabel='daty_avg_int_O1')  

    ion_fig.savefig('figout_ion.svg', format='svg')
    ion_fig.savefig('figout_ion.png', format='png')
    h1_fig.savefig('figout_H1.svg', format='svg')
    h1_fig.savefig('figout_H1.png', format='png')
    he1_fig.savefig('figout_He1.svg', format='svg')
    he1_fig.savefig('figout_He1.png', format='png')
    o1_fig.savefig('figout_O1.svg', format='svg')
    o1_fig.savefig('figout_O1.png', format='png')


def main():
    """
    Main function to run the program via command line. 
    Functionality included to parse and validate command line arguments so 
    modifications to the program are not required.
    """    

    time_s_str = "2014-08-26/00:00:00"
    time_e_str = "2014-08-28/06:00:00"
    level = 'l3'
    probe = 'a'
    factor = 1
    low_energy = [1e2, 1e2]
    up_energy = [5.5e4, 5.5e4]
    pot_corr = 1
    relat = 1
    swindow = 60

    main_calculations = HopeCalculations(time_s_str, time_e_str, probe,
                                         level, factor, low_energy, up_energy, 
                                         pot_corr, relat, swindow)
    datetime_list = main_calculations.datetime_range_list()
    hope_cdf_objs, efw_cdf_objs = main_calculations.get_cdf_objs(datetime_list)
    cdf_objs_data = main_calculations.read_cdf_data(hope_cdf_objs)
    plot_data(cdf_objs_data)


if (__name__ == "__main__"):
    main()
