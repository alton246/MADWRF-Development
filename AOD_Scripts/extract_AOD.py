import xarray as xr
import numpy as np
import pandas as pd
import glob
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from netCDF4 import Dataset
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx

def AODExtract(PATH, lat, lon):
    aod = np.array([])
    for file in sorted(glob.glob(PATH + "wrfout_d01*")):
        # print(file)
        ncfile = Dataset(str(file))
        slp_1 = ncfile.variables['aod'][0][:]
        # print(ncfile.variables)
        # slp_1 = getvar(ncfile, "aod")
        lons_1 = ncfile.variables['longitude'][0, :]
        lats_1 = ncfile.variables['latitude'][:, 0]

        # # Smooth the sea level pressure since it tends to be noisy near the
        # # mountains
        smooth_slp_1 = smooth2d(slp_1, 3, cenweight=4)

        # Get the latitude and Longitude
        lati = lat
        long = lon
        y = find_nearest(lons_1,long)
        x = find_nearest(lats_1,lati)
        # print(smooth_slp_1[x,y].values)
        aod = np.append(aod,smooth_slp_1[x,y].values)
    return aod

PATH = '/home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/AOD/WRF_CHEM-16062021-17062021/16/1200/'
file = 'wrfout_d01_2021-06-16_15:00:00.nc'

PNG = '/home/alton/Github/MADWRF-Development/AOD_Scripts/Plots/'
# ncfile = Dataset(PATH + file)
# # print(ncfile.variables)
# print(ncfile.variables['aod'][0])
# slp_1 = getvar(ncfile, "aod")
# lons_1 = ncfile.variables['longitude'][0, :]
# # print(lons_1)
# lats_1 = ncfile.variables['latitude'][:, 0]
# print(lats_1)
legend_properties = {'weight':'semibold','size':'7'}
AOD_RP= AODExtract(PATH, 13.165, -59.432)
AOD_Gu= AODExtract(PATH, 16.225, -61.528)

mask1 = pd.date_range("2021-06-16 12:00:00", freq="3H", periods=57)

fig = plt.figure(figsize=(10,4))

ax2 = plt.subplot(1, 1, 1)
# plt.title()
plt.plot(mask1, AOD_RP, color='b', label='Ragged_Point', linestyle='--', marker='*')
plt.plot(mask1, AOD_Gu, color='g', label='Guadeloupe', linestyle='--', marker='*')
# plt.plot(mask1, swdwn2_brtemp, color='c', label='brtemp', linestyle='--', marker='*')
# plt.plot(mask1, cimh_obs['Average W/m2'][48:73], color='r',label='obs', linestyle='-', marker='*')
# ax2.axvspan("2021-06-16 15:28:00", "2021-06-16 16:02:00", color='grey', alpha=0.3)
# plt.text(mask1[0], 200, 'b)', color='k', style='normal',fontsize='9')
plt.ylabel('AOD', fontsize=9, fontweight='semibold')
plt.xlabel('Date (hh:mm)', fontsize=9, fontweight='semibold')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=9)
plt.yticks(fontweight='semibold', fontsize=9)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best', prop=legend_properties)

plt.savefig(PNG + "AOD_Model_20210616.png", dpi=300, facecolor='w', 
            edgecolor='w', orientation='lanscape', papertype=None, format='png',
            bbox_inches='tight', pad_inches=0.1)

plt.show()