import pandas as pd
import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import HourLocator, DateFormatter, DayLocator



def GetIrradianceTslist(filepath, filename):
    df = pd.read_csv(os.path.join(filepath, filename), header=None, sep="\s+", skiprows=1)

    interval = 0.25
    num_columns = 5

    #TS_List_Processing
    df['groups'] = df.iloc[:, 1] // interval
    grps = df.groupby(['groups'])

    #for i in range(num_columns):
    #    print(i)
    swdwn2 = []
    swdif2 = []
    swdni2 = []
    for groups in grps:
#        print(groups)
        data = groups[1]
        # print(data.iloc[:, -6])
        swdwn2.append(data.iloc[:, -6].mean())
        swdif2.append(data.iloc[:, -4].mean())
        swdni2.append(data.iloc[:, -5].mean())
        # swdwn.append(data.iloc[:, -11].mean())
    return swdwn2, swdni2, swdif2

def GetObservedIrradiance(mask1, mask):
    swdown = []
    swdni = []
    swdif = []
    for i in range(len(mask1)):
#    print(mask1[i])
        for j in range(len(mask.time.values)):
#        print(mask.time.values[j])
            if mask1[i] == mask.time.values[j]:
                print(mask.time.values[j], mask.SWdown_global.values[j], i)
                swdown.append(mask.SWdown_global.values[j])
                swdni.append(mask.SWdown_direct.values[j])
                swdif.append(mask.SWdown_diffuse.values[j])

    # Need to figure out a cleaner way to put nans in if a vaule is missing            
    # values.insert(27,np.nan) 
    return swdown, swdni, swdif

##################################################
####       Progam Begins Here                 ####  
##################################################

# /home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/12Z/AOD_2_Ang_0_0034/CLDMASK/Ts_list

BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/'
TS_DIR = BASE_DIR + '18Z/CLDMASK/Ts_List/'
TS_DIR_CLDTOP = BASE_DIR + '18Z/CLDTOPZ_CLDBASEZ/Ts_List/'
TS_DIR_BR = BASE_DIR + '18Z/BRTEMP_CLDMASK_CLDBASEZ/Ts_List/'
ts_file = 'Bco.d04.TS'
# 
# print(TS_DIR[-9:-1])
PYRA_DIR = BASE_DIR + 'Observed_Data/'
pyr_file = 'Radiation__Deebles_Point__DownwellingRadiation__1s__20210616.nc'

PNG = BASE_DIR + '18Z/PNG/'

swdwn2, swdni2, swdif2 = GetIrradianceTslist(TS_DIR, ts_file )
swdwnctop, swdni2ctop, swdif2ctop = GetIrradianceTslist(TS_DIR_CLDTOP, ts_file )
swdwnbr, swdni2br, swdif2br = GetIrradianceTslist(TS_DIR_BR, ts_file )

# ds = xr.open_dataset(PYRA_DIR + pyr_file)
ds = xr.open_mfdataset(PYRA_DIR + 'Radiation__*.nc', combine='by_coords')
mask1 = pd.date_range("2021-06-16 18:00:00", freq="15T", periods=25)
mask = ds.sel(time=slice('2021-06-16 18:00:00', '2021-06-17 00:00:00'))
# print(len(mask1), len(swdwn2))
obs_swdwn, obs_swdni, obs_swdif = GetObservedIrradiance(mask1, mask)

# print(len(mask1), len(obs_swdwn))
# print(len(mask1), len(swdwn2))

d_fmt = DateFormatter("%m-%d")

fig = plt.figure(figsize=(12,10))
ax = plt.subplot(3,1,1)
# plt.title()
plt.plot(mask1, swdwn2, color='b', label='cldmask', linestyle='--', marker='*')
plt.plot(mask1, swdwnctop, color='g', label='ctoph', linestyle='--', marker='*')
plt.plot(mask1, swdwnbr, color='c', label='brtemp', linestyle='--', marker='*')
plt.plot(mask1, obs_swdwn, color='r',label='obs', linestyle='-', marker='*')
# plt.text(mask1[-8], max(swdwn2) - 100, 'AOD = 2',color='k',style='italic')
# plt.text(mask1[-8], max(swdwn2) - 150, 'Ang_Exp = 0.034' ,color='k', style='italic')
# plt.plot(mask1, swdwn2, color='g',label='swdwn2', linestyle=':', marker='*')
# plt.xlabel('Time (HH:MM)', fontsize=15)
plt.ylabel('Global Horizontal \n Irradiance $W/{m}^2$', fontsize=15)
# plt.xaxis.set_major_formatter(d_fmt)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=15)
plt.yticks(fontweight='semibold', fontsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best')

ax2 = plt.subplot(3,1,2)
plt.plot(mask1, swdni2, color='b', label='cldmask', linestyle='--', marker='*')
plt.plot(mask1, swdni2ctop, color='g', label='ctoph', linestyle='--', marker='*')
plt.plot(mask1, swdni2br, color='c', label='brtemp', linestyle='--', marker='*')
plt.plot(mask1, obs_swdni, color='r',label='obs', linestyle='-', marker='*')
# plt.text(mask1[-8], max(swdwn2) - 100, 'AOD = 2',color='k',style='italic')
# plt.text(mask1[-8], max(swdwn2) - 150, 'Ang_Exp = 0.034' ,color='k', style='italic')
# plt.plot(mask1, swdwn2, color='g',label='swdwn2', linestyle=':', marker='*')
# plt.xlabel('Time (HH:MM)', fontsize=15)
plt.ylabel('Direct Normal \n Irradiance $W/{m}^2$', fontsize=15)
# plt.xaxis.set_major_formatter(d_fmt)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=15)
plt.yticks(fontweight='semibold', fontsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best')

ax2 = plt.subplot(3,1,3)
plt.plot(mask1, swdif2, color='b', label='cldmask', linestyle='--', marker='*')
plt.plot(mask1, swdif2ctop, color='g', label='ctoph', linestyle='--', marker='*')
plt.plot(mask1, swdif2br, color='c', label='brtemp', linestyle='--', marker='*')
plt.plot(mask1, obs_swdif, color='r',label='obs', linestyle='-', marker='*')
# plt.text(mask1[-8], max(swdwn2) - 100, 'AOD = 2',color='k',style='italic')
# plt.text(mask1[-8], max(swdwn2) - 150, 'Ang_Exp = 0.034' ,color='k', style='italic')
# plt.plot(mask1, swdwn2, color='g',label='swdwn2', linestyle=':', marker='*')
plt.xlabel('Time (HH:MM)', fontsize=15)
plt.ylabel('Diffused \n Irradiance $W/{m}^2$', fontsize=15)
# plt.xaxis.set_major_formatter(d_fmt)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=15)
plt.yticks(fontweight='semibold', fontsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best')



plt.savefig(PNG + "BCO_18Z_Run.png", dpi=300, facecolor='w', 
            edgecolor='w', orientation='lanscape', papertype=None, format='png',
            bbox_inches='tight', pad_inches=0.1)

plt.show()