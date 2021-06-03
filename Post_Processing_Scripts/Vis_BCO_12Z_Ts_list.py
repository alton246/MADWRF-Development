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
    swdwn = []
    for groups in grps:
#        print(groups)
        data = groups[1]
        swdwn2.append(data.iloc[:, -6].mean())
        # swdwn.append(data.iloc[:, -11].mean())
    return swdwn2

def GetObservedIrradiance(mask1, mask):
    values = []
    for i in range(len(mask1)):
#    print(mask1[i])
        for j in range(len(mask.time.values)):
#        print(mask.time.values[j])
            if mask1[i] == mask.time.values[j]:
                print(mask.time.values[j], mask.SWdown_global.values[j], i)
                values.append(mask.SWdown_global.values[j])

    # Need to figure out a cleaner way to put nans in if a vaule is missing            
    values.insert(27,np.nan) 
    return values

##################################################
####       Progam Begins Here                 ####  
##################################################

BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/20200622/'
TS_DIR = BASE_DIR + '12Z/Ts_List_Ang0_068/'
ts_file = 'Bco.d04.TS'

print(TS_DIR[-9:-1])
PYRA_DIR = BASE_DIR + 'Pyranometer_Data/'
pyr_file = 'Radiation__Deebles_Point__DownwellingRadiation__1s__20200622.nc'

PNG = TS_DIR + 'PNG/'
swdwn2 = GetIrradianceTslist(TS_DIR, ts_file )


# irr = xr.open_dataset(PATH + file)
ds = xr.open_mfdataset(PYRA_DIR + 'Radiation__*.nc', concat_dim="time")
mask1 = pd.date_range("2020-06-22 12:00:00", freq="15T", periods=49)
mask = ds.sel(time=slice('2020-06-22 12:00:00', '2020-06-23 00:00:00'))

obs_swdwn = GetObservedIrradiance(mask1, mask)

print(len(swdwn2), len(obs_swdwn))

d_fmt = DateFormatter("%m-%d")

fig = plt.figure(figsize=(12,5))
#ax = fig.gca()
#ax.set_xticks(np.arange(0, 48, 1))
plt.plot(mask1, swdwn2, color='b', label='swdwn', linestyle='--', marker='*')
plt.plot(mask1, obs_swdwn, color='r',label='ghi_obs', linestyle='-', marker='*')
plt.text(mask1[-8], max(swdwn2) - 100, 'AOD = 2',color='k',style='italic')
plt.text(mask1[-8], max(swdwn2) - 150, 'Ang_Exp = ' + TS_DIR[-6:-1]  ,color='k',style='italic')
# plt.plot(mask1, swdwn2, color='g',label='swdwn2', linestyle=':', marker='*')
plt.xlabel('Time (HH:MM)', fontsize=15)
plt.ylabel('Irradiance $W/{m}^2$', fontsize=15)
# plt.xaxis.set_major_formatter(d_fmt)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=15)
plt.yticks(fontweight='semibold', fontsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best')

plt.savefig(PNG + "BCO_06Z_" + TS_DIR[-9:-1] + "_Run.png", dpi=300, facecolor='w', 
            edgecolor='w', orientation='lanscape', papertype=None, format='png',
            bbox_inches='tight', pad_inches=0.1)

plt.show()