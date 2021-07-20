import pandas as pd
import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


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

def GetSwdni2Tslist(filepath, filename):
        df = pd.read_csv(os.path.join(filepath, filename), header=None, sep="\s+", skiprows=1)

        interval = 0.25
        num_columns = 5

        #TS_List_Processing
        df['groups'] = df.iloc[:, 1] // interval
        grps = df.groupby(['groups'])

        #for i in range(num_columns):
        #    print(i)
        swdni2 = []
        for groups in grps:
    #        print(groups)
            data = groups[1]
            swdni2.append(data.iloc[:, -5].mean())
            # swdwn.append(data.iloc[:, -11].mean())
        return swdni2

def GetSwdif2Tslist(filepath, filename):
    df = pd.read_csv(os.path.join(filepath, filename), header=None, sep="\s+", skiprows=1)

    interval = 0.25
    num_columns = 5

    #TS_List_Processing
    df['groups'] = df.iloc[:, 1] // interval
    grps = df.groupby(['groups'])

    #for i in range(num_columns):
    #    print(i)
    swdif2 = []
    swdwn = []
    for groups in grps:
#        print(groups)
        data = groups[1]
        swdif2.append(data.iloc[:, -4].mean())
        # swdwn.append(data.iloc[:, -11].mean())
    return swdif2

def GetObservedIrradiance(mask1, mask):
    values = []
    for i in range(len(mask1)):
        # print(mask1[i])
        for j in range(len(mask.time.values)):
            # print(mask.time.values[j])
            if mask1[i] == mask.time.values[j]:
                # if mask1[i] == None:
                    # print(mask1[i])
                # print(mask.time.values[j], mask.SWdown_global.values[j], i)
                values.append(mask.SWdown_global.values[j])
# values.insert(3,np.nan)
    return values

##################################################
####                                          ####  
##################################################

BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/Experiment5/'
# TS_DIR_CUB = BASE_DIR + 'CLDBASEZ_Interp_Nearest/20210615/CLDMASK/Ts_List'
TS_DIR_NE =  BASE_DIR + 'CLDBASEZ_Interp_Nearest/20210616/06Z/AOD_0_12_Ang0_26/CLDMASK/Ts_list/'
ts_file = 'Bco.d04.TS'

# print(TS_DIR[-9:-1])
PYRA_DIR = BASE_DIR + 'CLDBASEZ_Interp_Nearest/20210616/Observed_Data/'
pyr_file = 'Radiation__Deebles_Point__DownwellingRadiation__1s__20210616.nc'

PNG = TS_DIR_NE + 'PNG/'
# swdwn2_cub = GetIrradianceTslist(TS_DIR_CUB, ts_file )
swdwn2_ne = GetIrradianceTslist(TS_DIR_NE, ts_file )


# irr = xr.open_dataset(PATH + file)
# ds = xr.open_mfdataset(PYRA_DIR + 'Radiation__*.nc', concat_dim="time")
ds = xr.open_dataset(PYRA_DIR + pyr_file)
mask1 = pd.date_range("2021-06-16 06:00:00", freq="15T", periods=49)
mask = ds.sel(time=slice('2021-06-16 06:00:00', '2021-06-16 18:00:00'))
# print(mask)
obs_swdwn = GetObservedIrradiance(mask1, mask)

obs_swdwn.insert(1,np.nan)
# print(len(swdwn2), len(obs_swdwn))

fig = plt.figure(figsize=(12,5))
#ax = fig.gca()
#ax.set_xticks(np.arange(0, 48, 1))
# plt.plot(mask1, swdwn2_cub, color='b', label='swdwn_cub', linestyle='--', marker='*')
plt.plot(mask1, swdwn2_ne, color='b', label='swdwn_near', linestyle='--', marker='*')
plt.plot(mask1, obs_swdwn, color='r',label='ghi_obs', linestyle='-', marker='*')
plt.text(mask1[0], max(swdwn2_ne) - 100, 'AOD = 0.12',color='k',style='italic')
# plt.text(mask1[0], max(swdwn2) - 150, 'Ang_Exp = ' + TS_DIR[-6:-1]  ,color='k',style='italic')
plt.text(mask1[0], max(swdwn2_ne) - 150, 'Ang_Exp = 0.26' ,color='k',style='italic')
# plt.plot(mask1, swdwn2, color='g',label='swdwn2', linestyle=':', marker='*')
plt.xlabel('Time (HH:MM)', fontsize=15)
plt.ylabel('Irradiance $W/{m}^2$', fontsize=15)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=15)
plt.yticks(fontweight='semibold', fontsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best')

# plt.savefig(PNG + "BCO_06Z_" + TS_DIR[-9:-1] + "_Run.png", dpi=300, facecolor='w',
plt.savefig(PNG + "BCO_06Z_Run_CLDMASK.png", dpi=300, facecolor='w', 
            edgecolor='w', orientation='lanscape', papertype=None, format='png',
            bbox_inches='tight', pad_inches=0.1)

plt.show()

