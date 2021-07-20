import pandas as pd
import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import HourLocator, DateFormatter, DayLocator

def GetSwdwnTslist(filepath, filename):
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
                # print(mask.time.values[j], mask.SWdown_global.values[j], i)
                swdown.append(mask.SWdown_global.values[j])
                swdni.append(mask.SWdown_direct.values[j])
                swdif.append(mask.SWdown_diffuse.values[j])

    # Need to figure out a cleaner way to put nans in if a vaule is missing            
    # values.insert(27,np.nan) 
    return swdown, swdni, swdif

##################################################
####       Progam Begins Here                 ####  
##################################################


BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/'
PNG = '/home/alton/Github/MADWRF-Development/MADWRF-Research_Paper/Figure/'

#Paths for tsfiles
TS_DIR = BASE_DIR + '12Z/AOD_0_091_Ang_1.323/CLDMASK/Ts_list/'
TS_DIR_CLDTOP = BASE_DIR + '12Z/AOD_0_091_Ang_1.323/CLDTOPZ_CLDBASEZ/Ts_List/'
TS_DIR_BR = BASE_DIR + '12Z/AOD_0_091_Ang_1.323/BRTEMP_CLDMASK_CLDBASEZ/Ts_List/'
bco_ts_file = 'Bco.d04.TS'
cimh_ts_file = 'Cimh.d04.TS'

#PATH to Observed Data
OBS_DIR = BASE_DIR + 'Observed_Data/'
bco_pyr_file = 'Radiation__Deebles_Point__DownwellingRadiation__1s__20210616.nc'
cimh_pyr_file = 'Solar_Request.xlsx'

ds = xr.open_dataset(OBS_DIR + bco_pyr_file)
cimh_obs = pd.read_excel(OBS_DIR + cimh_pyr_file, 
                    sheet_name='Sheet2', 
                    parse_dates=[['Date','Time']])
mask1 = pd.date_range("2021-06-16 12:00:00", freq="15T", periods=25)
mask = ds.sel(time=slice('2021-06-16 12:00:00', '2021-06-16 18:00:00'))

bco_swdwn, bco_swdni, bco_swdif = GetObservedIrradiance(mask1, mask)

#Retrieve Data from BCO tslist
swdwn2, swdni2, swdif2 = GetIrradianceTslist(TS_DIR, bco_ts_file )
swdwnctop, swdni2ctop, swdif2ctop = GetIrradianceTslist(TS_DIR_CLDTOP, bco_ts_file )
swdwnbr, swdni2br, swdif2br = GetIrradianceTslist(TS_DIR_BR, bco_ts_file )

swdwn2_cldmask = GetSwdwnTslist(TS_DIR, cimh_ts_file)
swdwn2_cldtopz = GetSwdwnTslist(TS_DIR_CLDTOP, cimh_ts_file)
swdwn2_brtemp = GetSwdwnTslist(TS_DIR_BR, cimh_ts_file)

d_fmt = DateFormatter("%m-%d")
plt.rcParams['font.weight']='semibold'

fig = plt.figure(figsize=(10,5))
ax = plt.subplot(2,1,1)
# plt.title()
plt.plot(mask1, swdwn2, color='b', label='cldmask', linestyle='--', marker='*')
plt.plot(mask1, swdwnctop, color='g', label='ctoph', linestyle='--', marker='*')
plt.plot(mask1, swdwnbr, color='c', label='brtemp', linestyle='--', marker='*')
plt.plot(mask1, bco_swdwn, color='r',label='obs', linestyle='-', marker='*')
plt.text(mask1[0], 900, 'a)', color='k', style='normal',fontsize='12')
plt.ylabel('Global Horizontal \n Irradiance $W/{m}^2$', fontsize=15, fontweight='semibold')
# plt.xaxis.set_major_formatter(d_fmt)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=12)
plt.yticks(fontweight='semibold', fontsize=12)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best')

ax2 = plt.subplot(2,1,2)
# plt.title()
plt.plot(mask1, swdwn2_cldmask, color='b', label='cldmask', linestyle='--', marker='*')
plt.plot(mask1, swdwn2_cldtopz, color='g', label='ctoph', linestyle='--', marker='*')
plt.plot(mask1, swdwn2_brtemp, color='c', label='brtemp', linestyle='--', marker='*')
plt.plot(mask1, cimh_obs['Average W/m2'][26:51], color='r',label='obs', linestyle='-', marker='*')
plt.text(mask1[0], 900, 'b)', color='k', style='normal',fontsize='12')
plt.ylabel('Global Horizontal \n Irradiance $W/{m}^2$', fontsize=15, fontweight='semibold')
plt.xlabel('Date (hh:mm)', fontsize=15, fontweight='semibold')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=12)
plt.yticks(fontweight='semibold', fontsize=12)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best')

plt.savefig(PNG + "BCO_CIMH_12Z_6HR_Run.png", dpi=300, facecolor='w', 
            edgecolor='w', orientation='lanscape', papertype=None, format='png',
            bbox_inches='tight', pad_inches=0.1)

plt.show()
