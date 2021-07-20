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

##################################################
####          Program begins here             ####  
##################################################

BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/12Z/AOD_0_091_Ang_1.323/'
TS_DIR_CLDMASK = BASE_DIR + 'CLDMASK/Ts_list/'
TS_DIR_CLDTOPZ = BASE_DIR + 'CLDTOPZ_CLDBASEZ/Ts_List/'
TS_DIR_BRTEMP = BASE_DIR + 'BRTEMP_CLDMASK_CLDBASEZ/Ts_List/'

PNG = BASE_DIR + 'Summary/'
ts_file = 'Cimh.d04.TS'

PATH = '/home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/Observed_Data/'
file = 'Solar_Request.xlsx'

#Extracting data from tslist files
swdwn2_cldmask = GetIrradianceTslist(TS_DIR_CLDMASK, ts_file)
swdwn2_cldtopz = GetIrradianceTslist(TS_DIR_CLDTOPZ, ts_file)
swdwn2_brtemp = GetIrradianceTslist(TS_DIR_BRTEMP, ts_file)

#Reading Excel File
obs = pd.read_excel(PATH+file, 
                    sheet_name='Sheet2', 
                    parse_dates=[['Date','Time']])

#Creating Time Periods of Interest
mask1 = pd.date_range("2021-06-16 12:00:00", freq="15T", periods=25)
# print(len(obs['Average W/m2'][26:75]),len(mask1))


fig = plt.figure(figsize=(10,5))
ax = plt.subplot(1,1,1)
# plt.title()
plt.plot(mask1, swdwn2_cldmask, color='b', label='cldmask', linestyle='--', marker='*')
plt.plot(mask1, swdwn2_cldtopz, color='g', label='ctoph', linestyle='--', marker='*')
plt.plot(mask1, swdwn2_brtemp, color='c', label='brtemp', linestyle='--', marker='*')
plt.plot(mask1, obs['Average W/m2'][26:51], color='r',label='obs', linestyle='-', marker='*')
plt.ylabel('Global Horizontal \n Irradiance $W/{m}^2$', fontsize=15)
# plt.xaxis.set_major_formatter(d_fmt)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=15)
plt.yticks(fontweight='semibold', fontsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best')

plt.savefig(PNG + "BCO_12Z_Run.png", dpi=300, facecolor='w', 
            edgecolor='w', orientation='lanscape', papertype=None, format='png',
            bbox_inches='tight', pad_inches=0.1)

plt.show()