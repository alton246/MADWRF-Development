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

BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/CLDMASK_BRTEMP_CLDBASEZ/'
TS_DIR = BASE_DIR + 'Ts_List/'

# BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/Experiment5/'
# TS_DIR = BASE_DIR + 'CLDBASEZ_Interp_Nearest/Run_with_CLDTOPZ_CLDBASEZ/Ts_list/'
ts_file = 'Cimh.d04.TS'

PATH = '/home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/Observed_Data/'
file = 'Solar_Request.xlsx'

PNG = TS_DIR + '/PNG/'

swdwn2 = GetIrradianceTslist(TS_DIR, ts_file)

#Reading Excel File
obs = pd.read_excel(PATH+file, 
                    sheet_name='Sheet2', 
                    parse_dates=[['Date','Time']])

#Creating Time Periods of Interest
mask1 = pd.date_range("2020-06-22 06:00:00", freq="15T", periods=49)
# print(len(obs['Average W/m2'][26:75]),len(mask1))
# print(len(swdwn2),len(mask1))

#Plotting Data
fig = plt.figure(figsize=(12,5))
#ax = fig.gca()
#ax.set_xticks(np.arange(0, 48, 1))
plt.plot(mask1, obs['Average W/m2'][0:49], color='r',label='ghi_obs', linestyle='-', marker='*')
plt.plot(mask1, swdwn2, color='g',label='swdwn2', linestyle=':', marker='*')
plt.text(mask1[0], max(swdwn2) - 100, 'AOD = 2',color='k',style='italic')
# plt.text(mask1[0], max(swdwn2) - 150, 'Ang_Exp = ' + TS_DIR[-6:-1]  ,color='k',style='italic')
plt.text(mask1[0], max(swdwn2) - 150, 'Ang_Exp = 0.034' ,color='k',style='italic')
plt.xlabel('Time (HH:MM)', fontsize=15)
plt.ylabel('Irradiance $W/{m}^2$', fontsize=15)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.xticks(fontweight='semibold', fontsize=15)
plt.yticks(fontweight='semibold', fontsize=15)
plt.grid(True, lw=0.5, ls=':')
plt.legend(loc='best')

plt.savefig(PNG + "CIMH_06Z_Run.png", dpi=300, facecolor='w', 
            edgecolor='w', orientation='lanscape', papertype=None, format='png',
            bbox_inches='tight', pad_inches=0.1)

plt.show()