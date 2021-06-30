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

# BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/'
# TS_DIR = BASE_DIR + '12Z/AOD_2_Ang_0_0034/CLDMASK/Ts_List/'

# BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/Experiment5/'
# TS_DIR = BASE_DIR + 'CLDBASEZ_Interp_Nearest/Run_with_CLDTOPZ_CLDBASEZ/Ts_list/'
# ts_file = 'Cimh.d04.TS'

PATH = '/home/alton/WRF_OUT/New_Experiments/Experiment5/CLDBASEZ_Interp_Nearest/20210616/Observed_Data/'
file = 'Solar_Request.xlsx'

# PNG = TS_DIR + '/PNG/'

# swdwn2 = GetIrradianceTslist(TS_DIR, ts_file)

#Reading Excel File
obs = pd.read_excel(PATH+file, 
                    sheet_name='Sheet2', 
                    parse_dates=[['Date','Time']])

#Creating Time Periods of Interest
mask1 = pd.date_range("2021-06-16 12:00:00", freq="15T", periods=49)
print(len(obs['Average W/m2'][26:51]),len(mask1))