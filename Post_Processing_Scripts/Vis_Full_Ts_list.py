import pandas as pd
import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates





BASE_DIR = '/home/alton/WRF_OUT/New_Experiments/Experiment6/CLDBASEZ_Interp_Nearest/20210616/12Z/AOD_0.616_Ang-0.054/'
TS_DIR_CLDMASK = BASE_DIR + 'CLDMASK_BRTEMP/Ts_List/'
TS_DIR_CLDTOPZ = BASE_DIR + 'CLDTOPZ_CLDBASEZ/Ts_List/'

filename = 'Cimh.d04.TS'
filename2 = 'Cimh.d04.TS'

df = pd.read_csv(os.path.join(TS_DIR_CLDMASK, filename), header=None, sep="\s+", skiprows=1)

df.iloc[:, -5].plot()

plt.show()

cmask = df.iloc[:, -6]

print(cmask)

