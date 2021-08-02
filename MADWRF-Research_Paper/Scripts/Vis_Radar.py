from datetime import datetime
from datetime import datetime
from datetime import timedelta
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import bz2
import pandas as pd
import numpy as np
import glob


PATH = '/home/alton/Github/MADWRF-Development/MADWRF-Research_Paper/Data/BCo_Radar/'
nc = Dataset(PATH + 'MMCR__MBR2__Spectral_Moments__10s__155m-18km__200622.nc')

date = nc.variables['YYYYMMDDHHMM'][:]
time = nc.variables['time'][:]
ref = nc.variables['Ze'][:]
hgt = nc.variables['range'][:]
# snr = nc.variables['SNRplank'][:]
# print(nc.variables.keys())

# print(nc.variables['Ze'])

DateTime = []
for i in range(len(time)):
    DateTime.append(datetime.strptime('20200622'+ ' ' + 
                                      str(timedelta(seconds=np.float64(time[i])))[12:20], 
                                      '%Y%m%d %H:%M:%S'))


plt.rcParams['font.weight']='semibold'
plt.rcParams['font.size']='9'
fig= plt.subplots(figsize=(10,4))
ax =plt.subplot(1,1,1)
plt.contourf(DateTime,hgt,ref.transpose(), c='jet')
plt.ylim(0,2500)
plt.xlim('2020-06-22 12:00','2020-06-22 18:00:00')
# ax.axvspan("2021-06-16 15:28:00", "2021-06-16 16:02:00", color='grey', alpha=0.3)
plt.xlabel('UTC Time (mm:ss)',fontsize=9, fontweight='semibold')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.ylabel('Height (m)', fontsize=9, fontweight='semibold')
plt.tick_params(axis='both', which='major', labelsize=9)
plt.tick_params(axis='both', which='minor', labelsize=9)
plt.grid(True, lw=1, ls=':')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Radar Reflectivity \n Factor Ze (dBZ)', fontweight='semibold',fontsize=12)


plt.savefig('/home/alton/Github/MADWRF-Development/MADWRF-Research_Paper/Figure/Radar_Ref_20200622.png', dpi=300, facecolor='w', edgecolor='w',
           orientation='lanscape', papertype=None, format='png',
             bbox_inches='tight', pad_inches=0.1)
plt.show()