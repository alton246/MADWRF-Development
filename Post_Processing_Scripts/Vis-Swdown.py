from numpy.core.fromnumeric import size
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics
import cartopy.crs as ccrs                   # import projections
import cartopy.feature as cf 
import cmocean
import cmocean.cm as cmo 
import os               # import features
from glob import glob
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset







#############################################
###                                       ###
#############################################

PATH = '/home/alton/WRF_OUT/New_Experiments/Experiment6/CLDBASEZ_Interp_Nearest/20210616/12Z/AOD_0.616_Ang-0.054/CLDTOPZ_CLDBASEZ/Output_Files/'
file = 'wrfout_ghi_d04_2021-06-16_12:00:00'
PNG = PATH + 'PNG/'
# PNG2 = '/home/alton/WRF_OUT/New_Experiments/20200622/06Z/Output_Files_Ang0_034/Plots/'

plt.rcParams['font.weight']='semibold'
plt.rcParams['font.size']='14'
irr = Dataset(PATH+file, 'r')

mask1 = pd.date_range("2021-06-16 12:00:00", freq="15T", periods=49)
# print(mask1[0])
for i in range(irr['Times'].shape[0]):
    print(mask1[i])
    # print(irr['Times'].shape[0])

    x1,x2=irr['XLONG'][i].data.min(),irr['XLONG'][i].data.max()
    y1,y2=irr['XLAT'][i].data.min(),irr['XLAT'][i].data.max()
    # print(irr['XLAT'][0].data.min())
    # print(x1,y1)
    # print(x2,y2)

    fig = plt.figure(figsize=(18,10))
    
    ax1 = plt.subplot(111,projection = ccrs.PlateCarree())

    v = np.linspace(0, 1100, 23, endpoint=True)
    # v = np.linspace(0, 1000, 21, endpoint=True)

    # swh = plt.contourf(irr['XLONG'][i][0,:], irr['XLAT'][i][:,0], irr['SWDOWN2'][i], v,
                #   transform=ccrs.PlateCarree(),cmap='gist_ncar')
    swh = plt.contourf(irr['XLONG'][i][0,:], irr['XLAT'][i][:,0], irr['SWDOWN2'][i], v,
                  transform=ccrs.PlateCarree(),cmap=plt.cm.Greys)

    plt.title(mask1[i], fontsize='18', fontweight='semibold')

    states = NaturalEarthFeature(category="cultural", scale="10m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    
    ax1.add_feature(states, linewidth=0.5, edgecolor="black")
    
    ax1.coastlines(resolution='10m')                
    ax1.set_extent([x1,x2,y1,y2])

    ax = plt.scatter(-59.4289, 13.1627, marker='*', color='k', s=100)
    ax2 = plt.scatter(-59.6245, 13.1499, marker='o', color='k', s=100)
    
    # cb = plt.colorbar(swh, orientation='horizontal', ticks=v)
    cbar_ax = fig.add_axes([0.09, 0.02, 0.84, 0.02])
    cb = fig.colorbar(swh, cax=cbar_ax, orientation='horizontal', ticks=v)
    cb.ax.set_title('Solar Irradiance ($W/{m}^2$)', fontweight='semibold', 
                    fontsize=14)

    gl = ax1.gridlines(draw_labels=True,
            linewidth=0.5, color='black', alpha=0.5, linestyle='--')
    
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabel_style = {'size': 14, 'color': 'black', 'weight': 'semibold'}
    gl.ylabel_style = {'size': 14, 'color': 'black', 'weight': 'semibold'}
    # print(mask1[1])
    # plt.title(mask1[i], fontsize='14', fontweight='semibold')

    plt.savefig(PNG+"MADWRF_Nowcast_" + str(i).zfill(2) + ".png", dpi=300, 
    facecolor='w', edgecolor='w',orientation='lanscape', 
    papertype=None, format='png',bbox_inches='tight', 
    transparent=True, pad_inches=0.1)

# plt.show()