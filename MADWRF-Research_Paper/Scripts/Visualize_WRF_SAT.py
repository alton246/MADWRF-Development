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
import metpy
from datetime import datetime
from glob import glob
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset

PATH = '/home/alton/Github/MADWRF-Development/MADWRF-Research_Paper/Data/WRF_OUT/20210616/CLDTOPZ_CLDBASEZ/'
file = 'wrfout_ghi_d04_2021-06-16_12:00:00'

PATH_SAT = '/home/alton/Github/MADWRF-Development/MADWRF-Research_Paper/Data/Satellite_Imagery/True_Color/20210616/'
PNG = '/home/alton/Github/MADWRF-Development/MADWRF-Research_Paper/Figure/'

plt.rcParams['font.weight']='semibold'
plt.rcParams['font.size']='9'

v = np.linspace(0, 1100, 23, endpoint=True)

irr = Dataset(PATH+file, 'r')
mask1 = pd.date_range("2021-06-16 12:00:00", freq="15T", periods=25)

# for i in range(len(mask1)):
#         print(i, mask1[i])

states = NaturalEarthFeature(category="cultural", scale="10m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
x1,x2=irr['XLONG'][0].data.min(),irr['XLONG'][0].data.max()
y1,y2=irr['XLAT'][0].data.min(),irr['XLAT'][0].data.max()

fig = plt.figure(figsize=(22,8))
    
ax1 = plt.subplot(2,3,1,projection = ccrs.PlateCarree())

swh1 = plt.contourf(irr['XLONG'][14][0,:], irr['XLAT'][14][:,0], irr['SWDOWN2'][14], v,
                  transform=ccrs.PlateCarree(),cmap=cmo.solar)

ax1.add_feature(states, linewidth=0.5, edgecolor="black")
    
ax1.coastlines(resolution='10m')                
ax1.set_extent([x1,x2,y1,y2])
ax = plt.scatter(-59.4289, 13.1627, marker='*', color='m', s=75)
ax = plt.scatter(-59.6245, 13.1499, marker='o', color='c', s=75)
ax = plt.text(-59.90, 13.45, 'a)', color='k', style='normal',fontsize='12')
gl = ax1.gridlines(draw_labels=True,
            linewidth=0.5, color='black', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}
gl.ylabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}
plt.title(mask1[14], fontsize='9', fontweight='semibold')

# # plt.show()

ax2 = plt.subplot(2,3,2,projection = ccrs.PlateCarree())

swh2 = plt.contourf(irr['XLONG'][15][0,:], irr['XLAT'][15][:,0], irr['SWDOWN2'][15], v,
                  transform=ccrs.PlateCarree(),cmap=cmo.solar)
ax2.add_feature(states, linewidth=0.5, edgecolor="black")
ax2.coastlines(resolution='10m')                
ax2.set_extent([x1,x2,y1,y2])
ax = plt.scatter(-59.4289, 13.1627, marker='*', color='m', s=75)
ax = plt.scatter(-59.6245, 13.1499, marker='o', color='c', s=75)
ax = plt.text(-59.90, 13.45, 'b)', color='k', style='normal',fontsize='12')
gl = ax2.gridlines(draw_labels=True,
            linewidth=0.5, color='black', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}
gl.ylabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}
plt.title(mask1[15], fontsize='9', fontweight='semibold')


ax3 = plt.subplot(2,3,3,projection = ccrs.PlateCarree())

swh3 = plt.contourf(irr['XLONG'][16][0,:], irr['XLAT'][16][:,0], irr['SWDOWN2'][16], v,
                  transform=ccrs.PlateCarree(),cmap=cmo.solar)
ax3.add_feature(states, linewidth=0.5, edgecolor="black")
ax3.coastlines(resolution='10m')                
ax3.set_extent([x1,x2,y1,y2])
ax = plt.scatter(-59.4289, 13.1627, marker='*', color='m', s=75)
ax = plt.scatter(-59.6245, 13.1499, marker='o', color='c', s=75)
ax = plt.text(-59.90, 13.45, 'c)', color='k', style='normal',fontsize='12')

gl = ax3.gridlines(draw_labels=True,
            linewidth=0.5, color='black', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}
gl.ylabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}
plt.title(mask1[16], fontsize='9', fontweight='semibold')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.01, 0.7])
cb = fig.colorbar(swh3, cax=cbar_ax, ticks =v)
cb.ax.set_ylabel('Global Horizontal Irradiance ($W/m^{2}$)', fontsize=9,
                 fontweight='semibold')
# cb.ax.set_xlabel('Global Horizontal Irradiance ($W/m^{2}$', fontsize=9,
#                  fontweight='semibold')  # Set Label of colorbar


os.chdir(PATH_SAT)
extension = 'nc'
all_filenames = [i for i in sorted(glob('*.{}'.format(extension)))]

print(all_filenames)

F = xr.open_dataset(all_filenames[1])

        # Scan's start time, converted to datetime object
scan_start = datetime.strptime(F.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

# Scan's end time, converted to datetime object
scan_end = datetime.strptime(F.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')

# File creation time, convert to datetime object
file_created = datetime.strptime(F.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')

# The 't' variable is the scan's midpoint time
midpoint = str(F['t'].data)[:-8]
scan_mid = datetime.strptime(midpoint, '%Y-%m-%dT%H:%M:%S.%f')

# We'll use the `CMI_C02` variable as a 'hook' to get the CF metadata.
dat = F.metpy.parse_cf('CMI_C02')

geos = dat.metpy.cartopy_crs


    # Load the RGB arrays
R = F['CMI_C02'][:].data
G = F['CMI_C03'][:].data
B = F['CMI_C01'][:].data

# Apply range limits for each channel. RGB values must be between 0 and 1
R = np.clip(R, 0, 1)
G = np.clip(G, 0, 1)
B = np.clip(B, 0, 1)

# Apply the gamma correction
gamma = 2.2
R = np.power(R, 1/gamma)
G = np.power(G, 1/gamma)
B = np.power(B, 1/gamma)

# Calculate the "True" Green
G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G
G_true = np.clip(G_true, 0, 1)

# The final RGB array :)
RGB = np.dstack([R, G_true, B])

cleanIR = F['CMI_C13'].data

# Normalize the channel between a range.
#       cleanIR = (cleanIR-minimumValue)/(maximumValue-minimumValue)
cleanIR = (cleanIR-90)/(313-90)

# Apply range limits to make sure values are between 0 and 1
cleanIR = np.clip(cleanIR, 0, 1)

# Invert colors so that cold clouds are white
cleanIR = 1 - cleanIR

# Lessen the brightness of the coldest clouds so they don't appear so bright
# when we overlay it on the true color image.
cleanIR = cleanIR/1.4

# Yes, we still need 3 channels as RGB values. This will be a grey image.
RGB_cleanIR = np.dstack([cleanIR, cleanIR, cleanIR])

RGB_ColorIR = np.dstack([np.maximum(R, cleanIR), np.maximum(G_true, cleanIR),
                            np.maximum(B, cleanIR)])
x = dat.x
y = dat.y

ax4 = fig.add_subplot(2, 3, 4, projection=ccrs.PlateCarree())
        # ax.set_extent([-62.5,-47.5,  6.5, 16.5], crs=pc)
# ax4.set_extent([-60.5,-58.5, 12.5, 14.5], crs=ccrs.PlateCarree())
ax4.set_extent([x1,x2, y1, y2], crs=ccrs.PlateCarree())
ax4.imshow(RGB_ColorIR, origin='upper',
extent=(x.min(), x.max(), y.min(), y.max()),
        transform=geos, interpolation='none')
ax4.coastlines(resolution='10m', color='yellow', linewidth=1)
ax4.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=1)
ax4.add_feature(states, linewidth=0.5, edgecolor="yellow")
ax = plt.scatter(-59.4289, 13.1627, marker='*', color='m', s=75)
ax = plt.scatter(-59.6245, 13.1499, marker='o', color='c', s=75)
ax = plt.text(-59.90, 13.45, 'd)', color='k', style='normal',fontsize='12')
gl = ax4.gridlines(draw_labels=True,
            linewidth=0.5, color='black', alpha=0.5, linestyle='--')
    
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}
gl.ylabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}

plt.title('GOES-16 True Color', fontweight='bold', fontsize=9, loc='left')
plt.title('Full Disk\n{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
                      fontsize=8, loc='right')


F = xr.open_dataset(all_filenames[2])

        # Scan's start time, converted to datetime object
scan_start = datetime.strptime(F.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

# Scan's end time, converted to datetime object
scan_end = datetime.strptime(F.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')

# File creation time, convert to datetime object
file_created = datetime.strptime(F.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')

# The 't' variable is the scan's midpoint time
midpoint = str(F['t'].data)[:-8]
scan_mid = datetime.strptime(midpoint, '%Y-%m-%dT%H:%M:%S.%f')

# We'll use the `CMI_C02` variable as a 'hook' to get the CF metadata.
dat = F.metpy.parse_cf('CMI_C02')

geos = dat.metpy.cartopy_crs


    # Load the RGB arrays
R = F['CMI_C02'][:].data
G = F['CMI_C03'][:].data
B = F['CMI_C01'][:].data

# Apply range limits for each channel. RGB values must be between 0 and 1
R = np.clip(R, 0, 1)
G = np.clip(G, 0, 1)
B = np.clip(B, 0, 1)

# Apply the gamma correction
gamma = 2.2
R = np.power(R, 1/gamma)
G = np.power(G, 1/gamma)
B = np.power(B, 1/gamma)

# Calculate the "True" Green
G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G
G_true = np.clip(G_true, 0, 1)

# The final RGB array :)
RGB = np.dstack([R, G_true, B])

cleanIR = F['CMI_C13'].data

# Normalize the channel between a range.
#       cleanIR = (cleanIR-minimumValue)/(maximumValue-minimumValue)
cleanIR = (cleanIR-90)/(313-90)

# Apply range limits to make sure values are between 0 and 1
cleanIR = np.clip(cleanIR, 0, 1)

# Invert colors so that cold clouds are white
cleanIR = 1 - cleanIR

# Lessen the brightness of the coldest clouds so they don't appear so bright
# when we overlay it on the true color image.
cleanIR = cleanIR/1.4

# Yes, we still need 3 channels as RGB values. This will be a grey image.
RGB_cleanIR = np.dstack([cleanIR, cleanIR, cleanIR])

RGB_ColorIR = np.dstack([np.maximum(R, cleanIR), np.maximum(G_true, cleanIR),
                            np.maximum(B, cleanIR)])
x = dat.x
y = dat.y

ax5 = fig.add_subplot(2, 3, 5, projection=ccrs.PlateCarree())
        # ax.set_extent([-62.5,-47.5,  6.5, 16.5], crs=pc)
# ax4.set_extent([-60.5,-58.5, 12.5, 14.5], crs=ccrs.PlateCarree())
ax5.set_extent([x1,x2, y1, y2], crs=ccrs.PlateCarree())
ax5.imshow(RGB_ColorIR, origin='upper',
extent=(x.min(), x.max(), y.min(), y.max()),
        transform=geos, interpolation='none')
ax5.coastlines(resolution='10m', color='yellow', linewidth=1)
ax5.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=1)
ax5.add_feature(states, linewidth=0.5, edgecolor="yellow")
ax = plt.scatter(-59.4289, 13.1627, marker='*', color='m', s=75)
ax = plt.scatter(-59.6245, 13.1499, marker='o', color='c', s=75)
ax = plt.text(-59.90, 13.45, 'e)', color='k', style='normal',fontsize='12')
gl = ax5.gridlines(draw_labels=True,
            linewidth=0.5, color='black', alpha=0.5, linestyle='--')
    
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}
gl.ylabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}

plt.title('GOES-16 True Color', fontweight='bold', fontsize=9, loc='left')
plt.title('Full Disk\n{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
                      fontsize=8, loc='right')

F = xr.open_dataset(all_filenames[3])

        # Scan's start time, converted to datetime object
scan_start = datetime.strptime(F.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

# Scan's end time, converted to datetime object
scan_end = datetime.strptime(F.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')

# File creation time, convert to datetime object
file_created = datetime.strptime(F.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')

# The 't' variable is the scan's midpoint time
midpoint = str(F['t'].data)[:-8]
scan_mid = datetime.strptime(midpoint, '%Y-%m-%dT%H:%M:%S.%f')

# We'll use the `CMI_C02` variable as a 'hook' to get the CF metadata.
dat = F.metpy.parse_cf('CMI_C02')

geos = dat.metpy.cartopy_crs


    # Load the RGB arrays
R = F['CMI_C02'][:].data
G = F['CMI_C03'][:].data
B = F['CMI_C01'][:].data

# Apply range limits for each channel. RGB values must be between 0 and 1
R = np.clip(R, 0, 1)
G = np.clip(G, 0, 1)
B = np.clip(B, 0, 1)

# Apply the gamma correction
gamma = 2.2
R = np.power(R, 1/gamma)
G = np.power(G, 1/gamma)
B = np.power(B, 1/gamma)

# Calculate the "True" Green
G_true = 0.48358168 * R + 0.45706946 * B + 0.06038137 * G
G_true = np.clip(G_true, 0, 1)

# The final RGB array :)
RGB = np.dstack([R, G_true, B])

cleanIR = F['CMI_C13'].data

# Normalize the channel between a range.
#       cleanIR = (cleanIR-minimumValue)/(maximumValue-minimumValue)
cleanIR = (cleanIR-90)/(313-90)

# Apply range limits to make sure values are between 0 and 1
cleanIR = np.clip(cleanIR, 0, 1)

# Invert colors so that cold clouds are white
cleanIR = 1 - cleanIR

# Lessen the brightness of the coldest clouds so they don't appear so bright
# when we overlay it on the true color image.
cleanIR = cleanIR/1.4

# Yes, we still need 3 channels as RGB values. This will be a grey image.
RGB_cleanIR = np.dstack([cleanIR, cleanIR, cleanIR])

RGB_ColorIR = np.dstack([np.maximum(R, cleanIR), np.maximum(G_true, cleanIR),
                            np.maximum(B, cleanIR)])
x = dat.x
y = dat.y

ax6 = fig.add_subplot(2, 3, 6, projection=ccrs.PlateCarree())
        # ax.set_extent([-62.5,-47.5,  6.5, 16.5], crs=pc)
# ax4.set_extent([-60.5,-58.5, 12.5, 14.5], crs=ccrs.PlateCarree())
ax6.set_extent([x1,x2, y1, y2], crs=ccrs.PlateCarree())
ax6.imshow(RGB_ColorIR, origin='upper',
extent=(x.min(), x.max(), y.min(), y.max()),
        transform=geos, interpolation='none')
ax6.coastlines(resolution='10m', color='yellow', linewidth=1)
ax6.add_feature(ccrs.cartopy.feature.BORDERS, linewidth=1)
ax6.add_feature(states, linewidth=0.5, edgecolor="yellow")
ax = plt.scatter(-59.4289, 13.1627, marker='*', color='m', s=75)
ax = plt.scatter(-59.6245, 13.1499, marker='o', color='c', s=75)
ax = plt.text(-59.90, 13.45, 'f)', color='k', style='normal',fontsize='12')
gl = ax6.gridlines(draw_labels=True,
            linewidth=0.5, color='black', alpha=0.5, linestyle='--')
    
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}
gl.ylabel_style = {'size': 9, 'color': 'black', 'weight': 'semibold'}

plt.title('GOES-16 True Color', fontweight='bold', fontsize=9, loc='left')
plt.title('Full Disk\n{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
                      fontsize=8, loc='right')

plt.savefig(PNG + 'MADWRF_Satellite_Images.png', facecolor='w', edgecolor='w', dpi=300,
#            plt.savefig(PATH_PNG + 'GOES-16_True_Color.png', facecolor='w', edgecolor='w', dpi=300,
                       orientation='lanscape', papertype=None, format='png',
                       bbox_inches='tight', pad_inches=0.1)

plt.show()
