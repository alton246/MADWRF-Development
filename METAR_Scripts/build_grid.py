import os
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy.interpolate import griddata
from glob import glob



def remove_unwanted_dates(PATH, filename, dataframe_name):
    

    """ 
    Removes unwanted time periods because we are only interested in hourly
    data.
    """
    
    dataframe_name = pd.read_csv(PATH + filename, skiprows=5,
                     sep=',')
    
    if len(dataframe_name['valid']) != 0:
             
        dataframe_name['DateTime'] = pd.to_datetime(dataframe_name['valid'])
    
        dataframe_name = dataframe_name.set_index('DateTime')
    
        lon = dataframe_name['lon'][0]
        lat = dataframe_name['lat'][0]
    
#    dataframe_name = dataframe_name.drop(['valid'], axis=1)
    
        dates_rem = []
        for i in range(len(dataframe_name.index)):  
            if str(dataframe_name.index[i])[14:16] != '00':
                dates_rem.append(dataframe_name.index[i])
    
        """
        NB Index can not be removed individually it is best to create a list then 
        remove the complete list from the dataframe  
        """ 
        
        updated_df = dataframe_name.drop(dates_rem)
        """        
        Test for a 24 hours of data if the file has less than 24 hours of data
        the index is reset, otherwise we continue
        """     
        if len(updated_df.index) != 24:
#        print('yes')
            idx = pd.date_range(start = '2020-06-22 00:00:00', freq="H", periods=24)
#        print('no')
            updated_df = dataframe_name.reindex(idx)
            
        """
        Removes old the column titled value as it was the old DateTime values
        """
        
        final_df = updated_df.drop(['valid'], axis=1)
        
#    elif len(len(updated_df.index) == 0):
#        print(True)
        """
        Changes nans and 'M' to -999.9 and fills out the missing lat, lon and 
        station id
        """
        final_df['station'] = final_df['station'].replace([np.NaN], 
                  filename[0:4])
        final_df['lon'] = final_df['lon'].replace([np.nan], lon)
        final_df['lat'] = final_df['lat'].replace([np.nan], lat)
        final_df['skyl1'] = final_df['skyl1'].replace(['M', np.NaN], 
                  np.nan)
        final_df['skyl2'] = final_df['skyl2'].replace(['M', np.NaN], 
                  np.nan)
        final_df['skyl3'] = final_df['skyl3'].replace(['M', np.NaN], 
                  np.nan)

#        comb_data.append(final_df)
        
        
    
        return final_df

#################################################################
###################### Program Begins Here ######################
#################################################################

PATH = '/home/alton/WRF_OUT/Sat_Reprojection/Data/METAR_CLDBASE_HGHT/'

##Quality Ceck Data and combine them by hour of the day
files = glob(PATH + '*.txt')
comb_data = []
for file in sorted(files):
    df = pd.read_csv(file, skiprows=5,sep=',')
    print('Working on station: ' + os.path.basename(file)[0:4])
    data = remove_unwanted_dates(PATH, os.path.basename(file), df)
    
    comb_data.append(data)
append_data = pd.concat(comb_data)

append_data = append_data.sort_index()

df_select = append_data.loc['2020-06-22 10:00:00']

# print(df_select)

#Building grid to interpolate data onto


y_grid = np.arange(-7.0000, 16.8996, 0.01938)
x_grid = np.arange(-100.0000, -53.5876, 0.02101)

X,Y = np.meshgrid(x_grid, y_grid)


points = np.random.rand(len(df_select), 2)
# print(len(df_select))
grid_z0 = griddata(points, df_select['skyl1'], (X, Y), method = 'nearest')
Z = griddata([(x,y) for x,y in zip(df_select['lon'],df_select['lat'])], 
              df_select['skyl1'], (X, Y), method='nearest')
print(Z)
# plt.subplot(111)

extent = [-53.5876, -100.0000, 7.0000, 16.8996]
#
fig = plt.figure(figsize=(12, 8))
ax = plt.subplot(111)
# ax.set_extent(extent)

# ax.add_feature(cartopy.feature.OCEAN)
# ax.add_feature(cartopy.feature.LAND, facecolor='tan')
# ax.add_feature(cartopy.feature.BORDERS, edgecolor='black')
# ax.coastlines(resolution='10m')
mesh = ax.pcolormesh(X,Y,Z)
# plt.scatter(X, Y, 
#            grid_z0, 'r', edgecolor='w')
# plt.imshow(Z, extent=(-100.0000,-53.5876,
#                      7.0000,16.8996), 
                    #    origin='lower')
# cb = plt.colorbar(orientation="horizontal", fraction=0.07,anchor=(1.0,0.0))
plt.show()
#cb.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#plt.savefig(PNG+'test_base.png', dpi=150, facecolor='w', edgecolor='w',
#            orientation='lanscape', papertype=None, format='png',
#              bbox_inches='tight', pad_inches=0.1)

# plt.title('Nearest')
# plt.show()
# print(grid_z0)
# print(len(lats), len(lons))