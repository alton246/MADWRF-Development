import os
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy.interpolate import griddata
from glob import glob



def QualityControlMETAR(PATH, filename, dataframe_name):
    

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
        # print (final_df.columns[final_df.isna().any()])

#        comb_data.append(final_df)
        return final_df


def InterpolateMETAR(df, start_lat, end_lat, delta_lat, start_lon, end_lon, delta_lon):
        """
        Interpolates the METAR data for each availalbe point within the domain
        """

        y_grid = np.arange(-7.0000, 16.8996, 0.01938)
        x_grid = np.arange(-100.0000, -53.5876, 0.02101)

        X,Y = np.meshgrid(x_grid, y_grid)


        points = np.random.rand(len(df), 2)
        # print(len(df_select))
        # grid_z0 = griddata(points, df['skyl1'], (X, Y), method = 'cubic')
        Z = griddata([(x,y) for x,y in zip(df['lon'],df['lat'])], 
                df['skyl1'], (X, Y), method='nearest')
        
        return X, Y, Z

def CombineMETARByTime(PATH):
        files = glob(PATH + '*.txt')
        comb_data = []
        for file in sorted(files):
                df = pd.read_csv(file, skiprows=5,sep=',')
                print('Working on station: ' + os.path.basename(file)[0:4])
                data = QualityControlMETAR(PATH, os.path.basename(file), df)
    
                comb_data.append(data)
        append_data = pd.concat(comb_data)

        append_data = append_data.sort_index()

        df_select = append_data.loc['2020-06-22 16:00:00']
        df_select['skyl1'] = df_select['skyl1'].astype(float)

        return df_select


def readCoastLine(file):
    '''
    Read NOAA coastline with -99.0 value for empty point
    Return two list of the two corrdinates X and Y
    '''
    f = open(file, 'r')
    cl = f.readlines(-1)
    f.close()
    CL = []
    for i in range(0, len(cl)):
        CL.append(np.array(cl[i].split(), dtype='f'))
    cl = np.reshape(CL, (len(cl), 2))
    cl = np.ma.masked_where(cl == -99, cl)
    return cl[:, 0], cl[:, 1]

#################################################################
###################### Program Begins Here ######################
#################################################################

PATH = '/home/alton/WRF_OUT/Sat_Reprojection/Data/METAR_CLDBASE_HGHT/'
PNG = '/home/alton/Github/MADWRF-Development/METAR_Scripts/Plots/'
coast_line_dir = '/home/alton/Github/MADWRF-Development/Data/Coastline_File/'
coast_line_file = 'merdescaraibe_m.dat'

DPI = 300
plt.rcParams['font.weight']='semibold'
plt.rcParams['font.size']='15'

#Quality Control and combine data
df_select = CombineMETARByTime(PATH)

#Interpolation

X, Y, Z = InterpolateMETAR(df_select, -7.0000, 16.8996, 0.01938, -100.0000, -53.5876, 0.02101)

#Visualizing the interpolation
extent = [-53.5876, -100.0000, -7.0000, 16.8996]
cl = readCoastLine(coast_line_dir + coast_line_file)
fig = plt.figure(figsize=(12, 10))
ax = plt.subplot(111)
mesh = ax.pcolormesh(X,Y,Z, cmap='rainbow')
ax.plot(cl[0], cl[1], color='black')
ax.set_xlim(extent[1], extent[0])
ax.set_ylim(extent[2], extent[3])
ax.set_aspect('equal')
ax.set_xlabel('Longitude', color='black', fontweight='demi', fontsize=12)
ax.set_ylabel('Latitude', color='black', fontweight='demi', fontsize=12)

#Visualize with cartopy options
# ax = plt.subplot(111,projection = ccrs.PlateCarree())

# ax.set_extent(extent)

# ax.add_feature(cartopy.feature.OCEAN)
# ax.add_feature(cartopy.feature.LAND, facecolor='tan')
# ax.add_feature(cartopy.feature.BORDERS, edgecolor='black')
# ax.coastlines(resolution='10m')
# ax.set_xlabel('Longitude', color='black', fontweight='demi', fontsize=12)
# ax.set_ylabel('Latitude', color='black', fontweight='demi', fontsize=12)
# mesh = ax.pcolormesh(X,Y,Z)

# cb = plt.colorbar(mesh, orientation="horizontal", fraction=0.07,anchor=(1.0,0.0))

plt.savefig(PNG+'2020-06-22_16_00_cldbase_cubic.png', dpi=DPI, facecolor='w', edgecolor='w',
           orientation='lanscape', papertype=None, format='png',
             bbox_inches='tight', pad_inches=0.1)

plt.show()

#cb.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))



