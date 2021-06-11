import os
import statistics
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import cartopy
import statistics
import matplotlib.pyplot as plt
from math import cos, sin, asin, sqrt, radians
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
                  -999.9)
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
        grid_z0 = griddata(points, df['skyl1'], (X, Y), method = 'nearest')
        Z = griddata([(x,y) for x,y in zip(df['lon'],df['lat'])], 
                df['skyl1'], (X, Y), method='nearest')
        # print(X[0])
        # print(Y[0])
        # print(Y.shape)
        # print(Y[0][0])
        # print(Z)
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

def CalculateDistance(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    km = 6367 * c
    
    return km

def WriteSLAB(f, data):
        f_val = data.shape[0] * data.shape[1]
        count = 1
        for col in Z:                                         #loop through the array and collect each column
                for val in col:                               #loop through each column and get eac value in that column
                        f.write('%s,' % val)                  # Write value to file
                        if count % 10 == 0:                    #test if count is divisible by 5
                                f.write('%s, & \n' % val)     # if it s then terminate that line that it is written into and start a new
                        if count == f_val:                    # test if we are at the end of the data
                                f.write('%s' % val)           
                        count = count + 1


def FortranFileWriter(X, Y, Z, daynormal2, HH, MM):
        f = open(out_dir + 'CLDBASEZ' + daynormal2 + '_' + HH + ':' + MM + ':' + '00.f95', 'w')
    
        f.write('program sample_read \n')
        f.write('\n' )
        f.write('! Fortran 90 version. \n')
        f.write('!   This is a simple program to write data in the WPS intermediate \n')
        f.write('!   format.  It is included mainly as an aid to understanding said \n')
        f.write('!   format. \n')
        f.write('\n')
        f.write('  implicit none \n')
        f.write('\n')
        f.write('! Declarations: \n')
        f.write('\n')
        f.write('  integer, parameter :: IUNIT = 10 \n')
        f.write('  integer, parameter :: OUNIT = 11 \n')
        f.write('  integer :: ierr=0 \n')
        f.write('\n')
        f.write('  integer :: IFV=5 \n')
        f.write('\n')
        f.write("  character(len=24) :: HDATE = '" + daynormal2 + '_' + HH + ':' + MM + ':' + '00' +  "' \n")
        f.write('  real :: XFCST=0.000000 \n')
        f.write("  character(len=8) :: STARTLOC='SWCORNER' \n")
        f.write("  character(len=9) :: FIELD= 'CLDBASEZ'\n" )
        f.write("  character(len=25) :: UNITS='kilometers' \n")
        f.write("  character(len=46) :: DESC='Cloud Base Height' \n")
        f.write("  character(len=32) :: MAP_SOURCE='IOWA State METAR Network' \n")
        f.write('  real :: XLVL=100 \n')
        f.write('  integer :: NX=' + str(X.shape[1]) +' \n')
        f.write('  integer :: NY=' + str(Y.shape[0]) +' \n')
        f.write('  integer :: IPROJ=' + str(0) +' \n')
        f.write('  real :: STARTLAT=' + str(Y[0][0]) +' \n')
        f.write('  real :: STARTLON=' + str(X[0][0]) +' \n')
        f.write('  real :: DELTALAT=' + ' ' + str(0.01938) +' \n')
        f.write('  real :: DELTALON=' + ' ' + str(0.02101) +' \n')
        f.write('  real :: DX=' + str(round(CalculateDistance(0,0,0,0.018), 5)) +' \n')
        f.write('  real :: DY=' + str(round(CalculateDistance(0,0,0.019,0), 5)) +' \n')
        f.write('  real :: XLONC=' + str(round(statistics.median(X[0]), 5)) +' \n')
        f.write('  real :: TRUELAT1= 12.941 \n')
        f.write('  real :: TRUELAT2= 0 \n')
        f.write('  real :: NLATS=' + str(Y.shape[0]) +' \n')
        f.write('  real :: EARTH_RADIUS = 6367470. * .001 \n')
        f.write('  logical :: IS_WIND_EARTH_REL = .FALSE. \n')
        f.write('! SLAB is an allocatable array, because we do not necessarily know in \n')
        f.write('! advance the size of the array we need to read. \n')
        f.write('  real, allocatable, dimension(:,:) :: SLAB \n')
        f.write('\n')
        f.write('  SLAB = reshape( (/ ')
        WriteSLAB(f, Z)
        f.write(' /), (/ '+ str(X.shape[1]) +' , ' + str(Y.shape[0]) + ' /) ) \n')
        f.write('\n')
        f.write("     OPEN(UNIT=IUNIT, file=" + "'CLDBASEZ:" +  daynormal2 + '_' + HH + "', status='new', & \n")
        f.write("                  FORM='UNFORMATTED', convert='big_endian') \n")
        f.write('\n')
        f.write('     WRITE (IUNIT) IFV')
        f.write('\n')
        f.write('     ! WRITE the second record, common to all projections: \n')
        f.write('\n')
        f.write('     WRITE (IUNIT) HDATE,      & \n')
        f.write('                   XFCST,      & \n')
        f.write('                   MAP_SOURCE, & \n')
        f.write('                   FIELD,      & \n')
        f.write('                   UNITS,      & \n')
        f.write('                   DESC,       & \n')
        f.write('                   XLVL,       & \n')
        f.write('                   NX,         & \n')
        f.write('                   NY,         & \n')
        f.write('                   IPROJ \n')
        f.write('\n')
        f.write('     ! WRITE the third record, which depends on the projection: \n')
        f.write('\n')
        f.write('     if (IPROJ == 0) then  \n')
        f.write('\n')
        f.write('       !  This is the Cylindrical Equidistant (lat/lon) projection: \n')
        f.write('       WRITE (IUNIT) STARTLOC,    & \n')
        f.write('                     STARTLAT,    & \n')
        f.write('                     STARTLON,    & \n')
        f.write('                     DELTALAT,    & \n')
        f.write('                     DELTALON,    & \n')
        f.write('                     EARTH_RADIUS  \n')
        f.write('\n')
        f.write('     elseif (IPROJ == 1) then  \n')
        f.write('       ! This is the Mercator projection: \n')
        f.write('       WRITE (IUNIT) STARTLOC,  & \n')
        f.write('                     STARTLAT,  & \n')
        f.write('                     STARTLON,  & \n')
        f.write('                     DX,        & \n')
        f.write('                     DY,        & \n')
        f.write('                     TRUELAT1,  & \n')
        f.write('                     EARTH_RADIUS \n')
        f.write('\n')
        f.write('     elseif (IPROJ == 3) then \n')
        f.write('       ! This is the Lambert Conformal projection: \n')
        f.write('       WRITE (IUNIT) STARTLOC,     & \n')
        f.write('                     STARTLAT,     & \n')
        f.write('                     STARTLON,     & \n')
        f.write('                     DX,           & \n')
        f.write('                     DY,           & \n')
        f.write('                     XLONC,        & \n')
        f.write('                     TRUELAT1,     & \n')
        f.write('                     TRUELAT2,     & \n')
        f.write('                     EARTH_RADIUS   \n')
        f.write('\n')
        f.write('     elseif (IPROJ == 4) then \n')
        f.write('       ! This is the Gaussian projection \n')
        f.write('       WRITE (IUNIT) STARTLOC,    & \n')
        f.write('                     STARTLAT,    & \n')
        f.write('                     STARTLON,    & \n')
        f.write('                     NLATS,       & \n')
        f.write('                     DELTALON,    & \n')
        f.write('                     EARTH_RADIUS \n')
        f.write('\n')
        f.write('     elseif (IPROJ == 5) then \n')
        f.write('        ! This is the Polar Stereographic projection: \n')
        f.write('       WRITE (IUNIT) STARTLOC,    & \n')
        f.write('                     STARTLAT,    & \n')
        f.write('                     STARTLON,    & \n')
        f.write('                     DX,          & \n')
        f.write('                     DY,          & \n')
        f.write('                     XLONC,       & \n')
        f.write('                     TRUELAT1,    & \n')
        f.write('                     EARTH_RADIUS \n')
        f.write('\n')
        f.write('     endif \n')
        f.write('\n')
        f.write('\n')
        f.write('     WRITE (IUNIT) IS_WIND_EARTH_REL \n')
        f.write('     WRITE (IUNIT) SLAB \n')
        f.write('\n')
        f.write('     CLOSE(UNIT=IUNIT, IOSTAT=IERR) \n')
        f.write('\n')
        #    f.write("  write(*,'(/,"End of read loop.  Program finished.")') \n")
        f.write('end program sample_read \n')
        f.close()

#################################################################
###################### Program Begins Here ######################
#################################################################

PATH = '/home/alton/WRF_OUT/Sat_Reprojection/Data/METAR_CLDBASE_HGHT/'
PNG = '/home/alton/Github/MADWRF-Development/METAR_Scripts/Plots/'
out_dir = '/home/alton/Github/MADWRF-Development/Satellite_Reprojection_scripts/WPS_Intermediate_Files/CLDBZIN/'
coast_line_dir = '/home/alton/Github/MADWRF-Development/Data/Coastline_File/'
coast_line_file = 'merdescaraibe_m.dat'

files = glob(PATH + '*.txt')
comb_data = []
for file in sorted(files):
        df = pd.read_csv(file, skiprows=5,sep=',')
        print('Working on station: ' + os.path.basename(file)[0:4])
        data = QualityControlMETAR(PATH, os.path.basename(file), df)
    
        comb_data.append(data)
append_data = pd.concat(comb_data)

append_data = append_data.sort_index()
append_data['skyl1'] = append_data['skyl1'].astype(float)

for i in append_data.index.unique():
        daynormal = str(i)[0:10]
        HH = str(i)[11:13]
        MM = str(i)[14:16]
        SS = str(i)[17:19]
        df_select = append_data.loc[i]
        print("Interpolating " + str(i) + "....." )
        # print(str(i)[17:19])
        X, Y, Z = InterpolateMETAR(df_select, -7.0000, 16.8996, 0.01938, -100.0000, -53.5876, 0.02101)

        FortranFileWriter(X, Y, Z, daynormal, HH, MM)
# print(append_data.index.unique())
# print(append_data)

#Quality Control and combine data
# df_select = CombineMETARByTime(PATH)

#Interpolation

# X, Y, Z = InterpolateMETAR(df_select, -7.0000, 16.8996, 0.01938, -100.0000, -53.5876, 0.02101)

#Visualizing the interpolation
# extent = [-53.5876, -100.0000, -7.0000, 16.8996]
# cl = readCoastLine(coast_line_dir + coast_line_file)
# fig = plt.figure(figsize=(12, 10))
# ax = plt.subplot(111)
# mesh = ax.pcolormesh(X,Y,Z, cmap='rainbow')
# ax.plot(cl[0], cl[1], color='black')
# ax.set_xlim(extent[1], extent[0])
# ax.set_ylim(extent[2], extent[3])
# ax.set_aspect('equal')
# ax.set_xlabel('Longitude', color='black', fontweight='demi', fontsize=12)
# ax.set_ylabel('Latitude', color='black', fontweight='demi', fontsize=12)

#Visualize with cartopy options
# ax = plt.subplot(111,projection = ccrs.PlateCarree())

# ax.set_extent(extent)





