import netCDF4 as nc4
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import statistics
import datetime
from netCDF4 import Dataset
from math import cos, sin, asin, sqrt, radians
from mpl_toolkits.basemap import Basemap, cm
import textwrap
from glob import glob




def proj_area2(data, roi):
    """
    Retrieves the latitude, longitude and data of
    the projected area of interest as array
    lat and lon don't have NaN
    from dictionary
        lon = data['lon']; lat = data['lat']
        value = data['data']
    roi = [max_lat, min_lat, max_lon, min_lon]
    """
    lat_center = int(data['lat'].shape[0]/2)
    lon_center = int(data['lon'].shape[1]/2)
    ind_lat = np.where((np.isnan(data['lat'][:, lon_center])==False))[0]
    ind_lon = np.where((np.isnan(data['lon'][lat_center, :])==False))[0]
    lon_s = data['lon'][lat_center, ind_lon[0]:ind_lon[-1]]
    lat_s = data['lat'][ind_lat[0]:ind_lat[-1], lon_center]
    LON, LAT = np.meshgrid(lon_s, lat_s)
    VAL = data['data'][ind_lat[0]:ind_lat[-1],
                       ind_lon[0]:ind_lon[-1]]
    indxs, indys = np.where((LON >= roi[3]) * (LON <= roi[2]) *
                            (LAT >= roi[1]) * (LAT <= roi[0]))
    data['real_lat'] = LAT[indxs[0]:indxs[-1], indys[0]:indys[-1]]
    data['real_lon'] = LON[indxs[0]:indxs[-1], indys[0]:indys[-1]]
    data['real_data'] = VAL[indxs[0]:indxs[-1], indys[0]:indys[-1]]
    return data

def calc_distance(lat1, lon1, lat2, lon2):
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

def write_SLAB(f,data):
    f_val = data.shape[1] * data.shape[0]
    count = 1
    for col in data:
        for val in col:
#            np.float(val)
#            print(type(val))
            if type(val) != np.float32:
                val = -999.9
                
            f.write('%s,' % val)
            if count % 5 == 0:
                f.write('%s, & \n' % val)
            if count == f_val:
                f.write('%s' % val)                
#            print(count,val)
            count = count + 1

def file_writer(data, daynormal2, HH, MM):
    f = open(out_dir + 'BRTEMP_' + daynormal2 + '_' + HH + ':' + MM + ':' + '00.f95', 'w')
    
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
    f.write('  real :: XFCST=' + HH + '\n')
    f.write("  character(len=8) :: STARTLOC='SWCORNER' \n")
    f.write("  character(len=9) :: FIELD= 'CLDTOPZ'\n" )
    f.write("  character(len=25) :: UNITS='km' \n")
    f.write("  character(len=46) :: DESC='Cloud Top Height' \n")
    f.write("  character(len=32) :: MAP_SOURCE='GOES Imagery' \n")
    f.write('  real :: XLVL=100 \n')
    f.write('  integer :: NX=' + str(data['real_lon'].shape[1]) +' \n')
    f.write('  integer :: NY=' + str(data['real_lat'].shape[0]) +' \n')
    f.write('  integer :: IPROJ=' + str(0) +' \n')
#    f.write('  real :: STARTLAT=' + str(round(min(data['real_lat'][:,0]), 5)) +' \n')
    f.write('  real :: STARTLAT=' + str(round(data['real_lat'][:, 0].min(), 5)) +' \n')
    f.write('  real :: STARTLON=' + str(round(min(data['real_lon'][0,:]), 5)) +' \n')
    f.write('  real :: DELTALAT=' + ' ' + str(round(data['real_lat'][:, 0][0] - data['real_lat'][:, 0][1], 5)) +' \n')
    f.write('  real :: DELTALON=' + ' ' + str(round(data['real_lon'][0,:][1] - data['real_lon'][0,:][0], 5)) +' \n')
    f.write('  real :: DX='+ str(round(calc_distance(0,0,0,0.018), 5)) +' \n')
    f.write('  real :: DY=' + str(round(calc_distance(0,0,0.019,0), 5)) +' \n')
    f.write('  real :: XLONC=' + str(round(statistics.median(data['real_lon'][0,:]), 5)) +' \n')
    f.write('  real :: TRUELAT1= 12.941 \n')
    f.write('  real :: TRUELAT2= 0 \n')
    f.write('  real :: NLATS='+ str(data['real_lat'].shape[0]) + '\n')
    f.write('  real :: EARTH_RADIUS = 6367470. * .001 \n')
    f.write('  logical :: IS_WIND_EARTH_REL = .FALSE. \n')
    f.write('! SLAB is an allocatable array, because we do not necessarily know in \n')
    f.write('! advance the size of the array we need to read. \n')
    f.write('  real, allocatable, dimension(:,:) :: SLAB \n')
    f.write('\n')
    f.write('  SLAB = reshape( (/ ')
    write_SLAB(f, data['real_data'][:])
##       f.write('  real, allocatable, dimension(' + str(data['real_lon'].shape[1]) + ',' +  str(data['real_lat'].shape[0]) + ')' +' :: SLAB \n')
    f.write(' /), (/ '+ str(data['real_lon'].shape[1]) +' , ' + str(data['real_lat'].shape[0]) + ' /) ) \n')
    f.write('\n')
    
#    f.write("     OPEN(UNIT=IUNIT, file=" + " 'GR_CTH:" +  daynormal2 + '_' + HH + ':' + MM + ':' + '00' + "', status='new', & \n")
    f.write("     OPEN(UNIT=IUNIT, file=" + "'BRTEMP:" +  daynormal2 + '_' + HH  + "', status='new', & \n")
#    f.write("     OPEN(UNIT=IUNIT, file=" + " 'CLDTOPZ:" +  daynormal2 + '_' + HH + ':' + MM + "', status='new', & \n")
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
    f.write('                      STARTLAT,    & \n')
    f.write('                      STARTLON,    & \n')
    f.write('                      DX,          & \n')
    f.write('                      DY,          & \n')
    f.write('                      XLONC,       & \n')
    f.write('                      TRUELAT1,    & \n')
    f.write('                      EARTH_RADIUS \n')
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


def gen_fortran_files(nc_folder, outdir, data):
    
    """
    Accepts a folder with GOES 16 CLOUDMASK files in netcdf.Reprojects radians 
    to lat lon return a dictionary from which a smaller area is extracted in 
    order to create the fortran files necessary to generate the Intermediate 
    WPS FILES
    created : January 2021
    last modification : 01/13/2020
    """
    
    for file in sorted(data):
#        print(os.path.basename(file))
        year = int(os.path.basename(file)[27:31])
        dayjulian = int(os.path.basename(file)[31:34]) - 1
        daynormal = datetime.datetime(year,1,1) + datetime.timedelta(dayjulian)
        daynormal2 = daynormal.strftime('%Y-%m-%d')
        HH = os.path.basename(file)[34:36] 
        MM = os.path.basename(file)[36:38]
        # print(year)
        # print(dayjulian)
        # print(daynormal2)
        # print(HH)
        # print(MM)
        if MM == '00':
            print('Working on ' + daynormal2 + '_' + HH + ':' + MM + ':' + '00......')
#        print('CLDMASK:' + daynormal2 + '_' + HH + ':' + MM + ':' + '00')
    
            ## designate dataset
            g16nc = Dataset(file, 'r')
            var_names = [ii for ii in g16nc.variables]
            var_name = var_names[0]
            try:
                band_id = g16nc.variables['band_id'][:]
                band_id = ' (Band: {},'.format(band_id[0])
                band_wavelength = g16nc.variables['band_wavelength']
                band_wavelength_units = band_wavelength.units
                band_wavelength_units = '{})'.format(band_wavelength_units)
                band_wavelength = ' {0:.2f} '.format(band_wavelength[:][0])
                # print('Band ID: {}'.format(band_id))
                # print('Band Wavelength: {} {}'.format(band_wavelength,
                                                #   band_wavelength_units))
            except:
                band_id = ''
                band_wavelength = ''
                band_wavelength_units = ''
        
            # GOES-R projection info and retrieving relevant constants
            proj_info = g16nc.variables['goes_imager_projection']
            lon_origin = proj_info.longitude_of_projection_origin
            H = proj_info.perspective_point_height+proj_info.semi_major_axis
            r_eq = proj_info.semi_major_axis
            r_pol = proj_info.semi_minor_axis
        
            # grid info
            lat_rad_1d = g16nc.variables['x'][:]
            lon_rad_1d = g16nc.variables['y'][:]
        
            # data info
            data = g16nc.variables[var_name][:]
            data_units = g16nc.variables[var_name].units
            data_time_grab = ((g16nc.time_coverage_end).replace('T',
                                                            ' ')).replace('Z', '')
            data_long_name = g16nc.variables[var_name].long_name
        
            # close file when finished
            g16nc.close()
            g16nc = None
        
            # create meshgrid filled with radian angles
            lat_rad, lon_rad = np.meshgrid(lat_rad_1d, lon_rad_1d)
        
            # lat/lon calc routine from satellite radian angle vectors
        
            lambda_0 = (lon_origin*np.pi)/180.0
        
            a_var = np.power(np.sin(lat_rad), 2.0) +\
                    (np.power(np.cos(lat_rad), 2.0) *
                    (np.power(np.cos(lon_rad), 2.0) + (((r_eq*r_eq)/(r_pol*r_pol)) *
                    np.power(np.sin(lon_rad), 2.0))))
        
            b_var = -2.0*H*np.cos(lat_rad)*np.cos(lon_rad)
            c_var = (H**2.0)-(r_eq**2.0)
        
            r_s = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)
        
            s_x = r_s*np.cos(lat_rad)*np.cos(lon_rad)
            s_y = - r_s*np.sin(lat_rad)
            s_z = r_s*np.cos(lat_rad)*np.sin(lon_rad)
        
            lat = (180.0/np.pi) *\
              (np.arctan(((r_eq*r_eq)/(r_pol*r_pol)) *
                         ((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
            lon = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)
        
            # print test coordinates
            #     print('{} N, {} W'.format(lat[318,1849],abs(lon[318,1849])))
        
            data = {'lon': lon, 'lat': lat, 'data': data, 'data_units': data_units,
                'data_time_grab': data_time_grab, 'data_long_name': data_long_name,
                'band_id': band_id, 'band_wavelength': band_wavelength,
                'band_wavelength_units': band_wavelength_units,
                'var_name': var_name}
#           print(data)
            #Specify sub area Here
            roi = [16.8996, -7.0000, -53.5876, -100.0000]
#           print('file')
            data_new = proj_area2(data, roi)
        
            fort_file = file_writer(data_new, daynormal2, HH, MM)





if __name__ == '__main__':

    for i in range(16,20):
        if len(str(i)) < 2:
            hour = str(i).zfill(2)
        else:
            hour = str(i)

        nc_folder = '/home/alton/WRF_OUT/Sat_Reprojection/Data/GOES/Channel13/20200622/' + hour + '/'
        out_dir = '/home/alton/Github/MADWRF-Development/Satellite_Reprojection_scripts/WPS_Intermediate_Files/BRTEMP/'
        files = glob(nc_folder+'OR_ABI-L2-CMIPF*.nc')


        
        gen_fortran_files(nc_folder, out_dir, files)