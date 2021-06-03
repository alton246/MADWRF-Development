#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 20:51:06 2020

@author: alton

note : jfdorville
    No semicolon needs at the end of a line


"""

from netCDF4 import Dataset
import netCDF4 as nc4
import numpy as np
import os
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import time


def lat_lon_reproj(nc_folder, nc_indx):
    """
    Reprojects radians to lat lon
    last modification : 10/19/2020
    """
    os.chdir(nc_folder)
    full_direc = os.listdir()
    nc_files = [ii for ii in full_direc if ii.endswith('.nc')]
    g16_data_file = nc_files[nc_indx] # select .nc file
    print(nc_files[nc_indx]) # print file name

    # designate dataset
    g16nc = Dataset(g16_data_file, 'r')
    var_names = [ii for ii in g16nc.variables]
    var_name = var_names[0]
    try:
        band_id = g16nc.variables['band_id'][:]
        band_id = ' (Band: {},'.format(band_id[0])
        band_wavelength = g16nc.variables['band_wavelength']
        band_wavelength_units = band_wavelength.units
        band_wavelength_units = '{})'.format(band_wavelength_units)
        band_wavelength = ' {0:.2f} '.format(band_wavelength[:][0])
        print('Band ID: {}'.format(band_id))
        print('Band Wavelength: {} {}'.format(band_wavelength,band_wavelength_units))
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
    data_time_grab = ((g16nc.time_coverage_end).replace('T',' ')).replace('Z','')
    data_long_name = g16nc.variables[var_name].long_name

    # close file when finished
    g16nc.close()
    g16nc = None

    # create meshgrid filled with radian angles
    lat_rad, lon_rad = np.meshgrid(lat_rad_1d, lon_rad_1d)

    # lat/lon calc routine from satellite radian angle vectors

    lambda_0 = (lon_origin*np.pi)/180.0

    a_var = np.power(np.sin(lat_rad), 2.0) +\
            (np.power(np.cos(lat_rad), 2.0) *\
             (np.power(np.cos(lon_rad), 2.0) + (((r_eq*r_eq)/(r_pol*r_pol)) *\
              np.power(np.sin(lon_rad), 2.0))))

    b_var = -2.0*H*np.cos(lat_rad)*np.cos(lon_rad)
    c_var = (H**2.0)-(r_eq**2.0)

    r_s = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)

    s_x = r_s*np.cos(lat_rad)*np.cos(lon_rad)
    s_y = - r_s*np.sin(lat_rad)
    s_z = r_s*np.cos(lat_rad)*np.sin(lon_rad)

    lat = (180.0/np.pi) *\
          (np.arctan(((r_eq*r_eq)/(r_pol*r_pol)) *\
                     ((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
    lon = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)

    # print test coordinates
#     print('{} N, {} W'.format(lat[318,1849],abs(lon[318,1849])))

    return lon,lat,data,data_units,data_time_grab,data_long_name,band_id,band_wavelength,band_wavelength_units,var_name


def lat_lon_reproj2(nc_folder, nc_indx):
    """
    Reprojects radians to lat lon
    return a dictionnary
    created : October 2020
    last modification : 10/19/2020
    """
    os.chdir(nc_folder)
    full_direc = os.listdir()
    nc_files = [ii for ii in full_direc if ii.endswith('.nc')]
    g16_data_file = nc_files[nc_indx] # select .nc file
    # print('Working on:' + nc_files[nc_indx])  # print file name
    year = nc_files[nc_indx][23:27]
    dayjulian = int(nc_files[nc_indx][27:30]) - 1
    daynormal = datetime.datetime(int(year),1,1) + datetime.timedelta(int(dayjulian))
    daynormal2 = daynormal.strftime('%Y-%m-%d')
    HH = nc_files[nc_indx][30:32] 
    MM = nc_files[nc_indx][32:34]
    print('Working on:' + daynormal2 + '_' + HH + '_' + MM)

    # designate dataset
    g16nc = Dataset(g16_data_file, 'r')
    var_names = [ii for ii in g16nc.variables]
    var_name = var_names[0]
    try:
        band_id = g16nc.variables['band_id'][:]
        band_id = ' (Band: {},'.format(band_id[0])
        band_wavelength = g16nc.variables['band_wavelength']
        band_wavelength_units = band_wavelength.units
        band_wavelength_units = '{})'.format(band_wavelength_units)
        band_wavelength = ' {0:.2f} '.format(band_wavelength[:][0])
        print('Band ID: {}'.format(band_id))
        print('Band Wavelength: {} {}'.format(band_wavelength,
                                              band_wavelength_units))
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
    return data, daynormal2, HH, MM,


def get_lat(lat,min_lat,max_lat):
    """
    Retrieves the desired Latitudes
    """
    lat_new = []
    value = lat.shape[0]
    assign = lat_new.append
    for i in range(value):
        for j in range(len(lat[i])):
            if lat[i][j] > min_lat and lat[i][j] < max_lat:
                assign(lat[i][j])
    return lat_new;


def get_lon(lon, min_lon, max_lon):
    """
    Retrieves the desired Longitudes
    """
    lon_new = []
    value = lon.shape[0]
    assign = lon_new.append
    for i in range(value):
        for j in range(len(lon[i])):
#         print(lat[i][j])
#     print(len(lat[i]))
            if lon[i][j] > min_lon and lon[i][j] < max_lon:
#             print(lat[i][j])
                assign(lon[i][j])
    return lon_new


def proj_area(lat, lon, data, max_lat, min_lat, max_lon, min_lon):

    """
    Retrieves the latitude, longitude and data of
    the projected area of interest
    """
    real_data = []
    real_lon = []
    real_lat = []
    indx, indy = np.where((np.isnan(data['lat'])==False) *
                          (np.isnan(data['lon'])==False))
    lat_s = data['lat'][indx[np.argmin(indx)]:indx[np.argmax(indx)],
                       indy[np.argmin(indx)]]
    lon_s = data['lon'][indx[np.argmin(indy)],
                        indy[np.argmin(indy)]:indx[np.argmax(indy)]]
    LON, LAT = np.meshgrid(lon_s, lat_s)
    VAL = data['data'][indx[np.argmin(indx)]:indx[np.argmax(indx)],
                       indy[np.argmin(indy)]:indx[np.argmax(indy)]]
    indxs, indys = np.where((LON >= roi[3]) * (LON <= roi[2]) *
                            (LAT >= roi[1]) * (LAT <= roi[0]))
    real_lat = LAT[indxs[0]:indxs[-1], indys[0]:indys[-1]]
    real_lon = LON[indxs[0]:indxs[-1], indys[0]:indys[-1]]
    real_data = VAL[indxs[0]:indxs[-1], indys[0]:indys[-1]]
    return real_lat, real_lon, real_data


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


if __name__ == '__main__':

    DPI = 250

    #################################################################
    #                 Programm begins here                          #
    #################################################################

    coast_line_file = 'merdescaraibe_m.dat'

    # tic = timeit.default_timer()
    base_dir ='/home/alton/Github/MADWRF-Development'
    nc_folder = base_dir + '/Data/Clear_Sky_Mask/'
    png = base_dir + '/Satellite_Visualization_Scripts/Plots/'
    coast_line_dir = base_dir + '/Data/Coastline_File/'
    reproj_dir = base_dir + '/Data/Clear_Sky_Mask/Reprojected_Files/'
#    nc_folder = '/home/jfdorville/Recherche/AltonDaley/project_cloud/'

    file_indx = 0
    # print(daynormal2)
    # print(HH)
    # print(MM)
    data, daynormal2, HH, MM, = lat_lon_reproj2(nc_folder, file_indx)
    roi = [16.8996, -7.0000, -53.5876, -100.0000]
#    roi = [16.8996, 8.91853, -53.5876, -65.5484]
#    roi = [17, 9, -52, -67]
#    roi = [13.6, 12.7, -58.95, -60.5]
    data = proj_area2(data, roi)

    # visu
    DPI = 150
#    DPI = 600
    plt.rcParams['font.weight']='semibold'
    plt.rcParams['font.size']='15'
    cl = readCoastLine(coast_line_dir + coast_line_file)
    fig = plt.figure(figsize=(12, 10))
    ax1 = fig.add_subplot(111)
    mesh = ax1.pcolormesh(data['real_lon'], data['real_lat'],
                          data['real_data'])
    ax1.plot(cl[0], cl[1], color='white')
    ax1.set_xlim([data['real_lon'].min(), data['real_lon'].max()])
    ax1.set_ylim([data['real_lat'].min(), data['real_lat'].max()])
    ax1.set_aspect('equal')
    ax1.set_xlabel('Longitude', color='black', fontweight='demi', fontsize=12)
    ax1.set_ylabel('Latitude', color='black', fontweight='demi', fontsize=12)
    daynormal2 + '_' + HH + '_' + MM
    plt.savefig(png + daynormal2 + '_' + HH + '_' + MM + '_cldmask_2.png',
                format='png', transparent=False, dpi=DPI,
                bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)

    real_data = np.ma.masked_where(data['real_data']==0 , data['real_data'])
    fig = plt.figure(figsize=(12, 10))
    ax1 = fig.add_subplot(111)
    mesh = ax1.pcolormesh(data['real_lon'], data['real_lat'],
                          real_data, cmap=plt.cm.hsv)
    ax1.plot(cl[0], cl[1], color='k')
    ax1.set_xlim([data['real_lon'].min(), data['real_lon'].max()])
    ax1.set_ylim([data['real_lat'].min(), data['real_lat'].max()])
    ax1.set_aspect('equal')
    ax1.set_xlabel('Longitude', color='black', fontweight='demi', fontsize=12)
    ax1.set_ylabel('Latitude', color='black', fontweight='demi', fontsize=12)
    plt.savefig(png + daynormal2 + '_' + HH + '_' + MM + '_cldmask_3.png',
                format='png', transparent=False, dpi=DPI,
                bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)


    fn = daynormal2 + '_' + HH + '_' + MM + '.nc'
    ds = nc4.Dataset(reproj_dir + fn, 'w', format='NETCDF4')
    
    #Create root groups
    bcm = ds.createGroup("GOES BCM")

    # Create Dimensions
#    time = ds.createDimension('time', 1)
    lats = ds.createDimension('lat', data['real_lat'].shape[0])
    lons = ds.createDimension('lon', data['real_lon'].shape[1])

    # Add NetCDF Variables
#    times = ds.createVariable('time', 'f4', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    BCM = ds.createVariable('value', 'f4', ('lat', 'lon',),
                              zlib=True, fill_value=-1)
    

    lats[:] = data['real_lat'][:, 0]
    lons[:] = data['real_lon'][0, :]
    BCM[:] = data['real_data']
    
    #Create Attributes
    
    ds.description = "Projection of Binary Cloud Mask Around the eastern Caribbean"
#    ds.history = "Created " + time.ctime(time.time())
    ds.source = "GOES Imagery"
    lats.units = "degrees north"
    lons.units = "degrees east"
    BCM.units = "Unknown"
    ds.close()
    #toc = timeit.default_timer()
    #print(toc-tic)

    #fig = plt.figure(figsize=(10,15),dpi=900)
    ###
    #n_add = 0
    #bbox = [np.min(lon1),np.min(lat1),np.max(lon1),np.max(lat1)]
    #m = Basemap(llcrnrlon=bbox[0]-n_add,llcrnrlat=bbox[1]-n_add,urcrnrlon=bbox[2]+n_add,urcrnrlat=bbox[3]+n_add,projection='merc', resolution='i')
    #m.drawcoastlines(linewidth=0.5)
    #m.drawcountries(linewidth=0.25)
    #m.pcolormesh(lon1, lat1, DATA, latlon=True)
    #parallels = np.linspace(np.min(lat),np.max(lat),5.)
    #m.drawparallels(parallels,labels=[True,False,False,False])
    #meridians = np.linspace(np.min(lon),np.max(lon),5.)
    #m.drawmeridians(meridians,labels=[False,False,False,True])
    #cb = m.colorbar()
    ##
    #plt.show()

#    exit()
