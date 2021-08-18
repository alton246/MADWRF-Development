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
    print(dlon, dlat)
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    print(a)
    c = 2 * asin(sqrt(a))
    km = 6367 * c

    return km

dist = calc_distance(12.9500, -59.4289, 13.2150, -59.4289)

print(dist)