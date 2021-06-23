#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 22:15:10 2021

@author: alton
"""

from __future__ import print_function
import json
import time
import datetime

# Python 2 and 3: alternative 4
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

# Number of attempts to download data
MAX_ATTEMPTS = 6
# HTTPS here can be problematic for installs that don't have Lets Encrypt CA
SERVICE = "http://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"


def download_data(uri):
    """Fetch the data from the IEM
    The IEM download service has some protections in place to keep the number
    of inbound requests in check.  This function implements an exponential
    backoff to keep individual downloads from erroring.
    Args:
      uri (string): URL to fetch
    Returns:
      string data
    """
    attempt = 0
    while attempt < MAX_ATTEMPTS:
        try:
            data = urlopen(uri, timeout=300).read().decode("utf-8")
            if data is not None and not data.startswith("ERROR"):
                return data
        except Exception as exp:
            print("download_data(%s) failed with %s" % (uri, exp))
            time.sleep(5)
        attempt += 1

    print("Exhausted attempts to download, returning empty data")
    return ""


def get_stations_from_filelist(filename):
    """Build a listing of stations from a simple file listing the stations.
    The file should simply have one station per line.
    """
    stations = []
    for line in open(filename):
        stations.append(line.strip())
    return stations


def get_stations_from_networks():
    """Build a station list by using a bunch of IEM networks."""
    stations = []
#    countries = """AK AL AR AZ CA CO CT DE FL GA HI IA ID IL IN KS KY LA MA MD ME
#     MI MN MO MS MT NC ND NE NH NJ NM NV NY OH OK OR PA RI SC SD TN TX UT VA VT
#     WA WI WV WY"""
    countries = """AI AG AW BB BZ BR VG CO CR DM DO EC SV GD GT GY HN JM MX NI 
    PA PR KN LC VC TT VE  """
    # IEM quirk to have Iowa AWOS sites in its own labeled network
    networks = []
    for country in countries.split():
#        print(country)
        networks.append("%s__ASOS" % (country,))
#        print(networks)

    for network in networks:
        # Get metadata
        uri = (
            "https://mesonet.agron.iastate.edu/geojson/network/%s.geojson"
        ) % (network,)
#        print(uri)
        data = urlopen(uri)
        jdict = json.load(data)
        for site in jdict["features"]:
            stations.append(site["properties"]["sid"])
    return stations


def download_alldata():
    """An alternative method that fetches all available data.
    Service supports up to 24 hours worth of data at a time."""
    # timestamps in UTC to request data for
    startts = datetime.datetime(2021, 6, 16, 00)
    endts = datetime.datetime(2021, 6, 17, 00)
    interval = datetime.timedelta(hours=24)
#    get
    service = SERVICE + "data=all&tz=Etc/UTC&format=comma&latlon=yes&"
    print(service)
    now = startts
    while now < endts:
        thisurl = service
        thisurl += now.strftime("year1=%Y&month1=%m&day1=%d&")
        thisurl += (now + interval).strftime("year2=%Y&month2=%m&day2=%d&")
        print(thisurl)
        print("Downloading: %s" % (now,))
        data = download_data(thisurl)
        outfn = PATH + "%s.txt" % (now.strftime("%Y%m%d"),)
        with open(outfn, "w") as fh:
            fh.write(data)
        now += interval


def main(PATH):
    """Our main method"""
    # timestamps in UTC to request data for
    startts = datetime.datetime(2021, 6, 16, 00)
    endts = datetime.datetime(2021, 6, 17, 00)

    service = SERVICE + \
    "data=skyl1&data=skyl2&data=skyl3&tz=Etc/UTC&format=comma&latlon=yes&"

    service += startts.strftime("year1=%Y&month1=%m&day1=%d&")
    service += endts.strftime("year2=%Y&month2=%m&day2=%d&")

    # Two examples of how to specify a list of stations
    stations = get_stations_from_networks()
    # stations = get_stations_from_filelist("mystations.txt")
    for station in stations:
        uri = "%s&station=%s" % (service, station)
#        print(uri)
        print("Downloading: %s" % (station,))
        data = download_data(uri)
        outfn = PATH + "%s_%s_%s.txt" % (
            station,
            startts.strftime("%Y%m%d%H%M"),
            endts.strftime("%Y%m%d%H%M"),
        )
        out = open(outfn, "w")
        out.write(data)
        out.close()


if __name__ == "__main__":
#    download_alldata()
#    stations = get_stations_from_networks()
    PATH = '/home/alton/WRF_OUT/Sat_Reprojection/Data/METAR_CLDBASE_HGHT/20210616/'
    main(PATH)
