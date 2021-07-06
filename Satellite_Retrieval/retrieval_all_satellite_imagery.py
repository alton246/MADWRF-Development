import os           # Miscellaneous operating system $
import subprocess   # The subprocess module allows yo$
import platform           # Access to underlying plat$
import sys



def GetSatCLDMASK(year, julian_day, outdir):
    # Get the O.S.
    osystem = platform.system()
    if osystem == "Windows":
        extension = '.exe'
    else:
        extension = ' '

    print ("GOES-R Cloud Data Downloader")
    print("Downloading Cloud Mask Data")
    # Desired Data
    BUCKET = 'noaa-goes16'    # For GOES-R the buckets ar$
    PRODUCT = 'ABI-L2-ACMF'  # Choose from ['ABI-L1b-RadC$
    #PRODUCT = 'ABI-L2-AODC'
    YEAR =  year         # Choose from ['2017', '201$
    JULIAN_DAY = julian_day        # Data available after juli$
    #HOUR = '15'               # Choose from 00 to 23


    for hr in range(0,24):
        print ("Current year, julian day and hour based on your local machine:")
        print("YEAR: ", YEAR)
        print("JULIAN DAY: ", JULIAN_DAY)

        if len(str(hr)) < 2:
            HOUR = str(hr).zfill(2)
        else:
            HOUR = str(hr)
        print("HOUR (UTC): ", HOUR)

        OUTDIR = outdir + str(HOUR) + '/'
        files = subprocess.check_output('rclone' + extension + " " + 'ls publicAWS:' + BUCKET + "/" + PRODUCT + "/" + YEAR + "/" + JULIAN_DAY + "/" + HOUR + "/", shell=True)

        files = files.decode()
        files = files.split('\n')
        files.remove('')

        files = [i.split(" ")[-1] for i in files]

        if not files:
            print("No files available yet... Exiting script")
            sys.exit()

        for i in range(len(files)):
    #    print(files[i])
            os.system('rclone' + extension + " " + 'copy publicAWS:' + BUCKET + "/" + PRODUCT + "/" + YEAR + "/" + JULIAN_DAY + "/" + str(HOUR) + "/" + str(files[i]) + " " + OUTDIR)
    print ("Finish Downloading Cloud Mask Data!")


def GetSatCTOPH(year, julian_day, outdir):
    # Get the O.S.
    osystem = platform.system()
    if osystem == "Windows":
        extension = '.exe'
    else:
        extension = ' '

    print ("GOES-R Cloud Data Downloader")
    print ("Downloading Cloud Top Height Data")
    # Desired Data
    BUCKET = 'noaa-goes16'    # For GOES-R the buckets are: ['noaa-goes16', 'noaa-goes17']
    PRODUCT = 'ABI-L2-ACHAF'  # Choose from ['ABI-L1b-RadC', 'ABI-L1b-RadF', 'ABI-L1b-RadM', 'ABI-L2-CMIPC', 'ABI-L2$
    
    YEAR = year             # Choose from ['2017', '2018', '2019']
    JULIAN_DAY = julian_day        # Data available after julian day 283, 2017 (October 10, 2017)

    for hr in range(0,24):
        print ("Current year, julian day and hour based on your local machine:")
        print("YEAR: ", YEAR)
        print("JULIAN DAY: ", JULIAN_DAY)

    if len(str(hr)) < 2:
        HOUR = str(hr).zfill(2)
    else:
	    HOUR = str(hr)
    print("HOUR (UTC): ", HOUR)

    OUTDIR = outdir + str(HOUR) + '/'

    files = subprocess.check_output('rclone' + extension + " " + 'ls publicAWS:' + BUCKET + "/" + PRODUCT + "/" + YEAR + "/" + JULIAN_DAY + "/" + HOUR + "/", shell=True)

    files = files.decode()
    ## Split files based on the new line and remove the empty item at the end.
    files = files.split('\n')
    files.remove('')
    ## Get only the file names for an specific channel

#    # Get only the file names, without the file sizes
    files = [i.split(" ")[-1] for i in files]


    for i in range(len(files)):
        os.system('rclone' + extension + " " + 'copy publicAWS:' + BUCKET + "/" + PRODUCT + "/" + YEAR + "/" + JULIAN_DAY + "/" + str(HOUR) + "/" + str(files[i]) + " " + OUTDIR)

    print ("Finish Downloading Cloud Top Height Data!")

    

def GetSatBRTEMP(year, julian_day, outdir):
    # Get the O.S.
    osystem = platform.system()
    if osystem == "Windows":
        extension = '.exe'
    else:
        extension = ' '

    print("GOES-R Cloud Data Downloader")
    print("Downloading Brightness Temperature Data")

    # Desired Data
    BUCKET = 'noaa-goes16'    # For GOES-R the buckets are: ['noaa-goes16', 'noaa-goes17']
    PRODUCT = 'ABI-L2-CMIPF'  # Choose from ['ABI-L1b-RadC', 'ABI-L1b-RadF', 'ABI-L1b-RadM', 'ABI-L2-CMIPC', 'ABI-L2$

    YEAR = year             # Choose from ['2017', '2018', '2019']
    JULIAN_DAY = julian_day        # Data available after julian day 283, 2017 (October 10, 2017)

    for hr in range(0,24):
        print ("Current year, julian day and hour based on your local machine:")
        print("YEAR: ", YEAR)
        print("JULIAN DAY: ", JULIAN_DAY)

    if len(str(hr)) < 2:
        HOUR = str(hr).zfill(2)
    else:
	    HOUR = str(hr)
    print("HOUR (UTC): ", HOUR)

    CHANNEL = 'C13'           # Choose from ['C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12 C13 C14 C15 C16']
    OUTDIR = outdir + str(HOUR) + '/'

    files = subprocess.check_output('rclone' + extension + " " + 'ls publicAWS:' + BUCKET + "/" + PRODUCT + "/" + YEAR + "/" + JULIAN_DAY + "/" + str(HOUR) + "/", shell=True)

    files = files.decode()
    ## Split files based on the new line and remove the empty item at the end.
    files = files.split('\n')
    files.remove('')
    ## Get only the file names for an specific channel
    files = [x for x in files if CHANNEL in x ]
    # Get only the file names, without the file sizes
    files = [i.split(" ")[-1] for i in files]

    if not files:
        print("No files available yet... Exiting script")
        sys.exit()

    print ("Downloading files ...")

    for i in range(len(files)):
#    print(files[i])
        os.system('rclone' + extension + " " + 'copy publicAWS:' + BUCKET + "/" + PRODUCT + "/" + YEAR + "/" + JULIAN_DAY + "/" + str(HOUR) + "/" + str(files[i]) + " " + OUTDIR)

print("Finish Downloading Brightness Temperature Data!")








#############################################################

#############################################################
JULIAN_DAY = '168'
YEAR = '2021'
OUTDIR = '/home/alton/Satellie_Imagery/Data/Channel13/20210616/'






# GetSatCLDMASK(YEAR, JULIAN_DAY, OUTDIR)
# GetSatCTOPH(YEAR, JULIAN_DAY, OUTDIR)
GetSatBRTEMP(YEAR, JULIAN_DAY, OUTDIR)