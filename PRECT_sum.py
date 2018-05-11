#!/usr/bin/env python

# Computes precipitation sums over a spatial region and plots a timeseries

import os, glob
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def find_indices(box, lats, lons):
    r = []
    for N, b in enumerate(box):
        if N < 2:
            r.append(abs(lats-b).tolist().index(min(abs(lats-b))))
        if N >= 2:
            r.append(abs(lons-b).tolist().index(min(abs(lons-b))))
    return r


start_year = 20
end_year = 20
start_month = 1
end_month = 12

# Contruct dates to look for
# All months will have 31 days, but that is ok because glob will ignore any files it doesn't find
years = []
if end_year < 10:
	years = years + ["000" + str(y) for y in range(start_year,(end_year+1))] 
elif end_year >= 10:
	if start_year < 10:
		years = years + ["000" + str(y) for y in range(start_year,10)]
		years = years + ["00" + str(y) for y in range(10,(end_year+1))]
	if start_year > 10:
		years = years + ["00" + str(y) for y in range(start_year,(end_year+1))]


months = []
if end_month < 10 :
	months = months + ["0" + str(y) for y in range(start_year,(end_month+1))] 
elif end_month >= 10:
	if start_month < 10:
		months = months + ["0" + str(y) for y in range(start_month,10)]
		months = months + [str(y) for y in range(10,(end_month+1))]
	if start_month > 10:
		months = months + [str(y) for y in range(start_month,(end_month+1))]


dates = []
for y in years:
	for m in months:
		dates = dates + [y + "-" + m + "-0" + str(d) for d in range(1,10)]
		dates = dates + [y + "-" + m + "-" + str(d) for d in range(10,32)]


print "The following dates will be used as search keywords for files:\n"
print dates

# Set region
region_name = "Sahara"
southern_lat = 25
northern_lat = 30
left_lon = -20
right_lon = 40

# region_name = "IndianMonsoon"
# southern_lat = 20
# northern_lat = 30
# left_lon = 70
# right_lon = 90

PRECT = []
files = []
used_dates = []
for n, d in enumerate(dates):
	# Find the file for this date
	f = glob.glob("*cam.h1."+d+"*.nc")
	#Progress tracker
	if n % (len(dates)/10) == 0:
		print str(round(float(n)/float(len(dates)) * 100)) + "% done..."
	if len(f) > 0: # If a file was found
		# Open the netcdf file for reading
		nc = Dataset(f[0], 'r')
		# Extract lat and lon
		lat = nc.variables['lat'][:]
		lon = nc.variables['lon'][:]
		# Find the indices for the region
		box = find_indices([southern_lat, northern_lat, left_lon, right_lon], lat, lon)		
		# Extract prect data for the region
		prect = nc.variables['PRECT'][:,box[0]:box[1],box[2]:box[3]]		
		# Keep track of prect and files
		PRECT.append(np.sum(prect) * 1000 * 60 * 60 * 24) # Convert to mm/day
		files.append(f[0])
		used_dates.append(d)


plt.figure()
ts = pd.Series(PRECT, index=used_dates)
ts.plot()
plt.show()
