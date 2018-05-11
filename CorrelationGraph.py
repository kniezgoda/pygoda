#!/usr/bin/env python

# Calculates area mean for two variables and plots them against each other on a simple correlation plot

import os, glob
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pygoda as pg

start_year = 1
end_year = 20
start_month = 1
end_month = 12

varx = 'PRECT'
vary = 'PRECT_d180'

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


########################
### Code starts here ###
########################

dates = pg.h1dates(start_year, end_year, start_month, end_month)

var1 = []
var2 = []
files = []
used_dates = []
first = True
for n, d in enumerate(dates):
	# Find the file for this date
	f = glob.glob("*cam.h1."+d+"*.nc")
	
	#Progress tracker
	if n % (len(dates)/10) == 0:
		print str(round(float(n)/float(len(dates)) * 100)) + "% done..."
	
	if len(f) > 0: # If a file was found
		# Open the netcdf file for reading
		nc = Dataset(f[0], 'r')
		
		while first:
			print "Extracting lat and lon"
			# Only have to do this once
			lat = nc.variables['lat'][:]
			lon = nc.variables['lon'][:]
			box = pg.find_indices([southern_lat, northern_lat, left_lon, right_lon], lat, lon)
			first = False
		
		# Extract the data
		hold1 = nc.variables[varx][:,box[0]:box[1],box[2]:box[3]]
		hold2 = pg.prectd18O(nc)[box[0]:box[1],box[2]:box[3]]

		# Remove crazy values
		keep = hold2 > -50 
		hold2 = hold2[keep]
		hold1 = hold1[keep]
		keep = hold2 < 50 
		hold2 = hold2[keep]
		hold1 = hold1[keep]
		
		# Track the data
		var1.append(np.sum(hold1))
		var2.append(np.nanmean(hold2))
		files.append(f[0])
		used_dates.append(d)

var1 = np.array(var1)*1000*60*60*24
var2 = np.array(var2)

keep = var1 > 5 # only days when rain rate was greater than 5 mm/day
var2 = var2[keep]
var1 = var1[keep]

plt.scatter(var1, var2)
plt.show()
