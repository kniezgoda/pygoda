#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python
from pygoda import camgoda, camdates, findClimoFile
import numpy as np
from scipy.stats import pearsonr as regress
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as bm
import argparse


'''
Author: Kyle Niezgoda, April 28, 2019

This code produces maps of correlation arrays between the local variable (-lv) at a point (-latlon) 
and the field variable (-fv) over the globe for a specified time (-dates)
'''
#########################
# Read in the arguments #
#########################

parser = argparse.ArgumentParser()
parser.add_argument('-cdir', '--control_directory', dest = 'cdir', default = 'fc5.2deg.wiso.piControl_kn028/atm/hist')
parser.add_argument('-tdir', '--test_directory', dest = 'tdir', default = 'fc5.2deg.wiso.mh6ka_kn032/atm/hist')
parser.add_argument('-r', '--region', dest = 'region', default = None)
parser.add_argument('-grep', dest = 'grep', default = 'cam.h0')
parser.add_argument('-lats', dest = 'lats', nargs = 2, default = [-90,90])
parser.add_argument('-lons', dest = 'lons', nargs = 2, default = [0,360])
parser.add_argument('-v', '--variables', dest = 'variables', nargs= 2, default = ["PRECT", "PRECT_d18O"])
parser.add_argument('-years', dest = 'years', nargs = 2, default = [10,30])
parser.add_argument('-months', dest = 'months', nargs = '*', default = [1,2,3,4,5,6,7,8,9,10,11,12])
parser.add_argument('-alpha', '--alpha', dest = 'alpha', default = .05)
parser.add_argument('-stiple', '--stiple', dest = 'stiple', action = "store_true")
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
cdir = ARGS.cdir
tdir = ARGS.tdir
start, end = [int(x) for x in ARGS.years]
months = [int(m) for m in ARGS.months]
v1, v2 = [v for v in ARGS.variables]
grep = ARGS.grep
region = ARGS.region
savefig = ARGS.savefig
showfig = ARGS.showfig
alpha = float(ARGS.alpha)
stiple = ARGS.stiple
if ARGS.developer_mode:
    print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
    savfig = False
    showfig = True


#############
# Main code #
#############

# Make the date array
dates = camdates(start, end, months)

# Read in each file, extract the local variable at the point and field variable on the globe, track them in a master array
# Correlate the local variable 1-d array to each grid cell of the field variable array
# Return the correlation coefficients 

# Set the lat bounds
# Default global tropics
region_name = "Box"
southern_lat, northern_lat = [int(l) for l in ARGS.lats]
left_lon, right_lon = [int(l) for l in ARGS.lons]

# Global
if (region == "GT") | (region == "GlobalTropics"):
	region_name = "GlobalTropics"
	southern_lat = -50
	northern_lat = 50
	left_lon = 0
	right_lon = 360
# Indian monsoon
if (region == "IM") | (region == "IndianMonsoon"):
	region_name = "IndianMonsoon"
	southern_lat = -5
	northern_lat = 45
	left_lon = 40
	right_lon = 110

# Greater maritime continent
if (region == "MC") | (region == "MaritimeContinent"):
	region_name = "MaritimeContinent"
	southern_lat = -20
	northern_lat = 20
	left_lon = 80
	right_lon = 160

# North Africa
if (region == "NA") | (region == "NorthAfrica"):
	region_name = "NorthAfrica"
	southern_lat = -20
	northern_lat = 45
	left_lon = -30
	right_lon = 70

#Central Equatorial Pacific to Western eq. Atlantic
if (region == "EP") | (region == "TropicalOceans"):
	region_name = "TropicalOceans"
	southern_lat = -25
	northern_lat = 25
	left_lon = 180
	right_lon = 355

#South America and the Carribean 
if (region == "SA") | (region == "SouthAmerica"):
	region_name = "SouthAmerica"
	southern_lat = -30
	northern_lat = 30
	left_lon = -120
	right_lon = -20

# Four proxy region in the warm pool
if (region == "4P") | (region == "FourProxies"):
	region_name = "FourProxies"
	southern_lat = -25
	northern_lat = 30
	left_lon = 85
	right_lon = 160

box = (southern_lat, northern_lat, left_lon, right_lon)

tloc_var_master = []
cloc_var_master = []

for n, d in enumerate(dates):
	# Find the file
	cpath, cfilename = findClimoFile(cdir + "/" + '*' + grep + "*" + d + '*')
	tpath, tfilename = findClimoFile(tdir + "/" + '*' + grep + "*" + d + '*')
	print tfilename
	print cfilename
	
	# Open the file
	cnc = camgoda(cpath)
	tnc = camgoda(tpath)

	cnc.ExtractData(v1, box)
	cv1_hold = np.expand_dims(cnc.data, 0)
	cnc.ExtractData(v2, box)
	cv2_hold = np.expand_dims(cnc.data, 0)

	tnc.ExtractData(v1, box)
	tv1_hold = np.expand_dims(tnc.data, 0)
	tnc.ExtractData(v2, box)
	tv2_hold = np.expand_dims(tnc.data, 0)


	if n == 0:
		cv1_master = cv1_hold
		cv2_master = cv2_hold
		tv1_master = tv1_hold
		tv2_master = tv2_hold

	else:
		cv1_master = np.concatenate((cv1_master, cv1_hold), axis = 0)
		cv2_master = np.concatenate((cv2_master, cv2_hold), axis = 0)

		tv1_master = np.concatenate((tv1_master, tv1_hold), axis = 0)
		tv2_master = np.concatenate((tv2_master, tv2_hold), axis = 0)

print "All data extracted, computing correlations..."

ntime, nlats, nlons = cv1_master.shape
ccorrArray = np.zeros((nlats, nlons))
tcorrArray = np.zeros((nlats, nlons))
for i in nlats:
	print i
	for j in nlons:
		print j
		ccorrArray[i,j] = np.corrcoef(cv1_master[:i,j], cv2_master[:i,j])[0,1]
		tcorrArray[i,j] = np.corrcoef(tv1_master[:i,j], tv2_master[:i,j])[0,1]

print ccorrArray