#!/home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python
'''
This code creates test, control, and difference maps for a collection of netcdf variables from iCAM5.
The three maps are plotted to a figure and saved to a user-entered directory.

Author: Kyle Niezgoda
Date last updated: June 14, 2017

####################
### Things to do ###
####################
1)
Figure out if there is a smart way to make color bars better. Sometimes outliers really stick out and 
make the maps look ugly. If I could find those outliers reliably and adjust the colorbar accordingly, that'd 
be nice and save a lot of manually hard-coding in the colorbar limits.

'''

import os, sys, glob, argparse
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap as bm 
from pygoda import ncgoda, findClimoFile, zeroCenterClev
root = os.getcwd()
if not os.path.exists("Maps"):
	os.mkdir("Maps")
	print "Created directory " + "Maps"
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
Probably will want to put the following variables in a text file
or some file structure and then have the user supply the file in a 
command line flag or argument.
'''
#------------- which variables to plot? (1 -> True, 0 -> False)
# This is all over-ridden if -v flag is supplied at the command line
use2d = {                \
	"PRECT" : 1,         \
	"PRECC" : 1,         \
	"PRECL" : 1,         \
	"PRECSC" : 1,        \
	"PRECSL" : 1,        \
	"PRECST" : 1,        \
	"QFLX" : 1,          \
	"CLDHGH" : 1,        \
	"PS" : 1,            \
	"PSL" : 1,           \
	"PRECT_d18O" : 1,    \
	"PRECT_dD" : 1,      \
	"PRECT_dxs" : 1      \
	}

use3d = {                \
	"Q850" : 1,          \
	"Q500" : 1,          \
	"Q200" : 1,          \
	"V850" : 1,          \
	"V500" : 1,          \
	"V200" : 1,          \
	"VT850" : 1,         \
	"VT500" : 1,         \
	"VT200" : 1,         \
	"VQ850" : 1,         \
	"VQ500" : 1,         \
	"VQ200" : 1,         \
	"U850" : 1,          \
	"U500" : 1,          \
	"U200" : 1,          \
	"UT850" : 1,         \
	"UT500" : 1,         \
	"UT200" : 1,         \
	"UQ850" : 1,         \
	"UQ500" : 1,         \
	"UQ200" : 1,         \
	"T850" : 1,          \
	"T500" : 1,          \
	"T200" : 1,          \
	"OMEGA850" : 1,      \
	"OMEGA500" : 1,      \
	"OMEGA200" : 1,      \
	"Z3850" : 1,         \
	"Z3500" : 1,         \
	"Z3200" : 1,         \
	"dDV850" : 1,        \
	"dDV500" : 1,        \
	"dDV200" : 1,        \
	"d18OV850" : 1,      \
	"d18OV500" : 1,      \
	"d18OV200" : 1,      \
	"dxsV850" : 1,       \
	"dxsV500" : 1,       \
	"dxsV200" : 1        \
	}

#------------- File type for saved image. Don't include any periods, just the suffix
ftype = 'ps'

#------------- Custom subplot title names 
# these are the names of the subplots, not the main figure --- the main figure always has the same name (the variable long name from the netcdf file)
# Set to None (no quotes) for default subplot names, the name of the file used
# Default 
# user_test_plot_title = None 
# user_control_plot_title = None 
# Regular 
user_test_plot_title = "Test: Mid-Holocene (6kya)"
user_control_plot_title = "Control: Pre-industrial (1850)"

##############################################
###### Don't change anything below here ######
##############################################


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#################################################
####### Parse arguments from command line #######
#################################################

'''
6 optional args:
-r (--region) REGION : sets the region, default to '' (global tropics)
-s (--season) SEASON : sets season, defaults to 'ANN'
-tdir (test_directory) TDIR : sets the directory to look for test files in, defaults to *_kn003 directory
-cdir (control_directory) CDIR : see above, for control, defaults to *_kn002 directory
-show : sets showfig to True, will print all plots to screen. Default is not to show figs
-nosave : sets savefig to False, will not save figures to current directory. Default is to save images to current directory
'''

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--region', dest = 'region', default = '')
parser.add_argument('-s', '--season', dest = 'season', default = 'ANN')
parser.add_argument('-cdir', '--control_directory', dest = 'controldir', default = "F.C5.2deg.wiso.defaultSSTICE_kn002")
parser.add_argument('-tdir', '--test_directory', dest = 'testdir', default = "F.C5.2deg.wiso.obs6kSST_kn003")
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variable', nargs= "*", default = None)
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
region = ARGS.region
print "Region is " + region
season = ARGS.season
print "Season is " + season
testdir = ARGS.testdir
controldir = ARGS.controldir
savefig = ARGS.savefig
showfig = ARGS.showfig
variable = ARGS.variable
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved and all plots will be printed to the screen."
	savfig = False
	showfig = True

# Set the new variable list if variable is not None
if variable is not None:
	use2d = {u : 0 for u in use2d}
	use3d = {u : 0 for u in use3d}
	for v in variable:
		if v in use2d:
			use2d[v] = 1
		elif v in use3d:
			use3d[v] = 1
		else:
			print "\nVariable " + v + " not in master variable list.\nWill not plot this variable."


# Set the lat bounds
# Default global tropics
region_name = "GlobalTropics"
southern_lat = -50
northern_lat = 50
left_lon = 0
right_lon = 355

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

# Create the region directoy if it doesn't already exist
if not os.path.exists("Maps/" + region_name):
	os.mkdir("Maps/" + region_name)
	print "Created directory " + "Maps/" + region_name

# Create season directory inside region directory
if not os.path.exists("Maps/" + region_name + "/" + season):
	os.mkdir("Maps/" + region_name + "/" + season)
	print "Created directory " + "Maps/" + region_name + "/" + season

# Look for the climo files in the root directory
print "\nLooking for control " + season + " files in " + controldir + "..."
controldatafname, controlfn = findClimoFile(season, root+"/"+controldir)
if not controldatafname:
	sys.exit()
else:
	print "Found file " + controlfn
print "\nLooking for test " + season + " files in " + testdir + "..."
testdatafname, testfn = findClimoFile(season, root+"/"+testdir)
if not testdatafname:
	sys.exit()
else:
	print "Found file " + testfn

# Read the data
controldata = ncgoda(controldatafname)
testdata = ncgoda(testdatafname)

# Set the boxlat and boxlon
box = (southern_lat, northern_lat, left_lon, right_lon)
controldata.variable("T", box) # this sets self.boxlon and self.boxlat

# Create bm coords from region bounds
# bm lonitude coords need to be 0 < coord < 360 
bmlon, bmlat = np.meshgrid(controldata.boxlon, controldata.boxlat)

# Reset the lat and lon bounds so that maps don't show grey areas 
southern_lat, northern_lat = np.array(controldata.boxlat)[[0,-1]]
# Change lons to be negative is 180 < lon < 360 because that's how bm works for 'cea' projection
left_lon, right_lon = np.array(controldata.boxlon)[[0,-1]]
if 0 in controldata.boxlon[1:-2]: # if we cross the gml
	left_lon = controldata.boxlon[0]-360


g = -9.8 # gravitational constant

# Change into figure directory (root/region/season/) for image creation
os.chdir("Maps/" + region_name + "/" + season)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------- Map creation ---------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for variable in use2d:

	if not use2d.get(variable):
		continue

	if variable == "PRECT_d18O":
		testdata.PRECT_d18O(box)
		controldata.PRECT_d18O(box)
	elif variable == "PRECT_dD":
		testdata.PRECT_dD(box)
		controldata.PRECT_dD(box)
	elif variable == "PRECT_dxs":
		testdata.PRECT_dxs(box)
		controldata.PRECT_dxs(box)
	else:
		try:
			testdata.variable(variable, box, verb = True)
			controldata.variable(variable, box)
		except KeyError:
			print "Not able to plot variable " + variable + "...\nSkipping this variable."
			continue

	print("\nPlotting " + season + " " + variable + " data...")
	fig = plt.figure()

	testdata.prep_map(season, region)
	controldata.prep_map(season, region)

	# test data
	plt.subplot(3,1,1)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	cs = m.contourf(bmlon, bmlat, testdata.data, testdata.clevs, shading = 'flat', latlon = True, cmap=testdata.cmap)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label(testdata.units, fontsize = 8)
	test_title = user_test_plot_title + ", " + season
	if user_test_plot_title is None:
		test_title = "test: " + testfn
	plt.title(test_title, fontsize = 8)
	
	# control data
	plt.subplot(3,1,2)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	cs = m.contourf(bmlon, bmlat, controldata.data, controldata.clevs, shading = 'flat', latlon = True, cmap=controldata.cmap)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label(testdata.units, fontsize = 8)
	control_title = user_control_plot_title + ", " + season
	if user_control_plot_title is None:
		control_title = "control: " + controlfn
	plt.title(control_title, fontsize = 8)
	
	# difference data 
	plt.subplot(3,1,3)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')

	diffclev = zeroCenterClev(testdata.data - controldata.data)
	cs = m.contourf(bmlon, bmlat, testdata.data - controldata.data, diffclev, shading = 'flat', latlon = True, cmap=controldata.diffcmap)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label("difference in " + testdata.units, fontsize = 8)
	plt.title("test - control difference", fontsize = 8)		
	
	# Make things pretty
	fig.suptitle(testdata.long_name, fontweight = 'bold', fontsize = 14)
	plt.subplots_adjust(hspace = .2)
	if savefig:
		plt.savefig(variable + "_" + testfn + "-" + controlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + variable + "_" + testfn + "-" + controlfn + "." + ftype + '\n')
	if showfig:
		plt.show()
	plt.clf()
	plt.cla()
	plt.close(fig)

#---------------------------------#

for variable in use3d:

	if not use3d.get(variable):
		continue

	var = variable[:-3]
	pressure = int(variable[-3:]) * 100
	
	if var == "d18OV":
		testdata.d18OV(pressure, box)
		controldata.d18OV(pressure, box)
	elif var == "dDV":
		testdata.dDV(pressure, box)
		controldata.dDV(pressure, box)
	elif var == "dxsV":
		testdata.dxsV(pressure, box)
		controldata.dxsV(pressure, box)
	else:
		try:
			testdata.isobar(var, pressure, box, verb = True)
			controldata.isobar(var, pressure, box)
		except KeyError:
			print "Not able to plot variable " + variable + "...\nSkipping this variable."
			continue

	print("\nPlotting " + season + " " + variable + " data...")
	fig = plt.figure()

	testdata.prep_map(season, region)
	controldata.prep_map(season, region)

	# test data
	plt.subplot(3,1,1)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	cs = m.contourf(bmlon, bmlat, testdata.data, 19, shading = 'flat', latlon = True, cmap=testdata.cmap)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label(testdata.units, fontsize = 8)
	test_title = user_test_plot_title + ", " + season
	if user_test_plot_title is None:
		test_title = "test: " + testfn
	plt.title(test_title, fontsize = 8)
	
	# control data
	plt.subplot(3,1,2)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	cs = m.contourf(bmlon, bmlat, controldata.data, 19, shading = 'flat', latlon = True, cmap=controldata.cmap)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label(controldata.units, fontsize = 8)
	control_title = user_control_plot_title + ", " + season
	if user_control_plot_title is None:
		control_title = "control: " + controlfn
	plt.title(control_title, fontsize = 8)
	
	# difference data 
	plt.subplot(3,1,3)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	diffclev = zeroCenterClev(testdata.data - controldata.data)
	cs = m.contourf(bmlon, bmlat, testdata.data - controldata.data, diffclev, shading = 'flat', latlon = True, cmap=controldata.diffcmap)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label("difference in " + testdata.units, fontsize = 8)
	plt.title("test - control difference", fontsize = 8)		
	
	# Make things pretty
	fig.suptitle(testdata.long_name + " at " + str(pressure/100) + " mb", fontweight = 'bold', fontsize = 14)
	plt.subplots_adjust(hspace = .2)
	if savefig:
		plt.savefig(variable + "_" + testfn + "-" + controlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + variable + "_" + testfn + "-" + controlfn + "." + ftype + '\n')
	if showfig:
		plt.show()
	plt.clf()
	plt.cla()
	plt.close(fig)
