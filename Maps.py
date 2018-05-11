#!/Users/kyleniezgoda/anaconda/bin/python
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
from pygoda import ncgoda, findClimoFile, niceClev
root = os.getcwd()
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
	"PRECC" : 0,         \
	"PRECL" : 0,         \
	"PRECSC" : 0,        \
	"PRECSL" : 0,        \
	"PRECST" : 0,        \
	"QFLX" : 1,          \
	"QFLX_d18O" : 0,	 \
	"QFLX_dD" : 0,		 \
	"fluxDelta" : 0,	 \
	"P_E" : 1,			 \
	"CLDHGH" : 0,        \
	"PS" : 0,            \
	"PSL" : 0,  		 \
	"Column_d18OV" : 01, \
	"Column_dDV" : 0,	 \
	"PRECT_d18O" : 0,    \
	"PRECT_dD" : 0,      \
	"PRECT_dxs" : 0      \
	}

use3d = {                \
	"Q850" : 0,          \
	"Q500" : 0,          \
	"Q200" : 0,          \
	"V850" : 0,          \
	"V500" : 0,          \
	"V200" : 0,          \
	"VT850" : 0,         \
	"VT500" : 0,         \
	"VT200" : 0,         \
	"VQ850" : 0,         \
	"VQ500" : 0,         \
	"VQ200" : 0,         \
	"U850" : 0,          \
	"U500" : 0,          \
	"U200" : 0,          \
	"UT850" : 0,         \
	"UT500" : 0,         \
	"UT200" : 0,         \
	"UQ850" : 0,         \
	"UQ500" : 0,         \
	"UQ200" : 0,         \
	"T850" : 0,          \
	"T500" : 0,          \
	"T200" : 0,          \
	"OMEGA850" : 0,      \
	"OMEGA500" : 0,      \
	"OMEGA200" : 0,      \
	"Z3850" : 0,         \
	"Z3500" : 0,         \
	"Z3200" : 0,         \
	"dDV850" : 0,        \
	"dDV500" : 0,        \
	"dDV200" : 0,        \
	"d18OV850" : 0,      \
	"d18OV500" : 0,      \
	"d18OV200" : 0,      \
	"d18OV999" : 0,		 \
	"dxsV850" : 0,       \
	"dxsV500" : 0,       \
	"dxsV200" : 0,       \
	"psi850" : 0,        \
	"psi500" : 0,        \
	"psi200" : 0,        \
	"UQ_d18O850" : 0,	 \
	"UQ_dD850" : 0,		 \
	"VQ_d18O850" : 0, 	 \
	"VQ_dD850" : 0,		 \
	"UQ_d18O500" : 0,    \
    "UQ_dD500" : 0,      \
    "VQ_d18O500" : 0,    \
	"VQ_dD500" : 0,		 \
	"UQ_d18O200" : 0,    \
    "UQ_dD200" : 0,      \
    "VQ_d18O200" : 0,    \
	"VQ_dD200" : 0,		 \
	"RH850" : 0, 		 \
	"RH500" : 0,		 \
	"RH200" : 0,		 \
	"RH999" : 0			 \
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
mkdir = True
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savfig = False
	showfig = True
	mkdir = False

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

if mkdir:
	# Create maps directory is it doesn't exist
	if not os.path.exists("Maps"):
		os.mkdir("Maps")
		print "Created directory " + "Maps"

	# Create the region directory if it doesn't already exist
	if not os.path.exists("Maps/" + region_name):
		os.mkdir("Maps/" + region_name)
		print "Created directory " + "Maps/" + region_name

	# Create season directory inside region directory
	if not os.path.exists("Maps/" + region_name + "/" + season):
		os.mkdir("Maps/" + region_name + "/" + season)
		print "Created directory " + "Maps/" + region_name + "/" + season

# Look for the climo files in the root directory
print "\nLooking for control " + season + " files in " + controldir + "..."
controldatafname, controlfn = findClimoFile("*" + season + "*", root+"/"+controldir)
if not controldatafname:
	sys.exit()
else:
	print "Found file " + controlfn
print "\nLooking for test " + season + " files in " + testdir + "..."
testdatafname, testfn = findClimoFile("*" + season + "*", root+"/"+testdir)
if not testdatafname:
	sys.exit()
else:
	print "Found file " + testfn

# Read the data
controldata = ncgoda(controldatafname)
testdata = ncgoda(testdatafname)

# Set the boxlat and boxlon
box = (southern_lat, northern_lat, left_lon, right_lon)
controldata.variable("T", box, setData = False) # this sets self.boxlon and self.boxlat

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
if mkdir:
	os.chdir("Maps/" + region_name + "/" + season)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------- Map creation ---------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for variable in use2d:

	if not use2d.get(variable):
		continue

	print("\nPlotting " + season + " " + variable + " data...")

	if variable == "PRECT_d18O":
		testdata.PRECT_d18O(box)
		controldata.PRECT_d18O(box)
	
	elif variable == "PRECT_dD":
		testdata.PRECT_dD(box)
		controldata.PRECT_dD(box)
	
	elif variable == "PRECT_dxs":
		testdata.PRECT_dxs(box)
		controldata.PRECT_dxs(box)
	
	elif variable == "QFLX_d18O":
		testdata.QFLX_d18O(box)
		controldata.QFLX_d18O(box)
	
	elif variable == "QFLX_dD":
		testdata.QFLX_dD(box)
		controldata.QFLX_dD(box)
	
	elif variable == "fluxDelta":
		testdata.fluxDelta(box)
		controldata.fluxDelta(box)
	
	elif variable == "Column_d18OV":
		#test
		testdata.variable('H2OV', box)
		denom = testdata.columnSum(box)
		testdata.variable('H218OV', box)
		num = testdata.columnSum(box)
		testdata.data = (num/denom - 1) * 1000
		#control
		controldata.variable('H2OV', box)
		denom = controldata.columnSum(box)
		controldata.variable('H218OV', box)
		num = controldata.columnSum(box)
		controldata.data = (num/denom - 1) * 1000
	
	elif variable == "Column_dDV":
		#test
		testdata.variable('H2OV', box)
		denom = testdata.columnSum(box)
		testdata.variable('HDOV', box)
		num = testdata.columnSum(box)
		testdata.data = (num/denom - 1) * 1000
		#control
		controldata.variable('H2OV', box)
		denom = controldata.columnSum(box)
		controldata.variable('HDOV', box)
		num = controldata.columnSum(box)
		controldata.data = (num/denom - 1) * 1000

	elif variable == "P_E":
		#test
		testdata.data = (testdata.variable('PRECT', box, math = False)*1000 - testdata.variable('QFLX', box, math = False)) * 60 * 60 * 24
		controldata.data = (controldata.variable('PRECT', box, math = False)*1000 - controldata.variable('QFLX', box, math = False)) * 60 * 60 * 24
		testdata.units = "kg/m2/day"
		testdata.long_name = "Moisture flux due to advection"
		controldata.units = "kg/m2/day"
		controldata.long_name = "Moisture flux due to advection"

	else:
		try:
			testdata.variable(variable, box, verb = True)
			controldata.variable(variable, box)
		except KeyError:
			print "Not able to plot variable " + variable + "...\nSkipping this variable."
			continue

	fig = plt.figure()

	testdata.prep_map(season, region)
	controldata.prep_map(season, region)
	testdata.clevs = np.linspace(-3,3,13) # temp
	controldata.clevs = np.linspace(-3,3,13) # temp
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

	diffclev = niceClev(testdata.data - controldata.data)
	diffclev = np.linspace(-30,30,21)
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

	print("\nPlotting " + season + " " + variable + " data...")

	var = variable[:-3]
	pressure = int(variable[-3:]) * 100

	if var == "d18OV":
		testdata.d18OV(box)
		controldata.d18OV(box)
	elif var == "dDV":
		testdata.dDV(box)
		controldata.dDV(box)
	elif var == "dxsV":
		testdata.dxsV(box)
		controldata.dxsV(box)
	elif var == "psi":
		testdata.psi(box)
		controldata.psi(box)
	elif var == "RH":
		testdata.RH(box)
		controldata.RH(box)
	elif var == "VQ_d18O":
		testdata.VQ_d18O(box)
		controldata.VQ_d18O(box)
	elif var == "VQ_dD":
		testdata.VQ_dD(box)
		controldata.VQ_dD(box)
	elif var == "UQ_d18O":
		testdata.UQ_d18O(box)
		controldata.UQ_d18O(box)
	elif var == "UQ_dD":
		testdata.UQ_dD(box)
		controldata.UQ_dD(box)
	elif var == "QFLX_d18O":
		testdata.QFLX_d18O(box)
		controldata.QFLX_d18O(box)
	else:
		try:
			testdata.variable(var, box, verb = True)
			controldata.variable(var, box)
		except KeyError:
			print "Not able to plot variable " + variable + "...\nSkipping this variable."
			continue

	testdata.data = testdata.isobar(pressure)
	controldata.data = controldata.isobar(pressure)
	
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
	diffclev = niceClev(testdata.data - controldata.data)
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
