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
from pygoda import popgoda, findClimoFile, niceClev, RegularClev
root = os.getcwd()
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
Probably will want to put the following variables in a text file
or some file structure and then have the user supply the file in a 
command line flag or argument.
'''

#------------- Custom subplot title names 
# these are the names of the subplots, not the main figure --- the main figure always has the same name (the variable long name from the netcdf file)
# Set to None (no quotes) for default subplot names, the name of the file used
# Default 
user_test_plot_title = None 
user_control_plot_title = None
# Regular 
#user_test_plot_title = "Test: Mid-Holocene (6kya)"
#user_control_plot_title = "Control: Pre-industrial (1850)"
#user_test_plot_title = "Test: Estimated MH"
#user_control_plot_title = "Control: Simulated MH"


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
parser.add_argument('-cdir', '--control_directory', dest = 'controldir', default = '.')
parser.add_argument('-tdir', '--test_directory', dest = 'testdir', default = '.')
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variable', nargs= "*", default = None)
parser.add_argument('-t', '--test', dest = 'testdatafname', default = None)
parser.add_argument('-c', '--control', dest = 'controldatafname', default = None)
parser.add_argument('-clev', dest = 'clev', type = float, nargs = 3, default = None)
parser.add_argument('-diffclev', dest = 'diffclev',type = float, nargs = 3, default = None)
parser.add_argument('-ft', '--file_type', dest = 'file_type', default = 'ps')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
region = ARGS.region
print "Region is " + region
testdatafname = ARGS.testdatafname
controldatafname = ARGS.controldatafname
findFile = True
if (testdatafname is not None) & (controldatafname is not None):
	findFile = False
season = ARGS.season
testdir = ARGS.testdir
controldir = ARGS.controldir
clev = ARGS.clev
diffclev = ARGS.diffclev= False
savefig = ARGS.savefig
showfig = ARGS.showfig
variable = ARGS.variable
ftype = ARGS.file_type
mkdir = True
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savefig = False
	showfig = True
	mkdir = False

# Set the lat bounds
# Default global tropics
region_name = "GlobalTropics"
southern_lat = -85
northern_lat = 85
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
	if not os.path.exists("popDiffMap"):
		os.mkdir("popDiffMap")
		print "Created directory " + "popDiffMap"

	# Create the region directory if it doesn't already exist
	if not os.path.exists("popDiffMap/" + region_name):
		os.mkdir("popDiffMap/" + region_name)
		print "Created directory " + "popDiffMap/" + region_name

	# Create season directory inside region directory
	if not os.path.exists("popDiffMap/" + region_name + "/" + season):
		os.mkdir("popDiffMap/" + region_name + "/" + season)
		print "Created directory " + "popDiffMap/" + region_name + "/" + season

if findFile:
	# Look for the climo files in the root directory
	print "\nLooking for control " + season + " files in " + controldir + "..."
	controldatafname, controlfn = findClimoFile("*" + season + "*", controldir)
	if not controldatafname:
		sys.exit()
	else:
		print "Found file " + controlfn
	print "\nLooking for test " + season + " files in " + testdir + "..."
	testdatafname, testfn = findClimoFile("*" + season + "*", testdir)
	if not testdatafname:
		sys.exit()
	else:
		print "Found file " + testfn

else:
	print "\nControl file is " + controldatafname
	print "\nTest file is " + testdatafname
	controlfn = os.path.splitext(os.path.split(controldatafname)[1])[0]
	testfn = os.path.splitext(os.path.split(testdatafname)[1])[0]

# Read the data
controldata = popgoda(controldatafname)
testdata = popgoda(testdatafname)

# Set the boxlat and boxlon
box = (southern_lat, northern_lat, left_lon, right_lon)

g = -9.8 # gravitational constant

# Change into figure directory (root/region/season/) for image creation
if mkdir:
	os.chdir("popDiffMap/" + region_name + "/" + season)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------- Map creation ---------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for V in variable:
	var = V
	vname = V
	pressure = None

	print("\nPlotting " + season + " " + vname  +  " data...")

	try:
		testdata.variable(var, box)
		controldata.variable(var, box)
	except KeyError:
		print "Not able to plot variable " + var + "...\nSkipping this variable."
		continue
	
	controldata.variable("QFLUX", box, setData = False) # this sets self.boxlon and self.boxlat

	# Create bm coords from region bounds
	# bm lonitude coords need to be 0 < coord < 360 
	bmlon, bmlat = np.meshgrid(controldata.boxlon, controldata.boxlat)

	# Reset the lat and lon bounds so that maps don't show grey areas 
	southern_lat, northern_lat = np.array(controldata.boxlat)[[0,-1]]
	# Change lons to be negative is 180 < lon < 360 because that's how bm works for 'cea' projection
	left_lon, right_lon = np.array(controldata.boxlon)[[0,-1]]
	if 0 in controldata.boxlon[1:-2]: # if we cross the gml
		left_lon = controldata.boxlon[0]-360
		fig = plt.figure()

	#----------------#
	# Create the map #
	#----------------#
	fig = plt.figure()

	# test data
	testdatadata = testdata.data
	plt.subplot(3,1,1)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	# cs = m.contourf(bmlon, bmlat, testdata.data, testdata.clevs, shading = 'flat', latlon = True, cmap=testdata.cmap)
	cs = m.contourf(bmlon, bmlat, testdatadata, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label("", fontsize = 8)
	if user_test_plot_title is None:
		test_title = "test: " + testfn
	else:
		test_title = user_test_plot_title
	plt.title(test_title, fontsize = 8)
	
	# control data
	controldatadata = controldata.data
	plt.subplot(3,1,2)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	# cs = m.contourf(bmlon, bmlat, controldata.data, controldata.clevs, shading = 'flat', latlon = True, cmap=controldata.cmap)
	cs = m.contourf(bmlon, bmlat, controldatadata, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label("", fontsize = 8)
	control_title = user_control_plot_title 
	if user_control_plot_title is None:
		control_title = "control: " + controlfn
	plt.title(control_title, fontsize = 8)
	
	# difference data 
	plt.subplot(3,1,3)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	cs = m.contourf(bmlon, bmlat, testdatadata - controldatadata, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label("", fontsize = 8)
	plt.title("test - control difference", fontsize = 8)		
	
	# Make things pretty
	suptitle = testdata.long_name
	savetitle = vname
	savetitle += "_" + testfn + "-" + controlfn + "." + ftype
	fig.suptitle(suptitle, fontweight = 'bold', fontsize = 14)
	plt.subplots_adjust(hspace = .2)
	if savefig:
		plt.savefig(savetitle, bbox_inches='tight', dpi = 500)
		print("Created " + savetitle + '\n')
	if showfig:
		plt.show()
	plt.clf()
	plt.cla()
	plt.close(fig)