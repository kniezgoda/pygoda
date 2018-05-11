#!/glade/u/apps/ch/opt/python/2.7.13/gnu/6.2.0/bin/python
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
from pygoda import camgoda, findClimoFile, niceClev, RegularClev
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
parser.add_argument('-barbs', '--wind_barb_pressure', dest = 'wind_barb_pressure', nargs = 1, type = float, default = None)
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
diffclev = ARGS.diffclev
barbs = ARGS.wind_barb_pressure
show_barbs = False
if barbs is not None:
	barb_pressure = ARGS.wind_barb_pressure[0] * 100
	show_barbs = True
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
controldata = camgoda(controldatafname)
testdata = camgoda(testdatafname)

# Set the boxlat and boxlon
box = (southern_lat, northern_lat, left_lon, right_lon)
controldata.variable(controldata.vars[0], box, setData = False) # this sets self.boxlon and self.boxlat

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

for V in variable:
	# Extract data for wind barbs if needed
	if show_barbs:
		test_v = testdata.variable("V", box)
		test_v = testdata.isobar(barb_pressure)
		test_u = testdata.variable("U", box)
		test_u = testdata.isobar(barb_pressure)
		control_v = controldata.variable("V", box)
		control_v = controldata.isobar(barb_pressure)
		control_u = controldata.variable("U", box)
		control_u = controldata.isobar(barb_pressure)
		diff_v = test_v - control_v
		diff_u = test_u - control_u		
		
		density = 5
		colsby4 = np.arange(0, bmlon.shape[1], density)
		rowsby4 = np.arange(0, bmlat.shape[0], density)
		barb_lon = bmlon[:,colsby4]
		barb_lon = barb_lon[rowsby4,:]
		barb_lat = bmlat[:,colsby4]
		barb_lat = barb_lat[rowsby4,:]
			
		test_v = test_v[rowsby4,:]
		test_v = test_v[:,colsby4]
		test_u = test_u[rowsby4,:]
		test_u = test_u[:,colsby4]
		control_v = control_v[rowsby4,:]	
		control_v = control_v[:,colsby4]
		control_u = control_u[rowsby4,:]	
		control_u = control_u[:,colsby4]
		diff_v = diff_v[rowsby4,:]
		diff_v = diff_v[:,colsby4]
		diff_u = diff_u[rowsby4,:]
		diff_u = diff_u[:,colsby4]

	# Extract variable info (sets var, vname, and pressure)
	controldata.ExtractData(V, box)
	testdata.ExtractDate(V, box)
	
	fig = plt.figure()
	
	testdata.clevs = RegularClev(testdata.data) # This sets the clevs
	testdata.prep_map(season, region) # This sets the cmap 
	controldata.clevs = RegularClev(controldata.data) # This sets the clevs
	controldata.prep_map(season, region) # This sets the cmap 
	
	dclev = niceClev(testdata.data - controldata.data)
	
	if clev is not None:
		testdata.clevs = np.linspace(clev[0], clev[1], clev[2])
		controldata.clevs = np.linspace(clev[0], clev[1], clev[2])	
		
	if diffclev is not None:
		dclev = np.linspace(diffclev[0], diffclev[1], diffclev[2]) 

	# test data
	plt.subplot(3,1,1)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	cs = m.contourf(bmlon, bmlat, testdata.data, testdata.clevs, shading = 'flat', latlon = True, cmap=testdata.cmap)
	if show_barbs:
		m.quiver(barb_lon, barb_lat, test_u, test_v, scale = 50, scale_units = "inches", cmap = plt.cm.autumn, latlon = True)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label(testdata.units, fontsize = 8)
	if user_test_plot_title is None:
		test_title = "test: " + testfn
	else:
		test_title = user_test_plot_title
	plt.title(test_title, fontsize = 8)
	
	# control data
	plt.subplot(3,1,2)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	cs = m.contourf(bmlon, bmlat, controldata.data, controldata.clevs, shading = 'flat', latlon = True, cmap=controldata.cmap)
	if show_barbs:
		m.quiver(barb_lon, barb_lat, control_u, control_v, scale = 50, scale_units = "inches", cmap = plt.cm.autumn, latlon = True)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label(controldata.units, fontsize = 8)
	control_title = user_control_plot_title 
	if user_control_plot_title is None:
		control_title = "control: " + controlfn
	plt.title(control_title, fontsize = 8)
	
	# difference data 
	plt.subplot(3,1,3)
	m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
	m.drawcoastlines()
	m.drawmapboundary(fill_color='0.3')
	cs = m.contourf(bmlon, bmlat, testdata.data - controldata.data, dclev, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
	if show_barbs:
		m.quiver(barb_lon, barb_lat, diff_u, diff_v, scale = 8, scale_units = "inches", cmap = plt.cm.autumn, latlon = True)
	cbar = m.colorbar(cs, location='right', pad="5%")
	cbar.set_label("difference in " + testdata.units, fontsize = 8)
	plt.title("test - control difference", fontsize = 8)		
	
	# Make things pretty
	suptitle = testdata.long_name
	savetitle = vname
	if var_is_3d:
		suptitle += " at " + str(pressure/100) + " mb"
	if show_barbs:
		suptitle += " with " + str(barb_pressure/100) + " hPa winds"
		savetitle += str(int(round(barb_pressure/100))) + "winds"
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
