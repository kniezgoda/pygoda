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

import os, glob
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap as bm 
from pygoda import ncgoda

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

######################################
####### User-entered data here #######
######################################

#------------- Set your control and test directory here:
# These directories should contain the ANN, DJF, and JJA climo files computed from the AMWG diagnostic package

#root = "/Users/kyleniezgoda/LocalDocuments/Yellowstone_Runs/model_runs/"
#controldir = root + "F.C5.2deg.wiso.defaultSSTICE_kn002/"
#testdir = root + "F.C5.2deg.wiso.obs6kSST_kn003/"

# COAS Climate machine
root = "/home/server/student/homes/kniezgod/model_runs" 
controldir = root +  "/F.C5.2deg.wiso.defaultSSTICE_kn002/"
testdir = root + "/F.C5.2deg.wiso.obs6kSST_kn003/"
#------------- 

#------------- If saving figure, what directory should the file be saved to?
# Three directories will be created inside this folder, called ANN, JJA, and DJF, if they don't already exist

# Kyle's work computer
#figdir = "/Users/kyleniezgoda/Library/Mobile Documents/com~apple~CloudDocs/Documents/Python/figures/F.C5.2deg.wiso.obs6kSST_kn003-F.C5.2deg.wiso.defaultSSTICE_kn002/"

# COAS Climate machine
figdir = "/home/server/student/homes/kniezgod/python/figures/F.C5.2deg.wiso.obs6kSST_kn003-F.C5.2deg.wiso.defaultSSTICE_kn002/"
#------------- 

#------------- Set region (Options: EP, IM, MC, NA, do not set to None. Set to "" for default global tropics)
region = "EP"

#------------- Set seasons 
# seasons = ["ANN", "DJF", "JJA"]
seasons = ["JJA"]

#------------- Which variables to plot? (1 -> True, 0 -> False)
PRECT=1
PRECC=1
PRECL=1
PRECSC=1
PRECSL=1
PRECST=1
QFLX=1
CLDHGH=1
PS=1
PSL=1
PRECT_d18O=1
PRECT_dD=1
PRECT_dxs=1
Q850=1
Q500=1
Q200=1
V850=1
V500=1
V200=1
VT850=1
VT500=1
VT200=1
VQ850=1
VQ500=1
VQ200=1
U850=1
U500=1
U200=1
UT850=1
UT500=1
UT200=1
UQ850=1
UQ500=1
UQ200=1
T850=1
T500=1
T200=1
OMEGA850=1
OMEGA500=1
OMEGA200=1
Z3850=1
Z3500=1
Z3200=1
dDV850=1
dDV500=1
dDV200=1
d18OV850=1
d18OV500=1
d18OV200=1
dxsV850=1
dxsV500=1
dxsV200=1


#------------- Show figures?
# This will take a while and it will show literally every plot along the way.
show_figure = False

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
# Create variable string list to loop through
var2dlist = []
if PRECT:
	var2dlist.append("PRECT")
if PRECC:
	var2dlist.append("PRECC")
if PRECL:
	var2dlist.append("PRECL")
if PRECSC:
	var2dlist.append("PRECSC")
if PRECSL:
	var2dlist.append("PRECSL")
if PRECST:
	var2dlist.append("PRECST")
if QFLX:
	var2dlist.append("QFLX")
if CLDHGH:
	var2dlist.append("CLDHGH")
if PS:
	var2dlist.append("PS")
if PSL:
	var2dlist.append("PSL")
if PRECT_d18O:
	var2dlist.append("PRECT_d18O")
if PRECT_dD:
	var2dlist.append("PRECT_dD")
if PRECT_dxs:
	var2dlist.append("PRECT_dxs")

var3dlist = []
if Q850:
	var3dlist.append("Q850")
if Q500:
	var3dlist.append("Q500")
if Q200:
	var3dlist.append("Q200")
if V850:
	var3dlist.append("V850")
if V500:
	var3dlist.append("V500")
if V200:
	var3dlist.append("V200")
if VT850:
	var3dlist.append("VT850")
if VT500:
	var3dlist.append("VT500")
if VT200:
	var3dlist.append("VT200")
if VQ850:
	var3dlist.append("VQ850")
if VQ500:
	var3dlist.append("VQ500")
if VQ200:
	var3dlist.append("VQ200")
if U850:
	var3dlist.append("U850")
if U500:
	var3dlist.append("U500")
if U200:
	var3dlist.append("U200")
if UT850:
	var3dlist.append("UT850")
if UT500:
	var3dlist.append("UT500")
if UT200:
	var3dlist.append("UT200")
if UQ850:
	var3dlist.append("UQ850")
if UQ500:
	var3dlist.append("UQ500")
if UQ200:
	var3dlist.append("UQ200")
if T850:
	var3dlist.append("T850")
if T500:
	var3dlist.append("T500")
if T200:
	var3dlist.append("T200")
if OMEGA850:
	var3dlist.append("OMEGA850")
if OMEGA500:
	var3dlist.append("OMEGA500")
if OMEGA200:
	var3dlist.append("OMEGA200")
if Z3850:
	var3dlist.append("Z3850")
if Z3500:
	var3dlist.append("Z3500")
if Z3200:
	var3dlist.append("Z3200")
if dDV850:
	var3dlist.append("dDV850")
if dDV500:
	var3dlist.append("dDV500")
if dDV200:
	var3dlist.append("dDV200")
if d18OV850:
	var3dlist.append("d18OV850")
if d18OV500:
	var3dlist.append("d18OV500")
if d18OV200:
	var3dlist.append("d18OV200")
if dxsV850:
	var3dlist.append("dxsV850")
if dxsV500:
	var3dlist.append("dxsV500")
if dxsV200:
	var3dlist.append("dxsV200")

# Set the lat bounds
# Default global tropics
region_name = "GlobalTropics"
southern_lat = -50
northern_lat = 50
left_lon = 0
right_lon = 355

# Indian monsoon
if region == "IM":
	region_name = "IndianMonsoon"
	southern_lat = -5
	northern_lat = 45
	left_lon = 40
	right_lon = 110

# Greater maritime continent
if region == "MC":
	region_name = "MaritimeContinent"
	southern_lat = -20
	northern_lat = 20
	left_lon = 80
	right_lon = 160

# North Africa
if region == "NA":
	region_name = "NorthAfrica"
	southern_lat = -20
	northern_lat = 45
	left_lon = -30
	right_lon = 70

#Central Equatorial Pacific to Western eq. Atlantic
if region == "EP":
	region_name = "TropicalOceans"
	southern_lat = -25
	northern_lat = 25
	left_lon = 180
	right_lon = 355

figdir += region_name

# Create the figdir if it doesn't already exist
if not os.path.exists(figdir):
	os.mkdir(figdir)
	print "Created directory " + figdir

os.chdir(figdir)

#------------- Beginning of code -------------#
for N, which in enumerate(seasons):
	# Create the directory for the season inside the figdir
	if not os.path.exists(which):
		os.mkdir(which)
		print "Created directory " + figdir + "/" + which

	# Look for the climo files
	print "Looking for control files in " + controldir
	controldatafname = glob.glob(controldir+"/*"+which+"*.nc")
	print "Looking for test files in " + testdir
	testdatafname = glob.glob(testdir+"/*"+which+"*.nc")

	if len(controldatafname) > 0:
		print "\nFound control " + which + " file: " + os.path.split(controldatafname[0])[1]
		controldatafname = controldatafname[0]

	else:
		print "Could not find control " + which + " file..."
		print "No " + which + " plots will be created..."
		continue # skip this iteration

	if len(testdatafname) > 0:
		print "Found test " + which + " file: " + os.path.split(testdatafname[0])[1]
		testdatafname = testdatafname[0]

	else:
		print "Could not find test " + which + " file..."
		print "No " + which + " plots will be created..."
		continue # skip this iteration

	# Assign file names for plotting later
	controlfn = os.path.splitext(os.path.split(controldatafname)[1])[0]
	testfn = os.path.splitext(os.path.split(testdatafname)[1])[0]

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
	left_lon, right_lon = np.array([l-360 if 180<l<360 else l for l in controldata.boxlon])[[0,-1]]
	
	g = -9.8 # gravitational constant

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------- Map creation ---------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	for variable in var2dlist:

		print("\nPlotting " + which + " " + variable + " data...")
		fig = plt.figure()

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

		testdata.set_clevs(which, region)
		controldata.set_clevs(which, region)

		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, testdata.data, testdata.clevs, shading = 'flat', latlon = True, cmap=testdata.cmap)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(testdata.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
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
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + controlfn
		plt.title(control_title, fontsize = 8)
		
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, testdata.data - controldata.data, testdata.diffclevs, shading = 'flat', latlon = True, cmap=testdata.diffcmap)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + testdata.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		fig.suptitle(testdata.long_name, fontweight = 'bold', fontsize = 14)
		plt.subplots_adjust(hspace = .2)
		plt.savefig(figdir + "/"  + which + "/" + variable + "_" + testfn + "-" + controlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + variable + "_" + testfn + "-" + controlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)

	#---------------------------------#

	for variable in var3dlist:

		print("\nPlotting " + which + " " + variable + " data...")
		fig = plt.figure()

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

		testdata.set_clevs(which, region)
		controldata.set_clevs(which, region)

		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, testdata.data, 19, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(testdata.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, controldata.data, 19, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(controldata.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + controlfn
		plt.title(control_title, fontsize = 8)
		
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, testdata.data - controldata.data, 19, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + testdata.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		fig.suptitle(testdata.long_name + " at " + str(pressure/100) + " mb", fontweight = 'bold', fontsize = 14)
		plt.subplots_adjust(hspace = .2)
		plt.savefig(figdir + "/"  + which + "/" + variable + "_" + testfn + "-" + controlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + variable + "_" + testfn + "-" + controlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)


# The following writes the directory that the figures are stored in to a hidden text file on the user's home directory
# This can be used in, for example, a shell script, to read where the figures are located without having to hard code it in
# E.g.: image_dir='cat ~/.tempfile.txt' in bash will store the image directory in the env variable called $image_dir
home = os.path.expanduser("~")
figdirtemp = open(home+"/.figdirtemp.txt", "w")
figdirtemp.write(figdir)

regiontemp = open(home+"/.regiontemp.txt", "w")
regiontemp.write(region_name)


