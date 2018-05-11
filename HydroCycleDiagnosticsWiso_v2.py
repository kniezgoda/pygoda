#!/opt/anaconda2/bin/python
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
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import shiftgrid 
from pygoda import ncgoda
from HydroCycleDiag_classes import *

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

######################################
####### User-entered data here #######
######################################

#------------- Set your control and test directory here:
# These directories should contain the ANN, DJF, and JJA climo files computed from the AMWG diagnostic package

# Kyle's work computer
# root = '/Users/kyleniezgoda/LocalDocuments/Yellowstone_Runs/model_runs/'
# cntrldir = root+"F.C5.2deg.wiso.defaultSSTICE_kn002/"
# testdir = root+"F.C5.2deg.wiso.obs6kSST_kn003/"

# COAS Climate machine
root = '/home/server/student/homes/kniezgod/model_runs/' 
cntrldir = root+"F.C5.2deg.wiso.defaultSSTICE_kn002/"
testdir = root+"F.C5.2deg.wiso.obs6kSST_kn003/"
#------------- 

#------------- If saving figure, what directory should the file be saved to?
# Three directories will be created inside this folder, called ANN, JJA, and DJF, if they don't already exist

# Kyle's work computer
# figdir = "/Users/kyleniezgoda/Library/Mobile Documents/com~apple~CloudDocs/Documents/Python/figures/F.C5.2deg.wiso.obs6kSST_kn003-F.C5.2deg.wiso.defaultSSTICE_kn002/"

# COAS Climate machine
figdir = "/home/server/student/homes/kniezgod/python/figures/F.C5.2deg.wiso.obs6kSST_kn003-F.C5.2deg.wiso.defaultSSTICE_kn002/"
#------------- 

#------------- Set region (Options: EP, IM, MC, NA, do not set to None. Set to "" for default global tropics)
region = ""

#------------- Set seasons 
# seasons = ["ANN", "DJF", "JJA"]
seasons = ["JJA"]

#------------- Which variables to plot? (1 -> True, 0 -> False)
plot_PRECT_H2O=0
plot_PRECT_H218O=0
plot_PRECT_HDO=0
plot_PRECC=0
plot_PRECL=0
plot_PRECSC=0
plot_PRECSL=0
plot_PRECST=0
plot_QFLX=0
plot_QFLX_H218O=0
plot_QFLX_HDO=0
plot_CLDHIGH=0
plot_PS=0
plot_PSL=0
plot_H2OV=0
plot_HDOV=0
plot_H218OV=0
plot_Qsurface=0
plot_Q850=0
plot_Q500=0
plot_Q200=0
plot_Vsurface=0
plot_V850=0
plot_V500=0
plot_V200=0
plot_VTsurface=0
plot_VT850=0
plot_VT500=0
plot_VT200=0
plot_VQsurface=0
plot_VQ850=0
plot_VQ500=0
plot_VQ200=0
plot_Usurface=0
plot_U850=0
plot_U500=0
plot_U200=0
plot_UTsurface=0
plot_UT850=0
plot_UT500=0
plot_UT200=0
plot_UQsurface=0
plot_UQ850=0
plot_UQ500=0
plot_UQ200=0
plot_Tsurface=0
plot_T850=0
plot_T500=0
plot_T200=0
plot_OMEGAsurface=0
plot_OMEGA850=0
plot_OMEGA500=0
plot_OMEGA200=0
plot_Z3surface=0
plot_Z3850=0
plot_Z3500=0
plot_Z3200=0
plot_dDVsurface=0
plot_dDV850=0
plot_dDV500=0
plot_dDV200=0
plot_d18OVsurface=0
plot_d18OV850=0
plot_d18OV500=0
plot_d18OV200=0
plot_dxsVsurface=0
plot_dxsV850=0
plot_dxsV500=0
plot_dxsV200=0
plot_PRECT_d18O=1
plot_PRECT_dD=1
plot_PRECT_dxs=1





#------------- Show figures?
# This will take a while and it will show literally every plot along the way.
show_figure = True

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
	print "Looking for control files in " + cntrldir
	cntrldatafname = glob.glob(cntrldir+"*"+which+"*.nc")
	print "Looking for test files in " + testdir
	testdatafname = glob.glob(testdir+"*"+which+"*.nc")

	if len(cntrldatafname) > 0:
		print "\nFound control " + which + " file: " + os.path.split(cntrldatafname[0])[1]
		cntrldatafname = cntrldatafname[0]

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
	cntrlfn = os.path.splitext(os.path.split(cntrldatafname)[1])[0]
	testfn = os.path.splitext(os.path.split(testdatafname)[1])[0]

	# Read the data
	cntrldata = Dataset(cntrldatafname, mode='r')
	testdata = Dataset(testdatafname, mode='r')

	# Extract lat, lon, and lev data
	lats = cntrldata.variables['lat'][:]
	lons = cntrldata.variables['lon'][:]
	box = find_indices([southern_lat, northern_lat, left_lon, right_lon], lats, lons)
	bmlon, bmlat = np.meshgrid(lons, lats)
	levs = cntrldata.variables['lev'][:]
	A_control = cntrldata.variables['hyam'][:]
	A_test = testdata.variables['hyam'][:]
	B_control = cntrldata.variables['hybm'][:]
	B_test = testdata.variables['hybm'][:]

	g = -9.8 # gravitational constant
	loaded_vars = []

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------- Map creation ---------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	#------------- PRECT_H2O -------------#
	if plot_PRECT_H2O:
		if 'prect_h2o' not in loaded_vars:
			prect_h2o = PRECT_H2O(cntrldata, testdata)
			loaded_vars.append('prect_h2o')
		print("\nPlotting " + which + " PRECT_H2O data...")
		fig = plt.figure()
		fig.suptitle(prect_h2o.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)x
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_h2o.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, prect_h2o.test, prect_h2o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_h2o.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, prect_h2o.control, prect_h2o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_h2o.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_h2o.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, prect_h2o.diff, prect_h2o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + prect_h2o.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT_H2O_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT_H2O_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECT_H218O -------------#
	if plot_PRECT_H218O:
		if 'prect_h218o' not in loaded_vars:
			prect_h218o = PRECT_H218O(cntrldata, testdata)
			loaded_vars.append('prect_h218o')
		print("\nPlotting " + which + " PRECT_H218O data...")
		fig = plt.figure()
		fig.suptitle(prect_h218o.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_h218o.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, prect_h218o.test, prect_h218o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_h218o.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, prect_h218o.control, prect_h218o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_h218o.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_h218o.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, prect_h218o.diff, prect_h218o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + prect_h218o.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT_H218O_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT_H218O_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECT_HDO -------------#
	if plot_PRECT_HDO:
		if 'prect_hdo' not in loaded_vars:
			prect_hdo = PRECT_HDO(cntrldata, testdata)
			loaded_vars.append('prect_hdo')
		print("\nPlotting " + which + " PRECT_HDO data...")
		fig = plt.figure()
		fig.suptitle(prect_hdo.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_hdo.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, prect_hdo.test, prect_hdo.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_hdo.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, prect_hdo.control, prect_hdo.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_hdo.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_hdo.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, prect_hdo.diff, prect_hdo.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + prect_hdo.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT_HDO_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT_HDO_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECC -------------#
	if plot_PRECC:
		if 'precc' not in loaded_vars:
			precc = PRECC(cntrldata, testdata)
			loaded_vars.append('precc')
		print("\nPlotting " + which + " PRECC data...")
		fig = plt.figure()
		fig.suptitle(precc.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precc.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, precc.test, precc.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precc.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, precc.control, precc.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precc.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precc.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, precc.diff, precc.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + precc.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECC_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECC_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECL -------------#
	if plot_PRECL:
		if 'precl' not in loaded_vars:
			precl = PRECL(cntrldata, testdata)
			loaded_vars.append('precl')
		print("\nPlotting " + which + " PRECL data...")
		fig = plt.figure()
		fig.suptitle(precl.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precl.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, precl.test, precl.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precl.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, precl.control, precl.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precl.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precl.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, precl.diff, precl.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + precl.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECL_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECL_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECSC -------------#
	if plot_PRECSC:
		if 'precsc' not in loaded_vars:
			precsc = PRECSC(cntrldata, testdata)
			loaded_vars.append('precsc')
		print("\nPlotting " + which + " PRECSC data...")
		fig = plt.figure()
		fig.suptitle(precsc.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precsc.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, precsc.test, precsc.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precsc.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, precsc.control, precsc.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precsc.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precsc.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, precsc.diff, precsc.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + precsc.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECSC_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECSC_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECSL -------------#
	if plot_PRECSL:
		if 'precsl' not in loaded_vars:
			precsl = PRECSL(cntrldata, testdata)
			loaded_vars.append('precsl')
		print("\nPlotting " + which + " PRECSL data...")
		fig = plt.figure()
		fig.suptitle(precsl.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precsl.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, precsl.test, precsl.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precsl.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, precsl.control, precsl.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precsl.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precsl.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, precsl.diff, precsl.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + precsl.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECSL_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECSL_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECST -------------#
	if plot_PRECST:
		if 'precst' not in loaded_vars:
			precst = PRECST(cntrldata, testdata)
			loaded_vars.append('precst')
		print("\nPlotting " + which + " PRECST data...")
		fig = plt.figure()
		fig.suptitle(precst.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precst.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, precst.test, precst.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precst.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, precst.control, precst.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(precst.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		precst.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, precst.diff, precst.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + precst.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECST_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECST_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- QFLX -------------#
	if plot_QFLX:
		if 'qflx' not in loaded_vars:
			qflx = QFLX(cntrldata, testdata)
			loaded_vars.append('qflx')
		print("\nPlotting " + which + " QFLX data...")
		fig = plt.figure()
		fig.suptitle(qflx.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qflx.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, qflx.test, qflx.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(qflx.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, qflx.control, qflx.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(qflx.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qflx.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, qflx.diff, qflx.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + qflx.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "QFLX_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "QFLX_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- QFLX_H218O -------------#
	if plot_QFLX_H218O:
		if 'qflx_h218o' not in loaded_vars:
			qflx_h218o = QFLX_H218O(cntrldata, testdata)
			loaded_vars.append('qflx_h218o')
		print("\nPlotting " + which + " QFLX_H218O data...")
		fig = plt.figure()
		fig.suptitle(qflx_h218o.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qflx_h218o.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, qflx_h218o.test, qflx_h218o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(qflx_h218o.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, qflx_h218o.control, qflx_h218o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(qflx_h218o.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qflx_h218o.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, qflx_h218o.diff, qflx_h218o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + qflx_h218o.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "QFLX_H218O_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "QFLX_H218O_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- QFLX_HDO -------------#
	if plot_QFLX_HDO:
		if 'qflx_hdo' not in loaded_vars:
			qflx_hdo = QFLX_HDO(cntrldata, testdata)
			loaded_vars.append('qflx_hdo')
		print("\nPlotting " + which + " QFLX_HDO data...")
		fig = plt.figure()
		fig.suptitle(qflx_hdo.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qflx_hdo.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, qflx_hdo.test, qflx_hdo.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(qflx_hdo.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, qflx_hdo.control, qflx_hdo.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(qflx_hdo.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qflx_hdo.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, qflx_hdo.diff, qflx_hdo.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + qflx_hdo.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "QFLX_HDO_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "QFLX_HDO_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- CLDHIGH -------------#
	if plot_CLDHIGH:
		if 'cldhigh' not in loaded_vars:
			cldhigh = CLDHIGH(cntrldata, testdata)
			loaded_vars.append('cldhigh')
		print("\nPlotting " + which + " CLDHGH data...")
		fig = plt.figure()
		fig.suptitle(cldhigh.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cldhigh.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, cldhigh.test, cldhigh.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(cldhigh.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, cldhigh.control, cldhigh.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(cldhigh.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cldhigh.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, cldhigh.diff, cldhigh.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + cldhigh.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "CLDHIGH_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "CLDHIGH_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PS -------------#
	if plot_PS:
		if 'ps' not in loaded_vars:
			ps = PS(cntrldata, testdata)
			loaded_vars.append('ps')
		print("\nPlotting " + which + " PS data...")
		fig = plt.figure()
		fig.suptitle(ps.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ps.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, ps.test, ps.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ps.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, ps.control, ps.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ps.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ps.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, ps.diff, ps.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + ps.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PS_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PS_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PSL -------------#
	if plot_PSL:
		if 'psl' not in loaded_vars:
			psl = PSL(cntrldata, testdata)
			loaded_vars.append('psl')
		print("\nPlotting " + which + " PSL data...")
		fig = plt.figure()
		fig.suptitle(psl.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		psl.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, psl.test, psl.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(psl.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, psl.control, psl.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(psl.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		psl.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, psl.diff, psl.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + psl.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PSL_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PSL_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- H2OV -------------#
	if plot_H2OV:
		if 'h2ov' not in loaded_vars:
			h2ov = H2OV(cntrldata, testdata)
			loaded_vars.append('h2ov')
		print("\nPlotting " + which + " H2OV data...")
		fig = plt.figure()
		fig.suptitle(h2ov.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		h2ov.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, h2ov.test, h2ov.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(h2ov.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, h2ov.control, h2ov.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(h2ov.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		h2ov.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, h2ov.diff, h2ov.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + h2ov.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "H2OV_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "H2OV_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- HDOV -------------#
	if plot_HDOV:
		if 'hdov' not in loaded_vars:
			hdov = HDOV(cntrldata, testdata)
			loaded_vars.append('hdov')
		print("\nPlotting " + which + " HDOV data...")
		fig = plt.figure()
		fig.suptitle(hdov.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		hdov.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, hdov.test, hdov.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(hdov.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, hdov.control, hdov.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(hdov.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		hdov.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, hdov.diff, hdov.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + hdov.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "HDOV_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "HDOV_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- H218OV -------------#
	if plot_H218OV:
		if 'h218ov' not in loaded_vars:
			h218ov = H218OV(cntrldata, testdata)
			loaded_vars.append('h218ov')
		print("\nPlotting " + which + " H218OV data...")
		fig = plt.figure()
		fig.suptitle(h218ov.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		h218ov.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, h218ov.test, h218ov.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(h218ov.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, h218ov.control, h218ov.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(h218ov.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		h218ov.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, h218ov.diff, h218ov.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + h218ov.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "H218OV_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "H218OV_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Qsurface -------------#
	if plot_Qsurface:
		if 'qsurface' not in loaded_vars:
			qsurface = Qsurface(cntrldata, testdata)
			loaded_vars.append('qsurface')
		print("\nPlotting " + which + " Qsurface data...")
		fig = plt.figure()
		fig.suptitle(qsurface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, qsurface.test, qsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(qsurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, qsurface.control, qsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(qsurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, qsurface.diff, qsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + qsurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Qsurface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Qsurface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Q850 -------------#
	if plot_Q850:
		if 'q850' not in loaded_vars:
			q850 = Q850(cntrldata, testdata)
			loaded_vars.append('q850')
		print("\nPlotting " + which + " Q850 data...")
		fig = plt.figure()
		fig.suptitle(q850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		q850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, q850.test, q850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(q850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, q850.control, q850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(q850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		q850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, q850.diff, q850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + q850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Q850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Q850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Q500 -------------#
	if plot_Q500:
		if 'q500' not in loaded_vars:
			q500 = Q500(cntrldata, testdata)
			loaded_vars.append('q500')
		print("\nPlotting " + which + " Q500 data...")
		fig = plt.figure()
		fig.suptitle(q500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		q500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, q500.test, q500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(q500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, q500.control, q500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(q500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		q500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, q500.diff, q500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + q500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Q500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Q500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Q200 -------------#
	if plot_Q200:
		if 'q200' not in loaded_vars:
			q200 = Q200(cntrldata, testdata)
			loaded_vars.append('q200')
		print("\nPlotting " + which + " Q200 data...")
		fig = plt.figure()
		fig.suptitle(q200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		q200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, q200.test, q200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(q200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, q200.control, q200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(q200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		q200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, q200.diff, q200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + q200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Q200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Q200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Vsurface -------------#
	if plot_Vsurface:
		if 'vsurface' not in loaded_vars:
			vsurface = Vsurface(cntrldata, testdata)
			loaded_vars.append('vsurface')
		print("\nPlotting " + which + " Vsurface data...")
		fig = plt.figure()
		fig.suptitle(vsurface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, vsurface.test, vsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vsurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, vsurface.control, vsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vsurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, vsurface.diff, vsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + vsurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Vsurface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Vsurface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- V850 -------------#
	if plot_V850:
		if 'v850' not in loaded_vars:
			v850 = V850(cntrldata, testdata)
			loaded_vars.append('v850')
		print("\nPlotting " + which + " V850 data...")
		fig = plt.figure()
		fig.suptitle(v850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		v850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, v850.test, v850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(v850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, v850.control, v850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(v850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		v850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, v850.diff, v850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + v850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "V850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "V850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- V500 -------------#
	if plot_V500:
		if 'v500' not in loaded_vars:
			v500 = V500(cntrldata, testdata)
			loaded_vars.append('v500')
		print("\nPlotting " + which + " V500 data...")
		fig = plt.figure()
		fig.suptitle(v500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		v500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, v500.test, v500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(v500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, v500.control, v500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(v500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		v500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, v500.diff, v500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + v500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "V500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "V500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- V200 -------------#
	if plot_V200:
		if 'v200' not in loaded_vars:
			v200 = V200(cntrldata, testdata)
			loaded_vars.append('v200')
		print("\nPlotting " + which + " V200 data...")
		fig = plt.figure()
		fig.suptitle(v200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		v200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, v200.test, v200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(v200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, v200.control, v200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(v200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		v200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, v200.diff, v200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + v200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "V200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "V200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- VTsurface -------------#
	if plot_VTsurface:
		if 'vtsurface' not in loaded_vars:
			vtsurface = VTsurface(cntrldata, testdata)
			loaded_vars.append('vtsurface')
		print("\nPlotting " + which + " VTsurface data...")
		fig = plt.figure()
		fig.suptitle(vtsurface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vtsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, vtsurface.test, vtsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vtsurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, vtsurface.control, vtsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vtsurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vtsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, vtsurface.diff, vtsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + vtsurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VTsurface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VTsurface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- VT850 -------------#
	if plot_VT850:
		if 'vt850' not in loaded_vars:
			vt850 = VT850(cntrldata, testdata)
			loaded_vars.append('vt850')
		print("\nPlotting " + which + " VT850 data...")
		fig = plt.figure()
		fig.suptitle(vt850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vt850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, vt850.test, vt850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vt850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, vt850.control, vt850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vt850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vt850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, vt850.diff, vt850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + vt850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VT850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VT850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- VT500 -------------#
	if plot_VT500:
		if 'vt500' not in loaded_vars:
			vt500 = VT500(cntrldata, testdata)
			loaded_vars.append('vt500')
		print("\nPlotting " + which + " VT500 data...")
		fig = plt.figure()
		fig.suptitle(vt500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vt500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, vt500.test, vt500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vt500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, vt500.control, vt500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vt500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vt500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, vt500.diff, vt500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + vt500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VT500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VT500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- VT200 -------------#
	if plot_VT200:
		if 'vt200' not in loaded_vars:
			vt200 = VT200(cntrldata, testdata)
			loaded_vars.append('vt200')
		print("\nPlotting " + which + " VT200 data...")
		fig = plt.figure()
		fig.suptitle(vt200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vt200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, vt200.test, vt200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vt200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, vt200.control, vt200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vt200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vt200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, vt200.diff, vt200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + vt200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VT200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VT200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- VQsurface -------------#
	if plot_VQsurface:
		if 'vqsurface' not in loaded_vars:
			vqsurface = VQsurface(cntrldata, testdata)
			loaded_vars.append('vqsurface')
		print("\nPlotting " + which + " VQsurface data...")
		fig = plt.figure()
		fig.suptitle(vqsurface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vqsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, vqsurface.test, vqsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vqsurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, vqsurface.control, vqsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vqsurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vqsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, vqsurface.diff, vqsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + vqsurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQsurface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQsurface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- VQ850 -------------#
	if plot_VQ850:
		if 'vq850' not in loaded_vars:
			vq850 = VQ850(cntrldata, testdata)
			loaded_vars.append('vq850')
		print("\nPlotting " + which + " VQ850 data...")
		fig = plt.figure()
		fig.suptitle(vq850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vq850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, vq850.test, vq850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vq850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, vq850.control, vq850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vq850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vq850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, vq850.diff, vq850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + vq850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQ850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQ850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- VQ500 -------------#
	if plot_VQ500:
		if 'vq500' not in loaded_vars:
			vq500 = VQ500(cntrldata, testdata)
			loaded_vars.append('vq500')
		print("\nPlotting " + which + " VQ500 data...")
		fig = plt.figure()
		fig.suptitle(vq500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vq500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, vq500.test, vq500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vq500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, vq500.control, vq500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vq500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vq500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, vq500.diff, vq500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + vq500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQ500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQ500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- VQ200 -------------#
	if plot_VQ200:
		if 'vq200' not in loaded_vars:
			vq200 = VQ200(cntrldata, testdata)
			loaded_vars.append('vq200')
		print("\nPlotting " + which + " VQ200 data...")
		fig = plt.figure()
		fig.suptitle(vq200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vq200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, vq200.test, vq200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vq200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, vq200.control, vq200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(vq200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		vq200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, vq200.diff, vq200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + vq200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQ200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQ200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Usurface -------------#
	if plot_Usurface:
		if 'qsurface' not in loaded_vars:
			usurface = Usurface(cntrldata, testdata)
			loaded_vars.append('usurface')
		print("\nPlotting " + which + " Usurface data...")
		fig = plt.figure()
		fig.suptitle(usurface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, usurface.test, usurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(usurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, usurface.control, usurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(usurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		qsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, usurface.diff, usurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + usurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Usurface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Usurface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- U850 -------------#
	if plot_U850:
		if 'u850' not in loaded_vars:
			u850 = U850(cntrldata, testdata)
			loaded_vars.append('u850')
		print("\nPlotting " + which + " U850 data...")
		fig = plt.figure()
		fig.suptitle(u850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		u850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, u850.test, u850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(u850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, u850.control, u850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(u850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		u850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, u850.diff, u850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + u850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "U850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "U850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- U500 -------------#
	if plot_U500:
		if 'u500' not in loaded_vars:
			u500 = U500(cntrldata, testdata)
			loaded_vars.append('u500')
		print("\nPlotting " + which + " U500 data...")
		fig = plt.figure()
		fig.suptitle(u500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		u500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, u500.test, u500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(u500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, u500.control, u500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(u500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		u500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, u500.diff, u500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + u500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "U500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "U500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- U200 -------------#
	if plot_U200:
		if 'u200' not in loaded_vars:
			u200 = U200(cntrldata, testdata)
			loaded_vars.append('u200')
		print("\nPlotting " + which + " U200 data...")
		fig = plt.figure()
		fig.suptitle(u200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		u200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, u200.test, u200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(u200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, u200.control, u200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(u200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		u200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, u200.diff, u200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + u200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "U200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "U200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- UTsurface -------------#
	if plot_UTsurface:
		if 'utsurface' not in loaded_vars:
			utsurface = UTsurface(cntrldata, testdata)
			loaded_vars.append('utsurface')
		print("\nPlotting " + which + " UTsurface data...")
		fig = plt.figure()
		fig.suptitle(utsurface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		utsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, utsurface.test, utsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(utsurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, utsurface.control, utsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(utsurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		utsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, utsurface.diff, utsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + utsurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UTsurface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UTsurface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- UT850 -------------#
	if plot_UT850:
		if 'ut850' not in loaded_vars:
			ut850 = UT850(cntrldata, testdata)
			loaded_vars.append('ut850')
		print("\nPlotting " + which + " UT850 data...")
		fig = plt.figure()
		fig.suptitle(ut850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ut850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, ut850.test, ut850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ut850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, ut850.control, ut850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ut850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ut850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, ut850.diff, ut850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + ut850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UT850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UT850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- UT500 -------------#
	if plot_UT500:
		if 'ut500' not in loaded_vars:
			ut500 = UT500(cntrldata, testdata)
			loaded_vars.append('ut500')
		print("\nPlotting " + which + " UT500 data...")
		fig = plt.figure()
		fig.suptitle(ut500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ut500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, ut500.test, ut500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ut500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, ut500.control, ut500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ut500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ut500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, ut500.diff, ut500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + ut500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UT500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UT500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- UT200 -------------#
	if plot_UT200:
		if 'ut200' not in loaded_vars:
			ut200 = UT200(cntrldata, testdata)
			loaded_vars.append('ut200')
		print("\nPlotting " + which + " UT200 data...")
		fig = plt.figure()
		fig.suptitle(ut200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ut200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, ut200.test, ut200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ut200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, ut200.control, ut200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ut200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ut200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, ut200.diff, ut200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + ut200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UT200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UT200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- UQsurface -------------#
	if plot_UQsurface:
		if 'uqsurface' not in loaded_vars:
			uqsurface = UQsurface(cntrldata, testdata)
			loaded_vars.append('uqsurface')
		print("\nPlotting " + which + " UQsurface data...")
		fig = plt.figure()
		fig.suptitle(uqsurface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		uqsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, uqsurface.test, uqsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(uqsurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, uqsurface.control, uqsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(uqsurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		uqsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, uqsurface.diff, uqsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + uqsurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQsurface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQsurface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- UQ850 -------------#
	if plot_UQ850:
		if 'uq850' not in loaded_vars:
			uq850 = UQ850(cntrldata, testdata)
			loaded_vars.append('uq850')
		print("\nPlotting " + which + " UQ850 data...")
		fig = plt.figure()
		fig.suptitle(uq850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		uq850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, uq850.test, uq850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(uq850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, uq850.control, uq850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(uq850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		uq850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, uq850.diff, uq850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + uq850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQ850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQ850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- UQ500 -------------#
	if plot_UQ500:
		if 'uq500' not in loaded_vars:
			uq500 = UQ500(cntrldata, testdata)
			loaded_vars.append('uq500')
		print("\nPlotting " + which + " UQ500 data...")
		fig = plt.figure()
		fig.suptitle(uq500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		uq500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, uq500.test, uq500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(uq500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, uq500.control, uq500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(uq500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		uq500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, uq500.diff, uq500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + uq500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQ500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQ500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- UQ200 -------------#
	if plot_UQ200:
		if 'uq200' not in loaded_vars:
			uq200 = UQ200(cntrldata, testdata)
			loaded_vars.append('uq200')
		print("\nPlotting " + which + " UQ200 data...")
		fig = plt.figure()
		fig.suptitle(uq200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		uq200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, uq200.test, uq200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(uq200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, uq200.control, uq200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(uq200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		uq200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, uq200.diff, uq200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + uq200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQ200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQ200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Tsurface -------------#
	if plot_Tsurface:
		if 'tsurface' not in loaded_vars:
			tsurface = Tsurface(cntrldata, testdata)
			loaded_vars.append('tsurface')
		print("\nPlotting " + which + " Tsurface data...")
		fig = plt.figure()
		fig.suptitle(tsurface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		tsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, tsurface.test, tsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(tsurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, tsurface.control, tsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(tsurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		tsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, tsurface.diff, tsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + tsurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Tsurface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Tsurface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- T850 -------------#
	if plot_T850:
		if 't850' not in loaded_vars:
			t850 = T850(cntrldata, testdata)
			loaded_vars.append('t850')
		print("\nPlotting " + which + " T850 data...")
		fig = plt.figure()
		fig.suptitle(t850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		t850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, t850.test, t850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(t850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, t850.control, t850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(t850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		t850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, t850.diff, t850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + t850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "T850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "T850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- T500 -------------#
	if plot_T500:
		if 't500' not in loaded_vars:
			t500 = T500(cntrldata, testdata)
			loaded_vars.append('t500')
		print("\nPlotting " + which + " T500 data...")
		fig = plt.figure()
		fig.suptitle(t500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		t500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, t500.test, t500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(t500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, t500.control, t500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(t500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		t500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, t500.diff, t500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + t500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "T500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "T500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- T200 -------------#
	if plot_T200:
		if 't200' not in loaded_vars:
			t200 = T200(cntrldata, testdata)
			loaded_vars.append('t200')
		print("\nPlotting " + which + " T200 data...")
		fig = plt.figure()
		fig.suptitle(t200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		t200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, t200.test, t200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(t200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, t200.control, t200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(t200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		t200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, t200.diff, t200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + t200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "T200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "T200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- OMEGAsurface -------------#
	if plot_OMEGAsurface:
		if 'omegasurface' not in loaded_vars:
			omegasurface = OMEGAsurface(cntrldata, testdata)
			loaded_vars.append('omegasurface')
		print("\nPlotting " + which + " OMEGAsurface data...")
		fig = plt.figure()
		fig.suptitle(omegasurface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		omegasurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, omegasurface.test, omegasurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(omegasurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, omegasurface.control, omegasurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(omegasurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		omegasurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, omegasurface.diff, omegasurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + omegasurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "OMEGAsurface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "OMEGAsurface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- OMEGA850 -------------#
	if plot_OMEGA850:
		if 'omega850' not in loaded_vars:
			omega850 = OMEGA850(cntrldata, testdata)
			loaded_vars.append('omega850')
		print("\nPlotting " + which + " OMEGA850 data...")
		fig = plt.figure()
		fig.suptitle(omega850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		omega850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, omega850.test, omega850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(omega850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, omega850.control, omega850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(omega850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		omega850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, omega850.diff, omega850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + omega850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "OMEGA850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "OMEGA850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- OMEGA500 -------------#
	if plot_OMEGA500:
		if 'omega500' not in loaded_vars:
			omega500 = OMEGA500(cntrldata, testdata)
			loaded_vars.append('omega500')
		print("\nPlotting " + which + " OMEGA500 data...")
		fig = plt.figure()
		fig.suptitle(omega500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		omega500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, omega500.test, omega500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(omega500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, omega500.control, omega500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(omega500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		omega500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, omega500.diff, omega500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + omega500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "OMEGA500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "OMEGA500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- OMEGA200 -------------#
	if plot_OMEGA200:
		if 'omega200' not in loaded_vars:
			omega200 = OMEGA200(cntrldata, testdata)
			loaded_vars.append('omega200')
		print("\nPlotting " + which + " OMEGA200 data...")
		fig = plt.figure()
		fig.suptitle(omega200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		omega200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, omega200.test, omega200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(omega200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, omega200.control, omega200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(omega200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		omega200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, omega200.diff, omega200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + omega200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "OMEGA200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "OMEGA200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Z3surface -------------#
	if plot_Z3surface:
		if 'z3surface' not in loaded_vars:
			z3surface = Z3surface(cntrldata, testdata)
			loaded_vars.append('z3surface')
		print("\nPlotting " + which + " Z3surface data...")
		fig = plt.figure()
		fig.suptitle(z3surface.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		z3surface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, z3surface.test, z3surface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(z3surface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, z3surface.control, z3surface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(z3surface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		z3surface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, z3surface.diff, z3surface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + z3surface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Z3surface_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Z3surface_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Z3850 -------------#
	if plot_Z3850:
		if 'z3850' not in loaded_vars:
			z3850 = Z3850(cntrldata, testdata)
			loaded_vars.append('z3850')
		print("\nPlotting " + which + " Z3850 data...")
		fig = plt.figure()
		fig.suptitle(z3850.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		z3850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, z3850.test, z3850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(z3850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, z3850.control, z3850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(z3850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		z3850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, z3850.diff, z3850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + z3850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Z3850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Z3850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Z3500 -------------#
	if plot_Z3500:
		if 'z3500' not in loaded_vars:
			z3500 = Z3500(cntrldata, testdata)
			loaded_vars.append('z3500')
		print("\nPlotting " + which + " Z3500 data...")
		fig = plt.figure()
		fig.suptitle(z3500.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		z3500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, z3500.test, z3500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(z3500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, z3500.control, z3500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(z3500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		z3500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, z3500.diff, z3500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + z3500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Z3500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Z3500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- Z3200 -------------#
	if plot_Z3200:
		if 'z3200' not in loaded_vars:
			z3200 = Z3200(cntrldata, testdata)
			loaded_vars.append('z3200')
		print("\nPlotting " + which + " Z3200 data...")
		fig = plt.figure()
		fig.suptitle(z3200.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		z3200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, z3200.test, z3200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(z3200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, z3200.control, z3200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(z3200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		z3200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, z3200.diff, z3200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + z3200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Z3200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Z3200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECT_d18O -------------#
	if plot_PRECT_d18O:
		if 'prect_h2o' not in loaded_vars:
			prect_h2o = PRECT_H2O(cntrldata, testdata)
			loaded_vars.append('prect_h2o')
		if 'prect_h218o' not in loaded_vars:
			prect_h218o = PRECT_H218O(cntrldata, testdata)
			loaded_vars.append('prect_h218o')
		prect_d18o = PRECT_d18O()
		prect_d18o.calculate(prect_h218o, prect_h2o)
		loaded_vars.append("prect_d18o")
		print("\nPlotting " + which + " PRECT_d18O data...")
		fig = plt.figure()
		fig.suptitle(prect_d18o.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_d18o.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, prect_d18o.test, prect_d18o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_d18o.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, prect_d18o.control, prect_d18o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_d18o.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_d18o.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, prect_d18o.diff, prect_d18o.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + prect_d18o.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT_d18O_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT_d18O_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECT_dD -------------#
	if plot_PRECT_dD:
		if 'prect_h2o' not in loaded_vars:
			prect_h2o = PRECT_H2O(cntrldata, testdata)
			loaded_vars.append('prect_h2o')
		if 'prect_hdo' not in loaded_vars:
			prect_hdo = PRECT_HDO(cntrldata, testdata)
			loaded_vars.append('prect_hdo')
		prect_dd = PRECT_dD()
		prect_dd.calculate(prect_hdo, prect_h2o)
		loaded_vars.append("prect_dd")
		print("\nPlotting " + which + " PRECT_dD data...")
		fig = plt.figure()
		fig.suptitle(prect_dd.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_dd.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, prect_dd.test, prect_dd.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_dd.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, prect_dd.control, prect_dd.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_dd.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_dd.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, prect_dd.diff, prect_dd.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + prect_dd.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT_dD_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT_dD_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- PRECT_dxs -------------#
	if plot_PRECT_dxs:
		if 'prect_d18o' not in loaded_vars:
			if 'prect_h2o' not in loaded_vars:
				prect_h2o = PRECT_H2O(cntrldata, testdata)
				loaded_vars.append('prect_h2o')
			if 'prect_h218o' not in loaded_vars:
				prect_h218o = PRECT_H218O(cntrldata, testdata)
				loaded_vars.append('prect_h218o')
			prect_d18o = PRECT_d18O()
			prect_d18o.calculate(prect_h218o, prect_h2o)
			loaded_vars.append("prect_d18o")
			
		if 'prect_dd' not in loaded_vars:
			if 'prect_hdo' not in loaded_vars:
				prect_hdo = PRECT_HDO(cntrldata, testdata)
				loaded_vars.append('prect_hdo')
			prect_dd = PRECT_dD()
			prect_dd.calculate(prect_hdo, prect_h2o)
			loaded_vars.append("prect_dd")

		prect_dxs = PRECT_dxs()
		prect_dxs.calculate(prect_d18o, prect_dd)

		print("\nPlotting " + which + " PRECT_dxs data...")
		fig = plt.figure()
		fig.suptitle(prect_dxs.long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_dxs.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, prect_dxs.test, prect_dxs.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_dxs.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, prect_dxs.control, prect_dxs.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(prect_dxs.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		prect_dxs.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, prect_dxs.diff, prect_dxs.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + prect_dxs.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT_dxs_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT_dxs_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- d18OVsurface -------------#
	if plot_d18OVsurface:
		if 'h2ovsurface' not in loaded_vars:
			h2ovsurface = H2OVsurface(cntrldata, testdata)
			loaded_vars.append('h2ovsurface')
		if 'h218ovsurface' not in loaded_vars:
			h218ovsurface = H218OVsurface(cntrldata, testdata)
			loaded_vars.append('h218ovsurface')
		d18ovsurface = d18OV()
		d18ovsurface.calculate(h218ovsurface, h2ovsurface)
		loaded_vars.append("d18ovsurface")
		print("\nPlotting " + which + " d18OVsurface data...")
		fig = plt.figure()
		fig.suptitle(d18ovsurface.long_name + " at surface", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		d18ovsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, d18ovsurface.test, d18ovsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18ovsurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18ovsurface.control, d18ovsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18ovsurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		d18ovsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, d18ovsurface.diff, d18ovsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + d18ovsurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "d18OV_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "d18OV_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- d18OV850 -------------#
	if plot_d18OV850:
		if 'h2ov850' not in loaded_vars:
			h2ov850 = H2OV850(cntrldata, testdata)
			loaded_vars.append('h2ov850')
		if 'h218ov850' not in loaded_vars:
			h218ov850 = H218OV850(cntrldata, testdata)
			loaded_vars.append('h218ov850')
		d18ov850 = d18OV()
		d18ov850.calculate(h218ov850, h2ov850)
		loaded_vars.append("d18ov850")
		print("\nPlotting " + which + " d18OV850 data...")
		fig = plt.figure()
		fig.suptitle(d18ov850.long_name + " at 850mb", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		d18ov850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, d18ov850.test, d18ov850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18ov850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18ov850.control, d18ov850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18ov850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		d18ov850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, d18ov850.diff, d18ov850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + d18ov850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "d18OV850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "d18OV850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- d18OV500 -------------#
	if plot_d18OV500:
		if 'h2ov500' not in loaded_vars:
			h2ov500 = H2OV500(cntrldata, testdata)
			loaded_vars.append('h2ov500')
		if 'h218ov500' not in loaded_vars:
			h218ov500 = H218OV500(cntrldata, testdata)
			loaded_vars.append('h218ov500')
		d18ov500 = d18OV()
		d18ov500.calculate(h218ov500, h2ov500)
		loaded_vars.append("d18ov500")
		print("\nPlotting " + which + " d18OV500 data...")
		fig = plt.figure()
		fig.suptitle(d18ov500.long_name + " at 500mb", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		d18ov500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, d18ov500.test, d18ov500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18ov500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18ov500.control, d18ov500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18ov500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		d18ov500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, d18ov500.diff, d18ov500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + d18ov500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "d18OV500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "d18OV500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- d18OV200 -------------#
	if plot_d18OV200:
		if 'h2ov200' not in loaded_vars:
			h2ov200 = H2OV200(cntrldata, testdata)
			loaded_vars.append('h2ov200')
		if 'h218ov200' not in loaded_vars:
			h218ov200 = H218OV200(cntrldata, testdata)
			loaded_vars.append('h218ov200')
		d18ov200 = d18OV()
		d18ov200.calculate(h218ov200, h2ov200)
		loaded_vars.append("d18ov200")
		print("\nPlotting " + which + " d18OV200 data...")
		fig = plt.figure()
		fig.suptitle(d18ov200.long_name + " at 200mb", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		d18ov200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, d18ov200.test, d18ov200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18ov200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18ov200.control, d18ov200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18ov200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		d18ov200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, d18ov200.diff, d18ov200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + d18ov200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "d18OV200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "d18OV200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- dDVsurface -------------#
	if plot_dDVsurface:
		if 'h2ovsurface' not in loaded_vars:
			h2ovsurface = H2OVsurface(cntrldata, testdata)
			loaded_vars.append('h2ovsurface')
		if 'hdovsurface' not in loaded_vars:
			hdovsurface = HDOVsurface(cntrldata, testdata)
			loaded_vars.append('hdovsurface')
		ddvsurface = dDV()
		ddvsurface.calculate(hdovsurface, h2ovsurface)
		loaded_vars.append("ddvsurface")
		print("\nPlotting " + which + " dDVsurface data...")
		fig = plt.figure()
		fig.suptitle(ddvsurface.long_name + " at surface", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ddvsurface.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, ddvsurface.test, ddvsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ddvsurface.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, ddvsurface.control, ddvsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ddvsurface.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ddvsurface.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, ddvsurface.diff, ddvsurface.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + ddvsurface.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "dDV_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "dDV_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- dDV850 -------------#
	if plot_dDV850:
		if 'h2ov850' not in loaded_vars:
			h2ov850 = H2OV850(cntrldata, testdata)
			loaded_vars.append('h2ov850')
		if 'hdov850' not in loaded_vars:
			hdov850 = HDOV850(cntrldata, testdata)
			loaded_vars.append('hdov850')
		ddv850 = dDV()
		ddv850.calculate(hdov850, h2ov850)
		loaded_vars.append("ddv850")
		print("\nPlotting " + which + " dDV850 data...")
		fig = plt.figure()
		fig.suptitle(ddv850.long_name + " at 850mb", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ddv850.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, ddv850.test, ddv850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ddv850.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, ddv850.control, ddv850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ddv850.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ddv850.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, ddv850.diff, ddv850.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + ddv850.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "dDV850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "dDV850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- dDV500 -------------#
	if plot_dDV500:
		if 'h2ov500' not in loaded_vars:
			h2ov500 = H2OV500(cntrldata, testdata)
			loaded_vars.append('h2ov500')
		if 'hdov500' not in loaded_vars:
			hdov500 = HDOV500(cntrldata, testdata)
			loaded_vars.append('hdov500')
		ddv500 = dDV()
		ddv500.calculate(hdov500, h2ov500)
		loaded_vars.append("ddv500")
		print("\nPlotting " + which + " dDV500 data...")
		fig = plt.figure()
		fig.suptitle(ddv500.long_name + " at 500mb", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ddv500.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, ddv500.test, ddv500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ddv500.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, ddv500.control, ddv500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ddv500.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ddv500.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, ddv500.diff, ddv500.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + ddv500.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "dDV500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "dDV500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#

	#------------- dDV200 -------------#
	if plot_dDV200:
		if 'h2ov200' not in loaded_vars:
			h2ov200 = H2OV200(cntrldata, testdata)
			loaded_vars.append('h2ov200')
		if 'hdov200' not in loaded_vars:
			hdov200 = HDOV200(cntrldata, testdata)
			loaded_vars.append('hdov200')
		ddv200 = dDV()
		ddv200.calculate(hdov200, h2ov200)
		loaded_vars.append("ddv200")
		print("\nPlotting " + which + " dDV200 data...")
		fig = plt.figure()
		fig.suptitle(ddv200.long_name + " at 200mb", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ddv200.set_clevs(which+"_"+region)
		cs = m.contourf(bmlon, bmlat, ddv200.test, ddv200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ddv200.units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, ddv200.control, ddv200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ddv200.units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		ddv200.set_clevs(which+"diff_"+region)
		cs = m.contourf(bmlon, bmlat, ddv200.diff, ddv200.clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + ddv200.units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "dDV200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "dDV200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
	#---------------------------------#



# The following writes the directory that the figures are stored in to a hidden text file on the user's home directory
# This can be used in, for example, a shell script, to read where the figures are located without having to hard code it in
# E.g.: image_dir='cat ~/.tempfile.txt' in bash will store the image directory in the env variable called $image_dir
home = os.path.expanduser("~")
figdirtemp = open(home+"/.figdirtemp.txt", "w")
figdirtemp.write(figdir)

regiontemp = open(home+"/.regiontemp.txt", "w")
regiontemp.write(region_name)


