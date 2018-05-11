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
from UsefulNetcdfFunctions import dR_dt, ExtractVarAtPressure

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
region = "NA"

#------------- Set seasons 
#seasons = ["ANN", "DJF", "JJA", "JAS"]
seasons = ["JAS"]

#------------- Which variables to plot?
PRECT850windbarbs = False
PRECT500windbarbs =  False
PRECT200windbarbs =  False
dRHDO_dt =  False
PRECT =  True
PRECC  =  False
PRECL  =  False
PRECT_dD =  False
PRECT_d18O =  False
PRECT_dxs =  False
PRECT_H2O =  False
PRECST =  False
PRECCfrac =  False
PRECLfrac =  False
Qlowest =  False
Q850 = False
Q500 = False
Q200 = False
dDVlowest =  False
dDVmiddle =  False
d18OVlowest =  False
d18OVmiddle =  False
dxsVlowest =  False
dxsVmiddle =  False
QFLX =  False
QFLX_HDO =  False
QFLX_H218O =  False
VQ_H2Olowest =  False
VQ_H2Omiddle =  False
VQ_dDlowest =  False
VQ_dDmiddle =  False
VQ_d18Olowest =  False
VQ_d18Omiddle =  False
UQ_H2Olowest =  False
UQ_H2Omiddle =  False
UQ_dDlowest =  False
UQ_dDmiddle =  False
UQ_d18Olowest =  False
UQ_d18Omiddle =  False
VTlowest =  False
VTmiddle =  False
UTlowest =  False
UTmiddle =  False
TTEND_TOTlowest =  False
TTEND_TOTmiddle =  False
Ulowest =  False
U850 = False
U500 = False
U200 = False
Vlowest =  False
V850 = False
V500 = False
V200 = False
Tlowest =  False
T850 = False
T500 = False
T200 = False
OMEGAlowest =  False
OMEGA850 = False
OMEGA500 = False
OMEGA200 = False
Z3850 = False
Z3500 = False
Z3200 = False
CLDHGH =  False
PS =  False
PSL =  False
PSI500 =  False
RH850 = False
RH500 = False
RH200 = False
RHlower =  False

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
	bmlon, bmlat = np.meshgrid(lons, lats)
	levs = cntrldata.variables['lev'][:]
	A = cntrldata.variables['hyam'][:]
	B = cntrldata.variables['hybm'][:]
	
	g = -9.8 # gravitational constant
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	print "\nExtracting 2D data from netcdf files..."
	#------------- Extract numerical data -------------#
	####################
	### 2D variables ###
	####################
	# Precipitation
	PRECT_control = cntrldata.variables['PRECT'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_test = testdata.variables['PRECT'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_units = "mm/day"
	PRECT_long_name = "Total precipitation rate"

	PRECT_H2O_control = cntrldata.variables['PRECT_H2O'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_H2O_test = testdata.variables['PRECT_H2O'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_H2O_units = "mm/day"
	PRECT_H2O_long_name = cntrldata.variables['PRECT_H2O'].long_name

	PRECT_H218O_control = cntrldata.variables['PRECT_H218O'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_H218O_test = testdata.variables['PRECT_H218O'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_H218O_units = "mm/day"
	PRECT_H218O_long_name = cntrldata.variables['PRECT_H218O'].long_name

	PRECT_d18O_control = (PRECT_H218O_control/PRECT_H2O_control - 1) * 1000
	PRECT_d18O_test = (PRECT_H218O_test/PRECT_H2O_test - 1) * 1000
	PRECT_d18O_units = "delta 18O (permil)"
	PRECT_d18O_long_name = "delta 18O for PRECT"

	PRECT_HDO_control = cntrldata.variables['PRECT_HDO'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_HDO_test = testdata.variables['PRECT_HDO'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_HDO_units = "mm/day"
	PRECT_HDO_long_name = cntrldata.variables['PRECT_HDO'].long_name
	
	PRECT_dD_control = (PRECT_HDO_control/PRECT_H2O_control - 1) * 1000
	PRECT_dD_test = (PRECT_HDO_test/PRECT_H2O_test - 1) * 1000
	PRECT_dD_units = "delta D (permil)"
	PRECT_dD_long_name = "delta D for PRECT"

	PRECT_dxs_control = PRECT_dD_control - 8*PRECT_d18O_control
	PRECT_dxs_test = PRECT_dD_test - 8*PRECT_d18O_test
	PRECT_dxs_units = "d-excess (permil)"
	PRECT_dxs_long_name = "d-excess (delta_D - 8*delta_18O) for PRECT"

	PRECC_control = cntrldata.variables['PRECC'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECC_test = testdata.variables['PRECC'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECC_units = "mm/day"
	PRECC_long_name = testdata.variables['PRECC'].long_name

	PRECL_control = cntrldata.variables['PRECL'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECL_test = testdata.variables['PRECL'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECL_units = "mm/day"
	PRECL_long_name = cntrldata.variables['PRECL'].long_name

	PRECSC_control = cntrldata.variables['PRECSC'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECSC_test = testdata.variables['PRECSC'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECSC_units = "mm/day"
	PRECSC_long_name = cntrldata.variables['PRECSC'].long_name

	PRECSL_control = cntrldata.variables['PRECSL'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECSL_test = testdata.variables['PRECSL'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECSL_units = "mm/day"
	PRECSL_long_name = cntrldata.variables['PRECSL'].long_name

	PRECST_control = PRECSC_control + PRECSL_control
	PRECST_test = PRECSC_test + PRECSL_test
	PRECST_units = "mm/day"
	PRECST_long_name = "Total (convective + large-scale) snow rate"
	
	QFLX_control = cntrldata.variables['QFLX'][0,:,:].squeeze() * 1000
	QFLX_test = testdata.variables['QFLX'][0,:,:].squeeze() * 1000
	QFLX_units = "g-H2O/m2/s"
	QFLX_long_name = cntrldata.variables['QFLX'].long_name

	QFLX_H218O_control = cntrldata.variables['QFLX_H218O'][0,:,:].squeeze()
	QFLX_H218O_test = testdata.variables['QFLX_H218O'][0,:,:].squeeze()
	QFLX_H218O_units = testdata.variables['QFLX_H218O'].units
	QFLX_H218O_long_name = cntrldata.variables['QFLX_H218O'].long_name

	QFLX_HDO_control = cntrldata.variables['QFLX_HDO'][0,:,:].squeeze()
	QFLX_HDO_test = testdata.variables['QFLX_HDO'][0,:,:].squeeze()
	QFLX_HDO_units = testdata.variables['QFLX_HDO'].units
	QFLX_HDO_long_name = cntrldata.variables['QFLX_HDO'].long_name

	CLDHGH_control = cntrldata.variables['CLDHGH'][0,:,:].squeeze()
	CLDHGH_test = testdata.variables['CLDHGH'][0,:,:].squeeze()
	CLDHGH_units = testdata.variables['CLDHGH'].units
	CLDHGH_long_name = testdata.variables['CLDHGH'].long_name

	PS_control = cntrldata.variables['PS'][0,:,:].squeeze()
	PS_test = testdata.variables['PS'][0,:,:].squeeze()
	PS_units = testdata.variables['PS'].units
	PS_long_name = cntrldata.variables['PS'].long_name

	PSL_control = cntrldata.variables['PSL'][0,:,:].squeeze()
	PSL_test = testdata.variables['PSL'][0,:,:].squeeze()
	PSL_units = testdata.variables['PSL'].units
	PSL_long_name = cntrldata.variables['PSL'].long_name

	print "Extracting 3D data from netcdf files..."
	print "This may take some time..."
	####################
	### 3D variables ###
	####################
	H2OV_control = cntrldata.variables['H2OV'][0,:,:,:].squeeze()
	H2OV_test = testdata.variables['H2OV'][0,:,:,:].squeeze()
	H2OV_units = cntrldata.variables['H2OV'].units
	H2OV_long_name = cntrldata.variables['H2OV'].long_name

	HDOV_control = cntrldata.variables['HDOV'][0,:,:,:].squeeze()
	HDOV_test = testdata.variables['HDOV'][0,:,:,:].squeeze()
	HDOV_units = cntrldata.variables['HDOV'].units
	HDOV_long_name = cntrldata.variables['HDOV'].long_name

	H218OV_control = cntrldata.variables['H218OV'][0,:,:,:].squeeze()
	H218OV_test = testdata.variables['H218OV'][0,:,:,:].squeeze()
	H218OV_units = cntrldata.variables['H218OV'].units
	H218OV_long_name = cntrldata.variables['H218OV'].long_name
	
	dDV_control = (HDOV_control/H2OV_control - 1) * 1000
	dDV_test = (HDOV_test/H2OV_test - 1) * 1000
	dDV_units = "permil"
	dDV_long_name = "delta D for VAPOR"

	d18OV_control = (H218OV_control/H2OV_control - 1) * 1000
	d18OV_test = (H218OV_test/H2OV_test - 1) * 1000
	d18OV_units = "permil"
	d18OV_long_name = "delta 18O for VAPOR"

	dxsV_control = dDV_control - 8*d18OV_control
	dxsV_test = dDV_test - 8*d18OV_test
	dxsV_units = "d-excess (permil)"
	dxsV_long_name = "d-excess (delta_D - 8*delta_18O) for VAPOR"

	Q_control = cntrldata.variables['Q'][0,:,:,:].squeeze() * 100
	Q_test = testdata.variables['Q'][0,:,:,:].squeeze() * 100
	Q_units = "g-H2O/kg"
	Q_long_name = cntrldata.variables['Q'].long_name

	Q850_control = ExtractVarAtPressure(cntrldata, "Q", 85000) * 100
	Q850_test = ExtractVarAtPressure(testdata, "Q", 85000) * 100
	Q850_units = Q_units
	Q850_long_name = Q_long_name + " at 850 mb"

	Q500_control = ExtractVarAtPressure(cntrldata, "Q", 50000) * 100
	Q500_test = ExtractVarAtPressure(testdata, "Q", 50000) * 100
	Q500_units = Q_units
	Q500_long_name = Q_long_name + " at 500 mb"

	Q200_control = ExtractVarAtPressure(cntrldata, "Q", 20000) * 100
	Q200_test = ExtractVarAtPressure(testdata, "Q", 20000) * 100
	Q200_units = Q_units
	Q200_long_name = Q_long_name + " at 200 mb"

	VQ_H2O_control = cntrldata.variables['VQ_H2O'][0,:,:,:].squeeze()
	VQ_H2O_test = testdata.variables['VQ_H2O'][0,:,:,:].squeeze()
	VQ_H2O_units = testdata.variables['VQ_H2O'].units
	VQ_H2O_long_name = cntrldata.variables['VQ_H2O'].long_name

	VQ_H218O_control = cntrldata.variables['VQ_H218O'][0,:,:,:].squeeze()
	VQ_H218O_test = testdata.variables['VQ_H218O'][0,:,:,:].squeeze()
	VQ_H218O_units = testdata.variables['VQ_H218O'].units
	VQ_H218O_long_name = cntrldata.variables['VQ_H218O'].long_name

	VQ_d18O_control = (VQ_H218O_control/VQ_H2O_control - 1) * 1000
	VQ_d18O_test = (VQ_H218O_test/VQ_H2O_test - 1) * 1000
	VQ_d18O_units = VQ_H218O_units
	VQ_d18O_long_name = VQ_H218O_long_name

	VQ_HDO_control = cntrldata.variables['VQ_HDO'][0,:,:,:].squeeze()
	VQ_HDO_test = testdata.variables['VQ_HDO'][0,:,:,:].squeeze()
	VQ_HDO_units = testdata.variables['VQ_HDO'].units
	VQ_HDO_long_name = cntrldata.variables['VQ_HDO'].long_name

	VQ_dD_control = (VQ_HDO_control/VQ_H2O_control - 1) * 1000
	VQ_dD_test = (VQ_HDO_test/VQ_H2O_test - 1) * 1000
	VQ_dD_units = VQ_HDO_units
	VQ_dD_long_name = VQ_HDO_long_name

	UQ_H2O_control = cntrldata.variables['UQ_H2O'][0,:,:,:].squeeze()
	UQ_H2O_test = testdata.variables['UQ_H2O'][0,:,:,:].squeeze()
	UQ_H2O_units = testdata.variables['UQ_H2O'].units
	UQ_H2O_long_name = "Zonal flux for H2O" # Had to hard code in these names because the netcdf files have the wrong long_name ---> need to correct this in the cam_diagnostics.f90 file

	UQ_H218O_control = cntrldata.variables['UQ_H218O'][0,:,:,:].squeeze()
	UQ_H218O_test = testdata.variables['UQ_H218O'][0,:,:,:].squeeze()
	UQ_H218O_units = testdata.variables['UQ_H218O'].units
	UQ_H218O_long_name = "Zonal flux for H218O"

	UQ_d18O_control = (UQ_H218O_control/UQ_H2O_control - 1) * 1000
	UQ_d18O_test = (UQ_H218O_test/UQ_H2O_test - 1) * 1000
	UQ_d18O_units = UQ_H218O_units
	UQ_d18O_long_name = "Zonal flux for H218O"

	UQ_HDO_control = cntrldata.variables['UQ_HDO'][0,:,:,:].squeeze()
	UQ_HDO_test = testdata.variables['UQ_HDO'][0,:,:,:].squeeze()
	UQ_HDO_units = testdata.variables['UQ_HDO'].units
	UQ_HDO_long_name = "Zonal flux for HDO"

	UQ_dD_control = (UQ_HDO_control/VQ_H2O_control - 1) * 1000
	UQ_dD_test = (UQ_HDO_test/UQ_H2O_test - 1) * 1000
	UQ_dD_units = UQ_HDO_units
	UQ_dD_long_name = "Zonal flux for HDO"

	V_control = cntrldata.variables['V'][0,:,:,:].squeeze()
	V_test = testdata.variables['V'][0,:,:,:].squeeze()
	V_units = testdata.variables['V'].units
	V_long_name = cntrldata.variables['V'].long_name

	V850_control = ExtractVarAtPressure(cntrldata, "V", 85000)
	V850_test = ExtractVarAtPressure(testdata, "V", 85000)
	V850_units = V_units
	V850_long_name = "Meridional wind speed at 850 mb"

	V500_control = ExtractVarAtPressure(cntrldata, "V", 50000)
	V500_test = ExtractVarAtPressure(testdata, "V", 50000)
	V500_units = V_units
	V500_long_name = "Meridional wind speed at 500 mb"

	V200_control = ExtractVarAtPressure(cntrldata, "V", 20000)
	V200_test = ExtractVarAtPressure(testdata, "V", 20000)
	V200_units = V_units
	V200_long_name = "Meridional wind speed at 200 mb"

	VT_control = cntrldata.variables['VT'][0,:,:,:].squeeze()
	VT_test = testdata.variables['VT'][0,:,:,:].squeeze()
	VT_units = testdata.variables['VT'].units
	VT_long_name = cntrldata.variables['VT'].long_name

	U_control = cntrldata.variables['U'][0,:,:,:].squeeze()
	U_test = testdata.variables['U'][0,:,:,:].squeeze()
	U_units = testdata.variables['U'].units
	U_long_name = cntrldata.variables['U'].long_name
	
	U850_control = ExtractVarAtPressure(cntrldata, "U", 85000)
	U850_test = ExtractVarAtPressure(testdata, "U", 85000)
	U850_units = U_units
	U850_long_name = "Zonal wind speed at 850 mb"

	U500_control = ExtractVarAtPressure(cntrldata, "U", 50000)
	U500_test = ExtractVarAtPressure(testdata, "U", 50000)
	U500_units = U_units
	U500_long_name = "Zonal wind speed at 500 mb"

	U200_control = ExtractVarAtPressure(cntrldata, "U", 20000)
	U200_test = ExtractVarAtPressure(testdata, "U", 20000)
	U200_units = U_units
	U200_long_name = "Zonal wind speed at 200 mb"

	UT_control = cntrldata.variables['UT'][0,:,:,:].squeeze()
	UT_test = testdata.variables['UT'][0,:,:,:].squeeze()
	UT_units = testdata.variables['UT'].units
	UT_long_name = cntrldata.variables['UT'].long_name

	T_control = cntrldata.variables['T'][0,:,:,:].squeeze()
	T_test = testdata.variables['T'][0,:,:,:].squeeze()
	T_units = testdata.variables['T'].units
	T_long_name = cntrldata.variables['T'].long_name
	
	T850_control = ExtractVarAtPressure(cntrldata, "T", 85000)
	T850_test = ExtractVarAtPressure(testdata, "T", 85000)
	T850_units = T_units
	T850_long_name = "Temperature at 850 mb"

	T500_control = ExtractVarAtPressure(cntrldata, "T", 50000)
	T500_test = ExtractVarAtPressure(testdata, "T", 50000)
	T500_units = T_units
	T500_long_name = "Temperature at 500 mb"

	T200_control = ExtractVarAtPressure(cntrldata, "T", 20000)
	T200_test = ExtractVarAtPressure(testdata, "T", 20000)
	T200_units = T_units
	T200_long_name = "Temperature at 200 mb"
	
	print "Almost done..."

	OMEGA_control = cntrldata.variables['OMEGA'][0,:,:,:].squeeze()
	OMEGA_test = testdata.variables['OMEGA'][0,:,:,:].squeeze()
	OMEGA_units = testdata.variables['OMEGA'].units
	OMEGA_long_name = cntrldata.variables['OMEGA'].long_name
	
	OMEGA850_control = ExtractVarAtPressure(cntrldata, "OMEGA", 85000)
	OMEGA850_test = ExtractVarAtPressure(testdata, "OMEGA", 85000)
	OMEGA850_units = OMEGA_units
	OMEGA850_long_name = OMEGA_long_name + " at 850 mb"

	OMEGA500_control = ExtractVarAtPressure(cntrldata, "OMEGA", 50000)
	OMEGA500_test = ExtractVarAtPressure(testdata, "OMEGA", 50000)
	OMEGA500_units = OMEGA_units
	OMEGA500_long_name = OMEGA_long_name + " at 500 mb"

	OMEGA200_control = ExtractVarAtPressure(cntrldata, "OMEGA", 20000)
	OMEGA200_test = ExtractVarAtPressure(testdata, "OMEGA", 20000)
	OMEGA200_units = OMEGA_units
	OMEGA200_long_name = OMEGA_long_name + " at 200 mb"
	
	Z3_control = cntrldata.variables['Z3'][0,:,:,:].squeeze() / 1000 # Convert to km
	Z3_test = testdata.variables['Z3'][0,:,:,:].squeeze() / 1000 
	Z3_units = "km"
	Z3_long_name = cntrldata.variables['Z3'].long_name

	Z3850_control = ExtractVarAtPressure(cntrldata, "Z3", 85000) / 1000 # Convert to km
	Z3850_test = ExtractVarAtPressure(testdata, "Z3", 85000) / 1000 # Convert to km
	Z3850_units = Z3_units
	Z3850_long_name = Z3_long_name + " at 850 mb"

	Z3500_control = ExtractVarAtPressure(cntrldata, "Z3", 50000) / 1000 # Convert to km
	Z3500_test = ExtractVarAtPressure(testdata, "Z3", 50000) / 1000 # Convert to km
	Z3500_units = Z3_units
	Z3500_long_name = Z3_long_name + " at 500 mb"

	Z3200_control = ExtractVarAtPressure(cntrldata, "Z3", 20000) / 1000 # Convert to km
	Z3200_test = ExtractVarAtPressure(testdata, "Z3", 20000) / 1000 # Convert to km
	Z3200_units = Z3_units
	Z3200_long_name = Z3_long_name + " at 200 mb"

	P_control = []
	P_test = []
	for n in range(len(B)):
		P_control.append(np.array(A[n] * 100000 + B[n] * PS_control))
		P_test.append(np.array(A[n] * 100000 + B[n] * PS_test))
	P_control = np.array(P_control)
	P_test = np.array(P_test)
	P_units = "Pa"
	P_long_name = "Pressure"

	E_control = Q_control/1000 * P_control / 0.622 # Vapor pressure, q*P/epsilon
	E_test = Q_test/1000 * P_test / 0.622
	E_units = "Pa"
	E_long_name = "Vapor pressure"
	
	E850_control = Q850_control/1000 * 85000 / 0.622 # Vapor pressure, q*P/epsilon
	E850_test = Q850_test/1000 * 85000 / 0.622 # Vapor pressure, q*P/epsilon
	E850_units = E_units
	E850_long_name = E_long_name + " at 850 mb"

	E500_control = Q500_control/1000 * 50000 / 0.622 # Vapor pressure, q*P/epsilon
	E500_test = Q500_test/1000 * 50000 / 0.622 # Vapor pressure, q*P/epsilon
	E500_units = E_units
	E500_long_name = E_long_name + " at 500 mb"

	E200_control = Q200_control/1000 * 20000 / 0.622 # Vapor pressure, q*P/epsilon
	E200_test = Q200_test/1000 * 20000 / 0.622 # Vapor pressure, q*P/epsilon
	E200_units = E_units
	E200_long_name = E_long_name + " at 200 mb"

	ES_control = 2.53e8 * np.exp(-5.42e3 / T_control) # saturation vapor pressure, Ae^(-B/T)
	ES_test = 2.53e8 * np.exp(-5.42e3 / T_test)
	ES_units = "Pa"
	ES_long_name = "Saturation vapor pressure"

	ES850_control = 2.53e8 * np.exp(-5.42e3 / T850_control) # saturation vapor pressure, Ae^(-B/T)
	ES850_test = 2.53e8 * np.exp(-5.42e3 / T850_test)
	ES850_units = "Pa"
	ES850_long_name = "Saturation vapor pressure"

	ES500_control = 2.53e8 * np.exp(-5.42e3 / T500_control) # saturation vapor pressure, Ae^(-B/T)
	ES500_test = 2.53e8 * np.exp(-5.42e3 / T500_test)
	ES500_units = "Pa"
	ES500_long_name = "Saturation vapor pressure"

	ES200_control = 2.53e8 * np.exp(-5.42e3 / T200_control) # saturation vapor pressure, Ae^(-B/T)
	ES200_test = 2.53e8 * np.exp(-5.42e3 / T200_test)
	ES200_units = "Pa"
	ES200_long_name = "Saturation vapor pressure"

	RH_control = E_control / ES_control
	RH_test = E_test / ES_test
	RH_units = "percent"
	RH_long_name = "Relative humidity"
	
	RH850_control = E850_control / ES850_control
	RH850_test = E850_test / ES850_test
	RH850_units = "percent"
	RH850_long_name = "Relative humidity at 850 mb"

	RH500_control = E500_control / ES500_control
	RH500_test = E500_test / ES500_test
	RH500_units = "percent"
	RH500_long_name = "Relative humidity at 500 mb"

	RH200_control = E200_control / ES200_control
	RH200_test = E200_test / ES200_test
	RH200_units = "percent"
	RH200_long_name = "Relative humidity at 200 mb"

	# Tendencies
	TTEND_TOT_control = cntrldata.variables['TTEND_TOT'][0,:,:,:].squeeze() * 60 * 60 * 24 * 365 # Convert to K/yr
	TTEND_TOT_test = testdata.variables['TTEND_TOT'][0,:,:,:].squeeze() * 60 * 60 * 24 * 365
	TTEND_TOT_units = "K/year"
	TTEND_TOT_long_name = cntrldata.variables['TTEND_TOT'].long_name

	DCQ_control = cntrldata.variables['DCQ'][0,:,:,:].squeeze()
	DCQ_test = testdata.variables['DCQ'][0,:,:,:].squeeze()
	DCQ_units = testdata.variables['DCQ'].units
	DCQ_long_name = cntrldata.variables['DCQ'].long_name

	DCQ_control = cntrldata.variables['DCQ'][0,:,:,:].squeeze()
	DCQ_test = testdata.variables['DCQ'][0,:,:,:].squeeze()
	DCQ_units = testdata.variables['DCQ'].units
	DCQ_long_name = cntrldata.variables['DCQ'].long_name

	CMFDT_control = cntrldata.variables['CMFDT'][0,:,:,:].squeeze()
	CMFDT_test = testdata.variables['CMFDT'][0,:,:,:].squeeze()
	CMFDT_units = testdata.variables['CMFDT'].units
	CMFDT_long_name = cntrldata.variables['CMFDT'].long_name

	CMFDQR_control = cntrldata.variables['CMFDQR'][0,:,:,:].squeeze()
	CMFDQR_test = testdata.variables['CMFDQR'][0,:,:,:].squeeze()
	CMFDQR_units = testdata.variables['CMFDQR'].units
	CMFDQR_long_name = cntrldata.variables['CMFDQR'].long_name

	DTCOND_control = cntrldata.variables['DTCOND'][0,:,:,:].squeeze()
	DTCOND_test = testdata.variables['DTCOND'][0,:,:,:].squeeze()
	DTCOND_units = testdata.variables['DTCOND'].units
	DTCOND_long_name = cntrldata.variables['DTCOND'].long_name
	
	dRHDO_dt_control = dR_dt(cntrldatafname, "HDO")
	dRHDO_dt_test = dR_dt(testdatafname, "HDO")
	dRHDO_dt_units = "R / s"
	dRHDO_dt_long_name = "Surface HDOV tendency due to horizontal advection"

	# Psi calcs (from psi500_grid.py - see this file for notes)
	A_control = cntrldata.variables['hyam'][:]
	A_test = testdata.variables['hyam'][:]
	B_control = cntrldata.variables['hybm'][:]
	B_test = testdata.variables['hybm'][:]
	
	def psi500calc(psi, P):
		pressures_above500mb = P - 50000 < 0
		pressures_below500mb = P - 50000 > 0		
		last_pressure_before500mb = P[pressures_above500mb][-1]
		first_pressure_after500mb = P[pressures_below500mb][0]
		fractional_level_index = np.argmax(P>50000)
		f = (50000 - last_pressure_before500mb) / abs(first_pressure_after500mb - last_pressure_before500mb)
		fractional_addition = (f * psi[fractional_level_index])
		return np.sum(psi[pressures_above500mb]) + fractional_addition

	# Control
	pressure_control = []
	for n in range(len(B_control)):
		pressure_control.append(np.array(A_control[n] * 100000 + B_control[n] * PS_control))
	pressure_control = np.array(pressure_control)
	psi_before_latnormalized = np.multiply(pressure_control, V_control) / g
	psi = []
	for n, lat in enumerate(lats):
		psi.append(psi_before_latnormalized[:,n,:] * np.cos(lat*np.pi/180))
	psi = np.swapaxes(np.array(psi), 0, 1)
	psi500 = []
	for i in range(len(lats)):
		for j in range(len(lons)):
			if j == 0:
				hold = []
			hold.append(psi500calc(psi[:,i,j], pressure_control[:,i,j]))
		psi500.append(hold)
		hold = []
	PSI500_control = np.array(psi500)
	# Test
	pressure_test = []
	for n in range(len(B_test)):
		pressure_test.append(np.array(A_test[n] * 100000 + B_test[n] * PS_test))
	pressure_test = np.array(pressure_test)
	psi_before_latnormalized = np.multiply(pressure_test, V_test) / g
	psi = []
	for n, lat in enumerate(lats):
		psi.append(psi_before_latnormalized[:,n,:] * np.cos(lat*np.pi/180))
	psi = np.swapaxes(np.array(psi), 0, 1)
	psi500 = []
	for i in range(len(lats)):
		for j in range(len(lons)):
			if j == 0:
				hold = []
			hold.append(psi500calc(psi[:,i,j], pressure_test[:,i,j]))
		psi500.append(hold)
		hold = []
	PSI500_test = np.array(psi500)
	PSI500_units = "meridional mass transport"
	PSI500_long_name = "Meridional mass stream function, integrated from model top to 500mb"

	print "All data extracted!"

	def find_clevs(mat, n = 21):
		MIN = np.min(mat)
		MAX = np.max(mat)
		return np.linspace(MIN, MAX, n)


	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------- Map creation ---------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if PRECT850windbarbs:
		#------------- PRECT with 850mb wind barbs -------------#
		density = 3

		print("\nPlotting " + which + " PRECT with 850mb winds data...")
		fig = plt.figure()
		fig.suptitle(PRECT_long_name + " with 850mb winds", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(0, 15, 19)
		clevs2 = np.linspace(-4,4,17)
		if which == "JJA":
			if region == "NA":
				clevs = np.linspace(0, 15, 19)
			if region == "IM":
				clevs = np.linspace(0, 30, 16)
				clevs2 = np.linspace(-8,8,17)
			if region == "MC":
				clevs2 = clevs2 = np.linspace(-6,6,13)
			if region == "EP":
				clevs = np.linspace(0,25,16)
				clevs2 = np.linspace(-6,6,13)
		if which == "DJF":
			if region == "NA":
				clevs = np.linspace(0, 15, 19)
			if region == "IM":
				clevs = np.linspace(0, 30, 16)
				clevs2 = np.linspace(-8,8,17)
			if region == "MC":
				clevs2 = clevs2 = np.linspace(-6,6,13)
			if region == "EP":
				clevs = np.linspace(0,25,16)
				clevs2 = np.linspace(-6,6,13)
		cs = m.contourf(bmlon, bmlat, PRECT_test, clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)

		colsby4 = np.arange(0, bmlon.shape[1], density)
		rowsby4 = np.arange(0, bmlat.shape[0], density)
		barb_lon = bmlon[:,colsby4]
		barb_lon = barb_lon[rowsby4,:]
		barb_lat = bmlat[:,colsby4]
		barb_lat = barb_lat[rowsby4,:]
		barbu_test = U850_test[rowsby4,:]
		barbu_test = barbu_test[:,colsby4] 
		barbv_test = V850_test[rowsby4,:]
		barbv_test = barbv_test[:,colsby4] 
		m.quiver(barb_lon, barb_lat, barbu_test, barbv_test, scale = 50, scale_units = "inches", cmap=plt.cm.autumn, latlon = True)

		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_control, clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)

		barbu_control = U850_control[rowsby4,:]
		barbu_control = barbu_control[:,colsby4] 
		barbv_control = V850_control[rowsby4,:]
		barbv_control = barbv_control[:,colsby4] 
		m.quiver(barb_lon, barb_lat, barbu_control, barbv_control, scale = 50, scale_units = "inches", cmap=plt.cm.autumn, latlon = True)

		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_test-PRECT_control, clevs2, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)

		barbu_diff = barbu_test - barbu_control
		barbv_diff = barbv_test - barbv_control
		m.quiver(barb_lon, barb_lat, barbu_diff, barbv_diff, scale = 12.5, scale_units = "inches", cmap=plt.cm.autumn, latlon = True)

		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT850windbarbs_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT850windbarb_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if PRECT500windbarbs:
		#------------- PRECT with 500mb wind barbs -------------#
		density = 3

		print("\nPlotting " + which + " PRECT with 500mb winds data...")
		fig = plt.figure()
		fig.suptitle(PRECT_long_name + " with 500mb winds", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(0, 15, 19)
		clevs2 = np.linspace(-4,4,17)
		if which == "JJA":
			if region == "NA":
				clevs = np.linspace(0, 15, 19)
			if region == "IM":
				clevs = np.linspace(0, 30, 16)
				clevs2 = np.linspace(-8,8,17)
			if region == "MC":
				clevs2 = clevs2 = np.linspace(-6,6,13)
			if region == "EP":
				clevs = np.linspace(0,25,16)
				clevs2 = np.linspace(-6,6,13)
		if which == "DJF":
			if region == "NA":
				clevs = np.linspace(0, 15, 19)
			if region == "IM":
				clevs = np.linspace(0, 30, 16)
				clevs2 = np.linspace(-8,8,17)
			if region == "MC":
				clevs2 = clevs2 = np.linspace(-6,6,13)
			if region == "EP":
				clevs = np.linspace(0,25,16)
				clevs2 = np.linspace(-6,6,13)
		cs = m.contourf(bmlon, bmlat, PRECT_test, clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)

		colsby4 = np.arange(0, bmlon.shape[1], density)
		rowsby4 = np.arange(0, bmlat.shape[0], density)
		barb_lon = bmlon[:,colsby4]
		barb_lon = barb_lon[rowsby4,:]
		barb_lat = bmlat[:,colsby4]
		barb_lat = barb_lat[rowsby4,:]
		barbu_test = U500_test[rowsby4,:]
		barbu_test = barbu_test[:,colsby4] 
		barbv_test = V500_test[rowsby4,:]
		barbv_test = barbv_test[:,colsby4] 
		m.quiver(barb_lon, barb_lat, barbu_test, barbv_test, scale = 50, scale_units = "inches", cmap=plt.cm.autumn, latlon = True)

		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_control, clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)

		barbu_control = U500_control[rowsby4,:]
		barbu_control = barbu_control[:,colsby4] 
		barbv_control = V500_control[rowsby4,:]
		barbv_control = barbv_control[:,colsby4] 
		m.quiver(barb_lon, barb_lat, barbu_control, barbv_control, scale = 50, scale_units = "inches", cmap=plt.cm.autumn, latlon = True)

		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_test-PRECT_control, clevs2, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)

		barbu_diff = barbu_test - barbu_control
		barbv_diff = barbv_test - barbv_control
		m.quiver(barb_lon, barb_lat, barbu_diff, barbv_diff, scale = 12.5, scale_units = "inches", cmap=plt.cm.autumn, latlon = True)

		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT500windbarbs_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT500windbarb_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if PRECT200windbarbs:
		#------------- PRECT with 200mb wind barbs -------------#
		density = 3

		print("\nPlotting " + which + " PRECT with 200mb winds data...")
		fig = plt.figure()
		fig.suptitle(PRECT_long_name + " with 200mb winds", fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(0, 15, 19)
		clevs2 = np.linspace(-4,4,17)
		if which == "JJA":
			if region == "NA":
				clevs = np.linspace(0, 15, 19)
			if region == "IM":
				clevs = np.linspace(0, 30, 16)
				clevs2 = np.linspace(-8,8,17)
			if region == "MC":
				clevs2 = clevs2 = np.linspace(-6,6,13)
			if region == "EP":
				clevs = np.linspace(0,25,16)
				clevs2 = np.linspace(-6,6,13)
		if which == "DJF":
			if region == "NA":
				clevs = np.linspace(0, 15, 19)
			if region == "IM":
				clevs = np.linspace(0, 30, 16)
				clevs2 = np.linspace(-8,8,17)
			if region == "MC":
				clevs2 = clevs2 = np.linspace(-6,6,13)
			if region == "EP":
				clevs = np.linspace(0,25,16)
				clevs2 = np.linspace(-6,6,13)
		cs = m.contourf(bmlon, bmlat, PRECT_test, clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)

		colsby4 = np.arange(0, bmlon.shape[1], density)
		rowsby4 = np.arange(0, bmlat.shape[0], density)
		barb_lon = bmlon[:,colsby4]
		barb_lon = barb_lon[rowsby4,:]
		barb_lat = bmlat[:,colsby4]
		barb_lat = barb_lat[rowsby4,:]
		barbu_test = U200_test[rowsby4,:]
		barbu_test = barbu_test[:,colsby4] 
		barbv_test = V200_test[rowsby4,:]
		barbv_test = barbv_test[:,colsby4] 
		m.quiver(barb_lon, barb_lat, barbu_test, barbv_test, scale = 50, scale_units = "inches", cmap=plt.cm.autumn, latlon = True)

		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_control, clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)

		barbu_control = U200_control[rowsby4,:]
		barbu_control = barbu_control[:,colsby4] 
		barbv_control = V200_control[rowsby4,:]
		barbv_control = barbv_control[:,colsby4] 
		m.quiver(barb_lon, barb_lat, barbu_control, barbv_control, scale = 50, scale_units = "inches", cmap=plt.cm.autumn, latlon = True)

		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_test-PRECT_control, clevs2, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)

		barbu_diff = barbu_test - barbu_control
		barbv_diff = barbv_test - barbv_control
		m.quiver(barb_lon, barb_lat, barbu_diff, barbv_diff, scale = 12.5, scale_units = "inches", cmap=plt.cm.autumn, latlon = True)

		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT200windbarbs_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT200windbarb_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if dRHDO_dt:
		#------------- dRHDO_dt -------------#
		print("\nPlotting " + which + " dRHDO_dt data...")
		fig = plt.figure()
		fig.suptitle(dRHDO_dt_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, dRHDO_dt_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dRHDO_dt_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, dRHDO_dt_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dRHDO_dt_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, dRHDO_dt_test-dRHDO_dt_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + dRHDO_dt_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		

		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "dRHDO_dt_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "dRHDO_dt_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if PRECT:
		#------------- PRECT -------------#
		print("\nPlotting " + which + " PRECT data...")
		fig = plt.figure()
		fig.suptitle(PRECT_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
        clevs = np.linspace(0, 15, 19)
        clevs2 = np.linspace(-4,4,17)
        if which == "JJA":
            if region == "NA":
                clevs = np.linspace(0, 15, 19)
            if region == "IM":
                clevs = np.linspace(0, 30, 16)
                clevs2 = np.linspace(-8,8,17)
            if region == "MC":
                clevs2 = clevs2 = np.linspace(-6,6,13)
            if region == "EP":
                clevs = np.linspace(0,25,16)
                clevs2 = np.linspace(-6,6,13)
        if which == "DJF":
            if region == "NA":
                clevs = np.linspace(0, 15, 19)
            if region == "IM":
                clevs = np.linspace(0, 30, 16)
                clevs2 = np.linspace(-8,8,17)
            if region == "MC":
                clevs2 = clevs2 = np.linspace(-6,6,13)
            if region == "EP":
                clevs = np.linspace(0,25,16)
                clevs2 = np.linspace(-6,6,13)
		cs = m.contourf(bmlon, bmlat, PRECT_test, clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_control, clevs, shading = 'flat', latlon = True, cmap=plt.cm.Reds)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_test-PRECT_control, clevs2, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECT_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECT_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if PRECC:
		#------------- PRECC -------------#
		print("\nPlotting " + which + " PRECC data...")
		fig = plt.figure()
		fig.suptitle(PRECC_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECC_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECC_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECC_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECC_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECC_test-PRECC_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECC_units, fontsize = 8)
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

	if PRECL:
		#------------- PRECL -------------#
		print("\nPlotting " + which + " PRECL data...")
		fig = plt.figure()
		fig.suptitle(PRECL_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECL_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECL_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECL_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECL_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECL_test-PRECL_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECL_units, fontsize = 8)
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

	if PRECT_dD:
		#------------- PRECT_dD -------------#
		print("\nPlotting " + which + " PRECT_dD data...")
		fig = plt.figure()
		fig.suptitle(PRECT_dD_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_dD_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_dD_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_dD_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_dD_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_dD_test-PRECT_dD_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECT_dD_units, fontsize = 8)
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


	if PRECT_d18O:
		#------------- PRECT_d18O -------------#
		print("\nPlotting " + which + " PRECT_d18O data...")
		fig = plt.figure()
		fig.suptitle(PRECT_d18O_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_d18O_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_d18O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_d18O_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_d18O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		if (region == 'NA') | (region == 'IM'):
			if which == 'JJA':
				clevs = np.linspace(-5,5,21)
				cs = m.contourf(bmlon, bmlat, PRECT_d18O_test-PRECT_d18O_control, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
			else:
				cs = m.contourf(bmlon, bmlat, PRECT_d18O_test-PRECT_d18O_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		else:
			cs = m.contourf(bmlon, bmlat, PRECT_d18O_test-PRECT_d18O_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECT_d18O_units, fontsize = 8)
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

	if PRECT_dxs:
		#------------- PRECT_dxs -------------#
		print("\nPlotting " + which + " PRECT_dxs data...")
		fig = plt.figure()
		fig.suptitle(PRECT_dxs_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(0,35,19)
		cs = m.contourf(bmlon, bmlat, PRECT_dxs_test, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_dxs_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(0,35,19)
		cs = m.contourf(bmlon, bmlat, PRECT_dxs_control, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_dxs_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(-10, 10, 19)
		cs = m.contourf(bmlon, bmlat, PRECT_dxs_test-PRECT_dxs_control, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECT_dxs_units, fontsize = 8)
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
		
	if PRECT_H2O:
		#------------- PRECT_H2O -------------#
		print("\nPlotting " + which + " PRECT_H2O data...")
		fig = plt.figure()
		fig.suptitle(PRECT_H2O_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_H2O_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_H2O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_H2O_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECT_H2O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECT_H2O_test-PRECT_H2O_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECT_H2O_units, fontsize = 8)
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

	if PRECST:
		#------------- PRECST -------------#
		print("\nPlotting " + which + " PRECST data...")
		fig = plt.figure()
		fig.suptitle(PRECST_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECST_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECST_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECST_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PRECST_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECST_test-PRECST_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PRECST_units, fontsize = 8)
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

	if PRECCfrac:
		#------------- PRECC / PRECT -------------#
		print("\nPlotting " + which + " PRECC / PRECT data...")
		fig = plt.figure()
		fig.suptitle(PRECC_long_name + " / \n" + PRECT_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECC_test / PRECT_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("percent", fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECC_control / PRECT_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("percent", fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, (PRECC_test/PRECT_test) - (PRECC_control/PRECT_control), 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in percent", fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECCfrac_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECCfrac_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if PRECLfrac:
		#------------- PRECL / PRECT -------------#
		print("\nPlotting " + which + " PRECL / PRECT data...")
		fig = plt.figure()
		fig.suptitle(PRECL_long_name + " / \n" + PRECT_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECL_test / PRECT_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("percent", fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PRECL_control / PRECT_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("percent", fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, (PRECL_test/PRECT_test) - (PRECL_control/PRECT_control), 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in percent", fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PRECLfrac_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PRECLfrac_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if Qlowest:
		#------------- Q at lowest level -------------#
		print("\nPlotting " + which + " Qlowest data...")
		fig = plt.figure()
		fig.suptitle(Q_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Q_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(Q_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Q_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(Q_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Q_test[-1,:,:] - Q_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + Q_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Qlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Qlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if dDVlowest:
		#------------- dDV at lowest level -------------#
		print("\nPlotting " + which + " dDVlowest data...")
		fig = plt.figure()
		fig.suptitle(dDV_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, dDV_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dDV_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, dDV_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dDV_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, dDV_test[-1,:,:] - dDV_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + dDV_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "dDVlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "dDVlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if dDVmiddle:
		#------------- dDV at middle level -------------#
		print("\nPlotting " + which + " dDVmiddle data...")
		fig = plt.figure()
		fig.suptitle(dDV_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, dDV_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dDV_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, dDV_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dDV_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, dDV_test[15,:,:] - dDV_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + dDV_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "dDVmiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "dDVmiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if d18OVlowest:
		#------------- d18OV at lowest level -------------#
		print("\nPlotting " + which + " d18OVlowest data...")
		fig = plt.figure()
		fig.suptitle(d18OV_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18OV_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18OV_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18OV_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18OV_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18OV_test[-1,:,:] - d18OV_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + d18OV_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "d18OVlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "d18OVlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	 
	if d18OVmiddle:
		#------------- d18OV at middle level -------------#
		print("\nPlotting " + which + " d18OVmiddle data...")
		fig = plt.figure()
		fig.suptitle(d18OV_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18OV_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18OV_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18OV_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(d18OV_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, d18OV_test[15,:,:] - d18OV_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + d18OV_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "d18OVmiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "d18OVmiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if dxsVlowest:
		#------------- dxsV at lowest level -------------#
		print("\nPlotting " + which + " dxsVlowest data...")
		fig = plt.figure()
		fig.suptitle(dxsV_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(0,35,19)
		cs = m.contourf(bmlon, bmlat, dxsV_test[-1,:,:], clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dxsV_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(0,35,19)
		cs = m.contourf(bmlon, bmlat, dxsV_control[-1,:,:], clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dxsV_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(-10,10,19)
		cs = m.contourf(bmlon, bmlat, dxsV_test[-1,:,:] - dxsV_control[-1,:,:], clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + dxsV_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "dxsVlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "dxsVlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if dxsVmiddle:
		#------------- dxsV at middle level -------------#
		print("\nPlotting " + which + " dxsVmiddle data...")
		fig = plt.figure()
		fig.suptitle(dxsV_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(0,35,19)
		cs = m.contourf(bmlon, bmlat, dxsV_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dxsV_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(0,35,19)
		cs = m.contourf(bmlon, bmlat, dxsV_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(dxsV_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		clevs = np.linspace(-10,10,19)
		cs = m.contourf(bmlon, bmlat, dxsV_test[15,:,:] - dxsV_control[15,:,:], clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + dxsV_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "dxsVmiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "dxsVmiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if QFLX:
		#------------- QFLX -------------#
		print("\nPlotting " + which + " QFLX data...")
		fig = plt.figure()
		fig.suptitle(QFLX_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, QFLX_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(QFLX_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, QFLX_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(QFLX_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, QFLX_test - QFLX_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + QFLX_units, fontsize = 8)
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

 	if QFLX_HDO:
	 	#------------- QFLX_HDO -------------#
		print("\nPlotting " + which + " QFLX_HDO data...")
		fig = plt.figure()
		fig.suptitle(QFLX_HDO_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, QFLX_HDO_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(QFLX_HDO_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, QFLX_HDO_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(QFLX_HDO_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, QFLX_HDO_test - QFLX_HDO_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + QFLX_HDO_units, fontsize = 8)
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

	if QFLX_H218O:
		#------------- QFLX_H218O -------------#
		print("\nPlotting " + which + " QFLX_H218O data...")
		fig = plt.figure()
		fig.suptitle(QFLX_H218O_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, QFLX_H218O_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(QFLX_H218O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, QFLX_H218O_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(QFLX_H218O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, QFLX_H218O_test - QFLX_H218O_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + QFLX_H218O_units, fontsize = 8)
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

	if VQ_H2Olowest:
		#------------- VQ_H2O lowest-------------#
		print("\nPlotting " + which + " VQ_H2Olowest data...")
		fig = plt.figure()
		fig.suptitle(VQ_H2O_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_H2O_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_H2O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_H2O_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_H2O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_H2O_test[-1,:,:] - VQ_H2O_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + VQ_H2O_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQ_H2Olowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQ_H2Olowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
		
	if VQ_H2Omiddle:	
		#------------- VQ_H2O middle-------------#
		print("\nPlotting " + which + " VQ_H2Omiddle data...")
		fig = plt.figure()
		fig.suptitle(VQ_H2O_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_H2O_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_H2O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_H2O_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_H2O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_H2O_test[15,:,:] - VQ_H2O_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + VQ_H2O_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQ_H2Omiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQ_H2Omiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if VQ_dDlowest:
		#------------- VQ_dD lowest-------------#
		print("\nPlotting " + which + " VQ_dDlowest data...")
		fig = plt.figure()
		fig.suptitle(VQ_dD_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_dD_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_dD_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_dD_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_dD_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_dD_test[-1,:,:] - VQ_dD_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + VQ_dD_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQ_dDlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQ_dDlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if VQ_dDmiddle:
		#------------- VQ_dD middle-------------#
		print("\nPlotting " + which + " VQ_dDmiddle data...")
		fig = plt.figure()
		fig.suptitle(VQ_dD_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_dD_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_dD_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_dD_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_dD_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_dD_test[15,:,:] - VQ_dD_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + VQ_dD_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQ_dDmiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQ_dDmiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if VQ_d18Olowest:

		#------------- VQ_d18O lowest-------------#
		print("\nPlotting " + which + " VQ_d18Olowest data...")
		fig = plt.figure()
		fig.suptitle(VQ_d18O_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_d18O_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_d18O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_d18O_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_d18O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_d18O_test[-1,:,:] - VQ_d18O_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + VQ_d18O_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQ_d18Olowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQ_d18Olowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if VQ_d18Omiddle:
		#------------- VQ_d18O middle-------------#
		print("\nPlotting " + which + " VQ_d18Omiddle data...")
		fig = plt.figure()
		fig.suptitle(VQ_d18O_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_d18O_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_d18O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_d18O_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VQ_d18O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VQ_d18O_test[15,:,:] - VQ_d18O_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + VQ_d18O_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VQ_d18Omiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VQ_d18Omiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if UQ_H2Olowest:
		#------------- UQ_H2O lowest-------------#
		print("\nPlotting " + which + " UQ_H2Olowest data...")
		fig = plt.figure()
		fig.suptitle(UQ_H2O_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_H2O_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_H2O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_H2O_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_H2O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_H2O_test[-1,:,:] - UQ_H2O_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + UQ_H2O_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQ_H2Olowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQ_H2Olowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if UQ_H2Omiddle:
		#------------- UQ_H2O middle-------------#
		print("\nPlotting " + which + " UQ_H2Omiddle data...")
		fig = plt.figure()
		fig.suptitle(UQ_H2O_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_H2O_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_H2O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_H2O_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_H2O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_H2O_test[15,:,:] - UQ_H2O_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + UQ_H2O_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQ_H2Omiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQ_H2Omiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	  
	if UQ_dDlowest:

		#------------- UQ_dD lowest-------------#
		print("\nPlotting " + which + " UQ_dDlowest data...")
		fig = plt.figure()
		fig.suptitle(UQ_dD_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_dD_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_dD_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_dD_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_dD_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_dD_test[-1,:,:] - UQ_dD_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + UQ_dD_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQ_dDlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQ_dDlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if UQ_dDmiddle:
		#------------- UQ_dD middle-------------#
		print("\nPlotting " + which + " UQ_dDmiddle data...")
		fig = plt.figure()
		fig.suptitle(UQ_dD_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_dD_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_dD_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_dD_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_dD_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_dD_test[15,:,:] - UQ_dD_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + UQ_dD_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQ_dDmiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQ_dDmiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if UQ_d18Olowest:

		#------------- UQ_d18O lowest-------------#
		print("\nPlotting " + which + " UQ_d18Olowest data...")
		fig = plt.figure()
		fig.suptitle(UQ_d18O_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_d18O_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_d18O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_d18O_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_d18O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_d18O_test[-1,:,:] - UQ_d18O_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + UQ_d18O_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQ_d18Olowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQ_d18Olowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if UQ_d18Omiddle:
		#------------- UQ_d18O middle-------------#
		print("\nPlotting " + which + " UQ_d18Omiddle data...")
		fig = plt.figure()
		fig.suptitle(UQ_d18O_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_d18O_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_d18O_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_d18O_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UQ_d18O_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UQ_d18O_test[15,:,:] - UQ_d18O_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + UQ_d18O_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UQ_d18Omiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UQ_d18Omiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	  
	if VTlowest:
		#------------- VT lowest-------------#
		print("\nPlotting " + which + " VTlowest data...")
		fig = plt.figure()
		fig.suptitle(VT_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VT_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VT_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VT_test[-1,:,:] - VT_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + VT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VTlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VTlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if VTmiddle:
		#------------- VT middle-------------#
		print("\nPlotting " + which + " VTmiddle data...")
		fig = plt.figure()
		fig.suptitle(VT_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VT_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VT_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(VT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, VT_test[15,:,:] - VT_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + VT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "VTmiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "VTmiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if UTlowest:

		#------------- UT lowest-------------#
		print("\nPlotting " + which + " UTlowest data...")
		fig = plt.figure()
		fig.suptitle(UT_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UT_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UT_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UT_test[-1,:,:] - UT_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + UT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UTlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UTlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if UTmiddle:
		#------------- UT middle-------------#
		print("\nPlotting " + which + " UTmiddle data...")
		fig = plt.figure()
		fig.suptitle(UT_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UT_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UT_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(UT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, UT_test[15,:,:] - UT_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + UT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "UTmiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "UTmiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if TTEND_TOTlowest:

		#------------- TTEND_TOT lowest-------------#
		print("\nPlotting " + which + " TTEND_TOTlowest data...")
		fig = plt.figure()
		fig.suptitle(TTEND_TOT_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, TTEND_TOT_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(TTEND_TOT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, TTEND_TOT_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(TTEND_TOT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, TTEND_TOT_test[-1,:,:] - TTEND_TOT_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + TTEND_TOT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "TTEND_TOTlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "TTEND_TOTlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if TTEND_TOTmiddle:
		#------------- TTEND_TOT middle-------------#
		print("\nPlotting " + which + " TTEND_TOTmiddle data...")
		fig = plt.figure()
		fig.suptitle(TTEND_TOT_long_name + ", middle model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, TTEND_TOT_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(TTEND_TOT_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, TTEND_TOT_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(TTEND_TOT_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, TTEND_TOT_test[15,:,:] - TTEND_TOT_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + TTEND_TOT_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "TTEND_TOTmiddle_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "TTEND_TOTmiddle_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if Ulowest:
		#------------- U lowest-------------#
		print("\nPlotting " + which + " U lowest data...")
		fig = plt.figure()
		fig.suptitle(U_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(U_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(U_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U_test[-1,:,:] - U_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + U_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Ulowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Ulowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if Vlowest:
		#------------- V lowest-------------#
		print("\nPlotting " + which + " V lowest data...")
		fig = plt.figure()
		fig.suptitle(V_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(V_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(V_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V_test[-1,:,:] - V_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + V_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Vlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Vlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#


	if Tlowest:
		#------------- T lowest-------------#
		print("\nPlotting " + which + " Tlowest data...")
		fig = plt.figure()
		fig.suptitle(T_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(T_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(T_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T_test[-1,:,:] - T_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + T_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "Tlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "Tlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	
	if T850:
		#------------- T850 -------------#
		print("\nPlotting " + which + " T850 data...")
		fig = plt.figure()
		fig.suptitle(T850_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T850_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(T850_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(T850_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T850_test[:,:] - T850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + T850_units, fontsize = 8)
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

	if T500:
		#------------- T500 -------------#
		print("\nPlotting " + which + " T500 data...")
		fig = plt.figure()
		fig.suptitle(T500_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T500_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(T500_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(T500_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T500_test[:,:] - T500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + T500_units, fontsize = 8)
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

	if T200:
		#------------- T200 -------------#
		print("\nPlotting " + which + " T200 data...")
		fig = plt.figure()
		fig.suptitle(T200_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T200_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(T200_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(T200_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, T200_test[:,:] - T200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + T200_units, fontsize = 8)
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
	
	if OMEGAlowest:
		#------------- OMEGA lowest-------------#
		print("\nPlotting " + which + " OMEGAlowest data...")
		fig = plt.figure()
		fig.suptitle(OMEGA_long_name + ", lowest model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA_test[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(OMEGA_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(OMEGA_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA_test[15,:,:] - OMEGA_control[15,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + OMEGA_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "OMEGAlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "OMEGAlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if CLDHGH:
		#------------- CLDHGH -------------#
		print("\nPlotting " + which + " CLDHGH data...")
		fig = plt.figure()
		fig.suptitle(CLDHGH_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, CLDHGH_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(CLDHGH_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, CLDHGH_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(CLDHGH_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, CLDHGH_test - CLDHGH_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + CLDHGH_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "CLDHGH_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "CLDHGH_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#
	if PS:
		#------------- PS -------------#
		print("\nPlotting " + which + " PS data...")
		fig = plt.figure()
		fig.suptitle(PS_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PS_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PS_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PS_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PS_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PS_test-PS_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PS_units, fontsize = 8)
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
	if PSL:

		#------------- PSL -------------#
		print("\nPlotting " + which + " PSL data...")
		fig = plt.figure()
		fig.suptitle(PSL_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PSL_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PSL_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PSL_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PSL_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PSL_test-PSL_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PSL_units, fontsize = 8)
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
	if PSI500:

		#------------- PSI500 -------------#
		print("\nPlotting " + which + " PSI500 data...")
		fig = plt.figure()
		fig.suptitle(PSI500_long_name, fontweight = 'bold', fontsize = 14)
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PSI500_test, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PSI500_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PSI500_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(PSI500_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, PSI500_test-PSI500_control, 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + PSI500_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		
		
		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "PSI500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "PSI500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if RHlower:
		#------------- RH lower-------------#
		print("\nPlotting " + which + " RHlowest data...")
		fig = plt.figure()
		fig.suptitle(RH_long_name + ", lower model level", fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH_test[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(RH_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(RH_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH_test[-1,:,:] - RH_control[-1,:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + RH_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		

		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "RHlowest_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "RHlowest_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if RH850:
		#------------- RH850 -------------#
		print("\nPlotting " + which + " RH850 data...")
		fig = plt.figure()
		fig.suptitle(RH850_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH850_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(RH850_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(RH850_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH850_test[:,:] - RH850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + RH850_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		

		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "RH850_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "RH850_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if RH500:
		#------------- RH500 -------------#
		print("\nPlotting " + which + " RH500 data...")
		fig = plt.figure()
		fig.suptitle(RH500_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH500_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(RH500_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(RH500_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH500_test[:,:] - RH500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + RH500_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		

		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "RH500_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "RH500_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if RH200:
		#------------- RH200 -------------#
		print("\nPlotting " + which + " RH200 data...")
		fig = plt.figure()
		fig.suptitle(RH200_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH200_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(RH200_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(RH200_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, RH200_test[:,:] - RH200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + RH200_units, fontsize = 8)
		plt.title("test - control difference", fontsize = 8)		

		# Make things pretty
		plt.subplots_adjust(hspace = .15)
		plt.savefig(figdir + "/"  + which + "/" + "RH200_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
		print("Created " + which + "/" + "RH200_" + testfn + "-" + cntrlfn + "." + ftype + '\n')
		if show_figure:
			plt.show()
		plt.clf()
		plt.cla()
		plt.close(fig)
		#---------------------------------#

	if U850:
		#------------- U850 -------------#
		print("\nPlotting " + which + " U850 data...")
		fig = plt.figure()
		fig.suptitle(U850_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U850_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(U850_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(U850_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U850_test[:,:] - U850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + U850_units, fontsize = 8)
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

	if U500:
		#------------- U500 -------------#
		print("\nPlotting " + which + " U500 data...")
		fig = plt.figure()
		fig.suptitle(U500_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U500_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(U500_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(U500_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U500_test[:,:] - U500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + U500_units, fontsize = 8)
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

	if U200:
		#------------- U200 -------------#
		print("\nPlotting " + which + " U200 data...")
		fig = plt.figure()
		fig.suptitle(U200_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U200_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(U200_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(U200_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, U200_test[:,:] - U200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + U200_units, fontsize = 8)
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

	if V850:
		#------------- V850 -------------#
		print("\nPlotting " + which + " V850 data...")
		fig = plt.figure()
		fig.suptitle(V850_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V850_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(V850_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(V850_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V850_test[:,:] - V850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + V850_units, fontsize = 8)
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

	if V500:
		#------------- V500 -------------#
		print("\nPlotting " + which + " V500 data...")
		fig = plt.figure()
		fig.suptitle(V500_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V500_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(V500_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(V500_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V500_test[:,:] - V500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + V500_units, fontsize = 8)
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

	if V200:
		#------------- V200 -------------#
		print("\nPlotting " + which + " V200 data...")
		fig = plt.figure()
		fig.suptitle(V200_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V200_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(V200_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(V200_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, V200_test[:,:] - V200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + V200_units, fontsize = 8)
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

	if OMEGA850:
		#------------- OMEGA850 -------------#
		print("\nPlotting " + which + " OMEGA850 data...")
		fig = plt.figure()
		fig.suptitle(OMEGA850_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA850_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(OMEGA850_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(OMEGA850_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA850_test[:,:] - OMEGA850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + OMEGA850_units, fontsize = 8)
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

	if OMEGA500:
		#------------- OMEGA500 -------------#
		print("\nPlotting " + which + " OMEGA500 data...")
		fig = plt.figure()
		fig.suptitle(OMEGA500_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA500_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(OMEGA500_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(OMEGA500_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA500_test[:,:] - OMEGA500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + OMEGA500_units, fontsize = 8)
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

	if OMEGA200:
		#------------- OMEGA200 -------------#
		print("\nPlotting " + which + " OMEGA200 data...")
		fig = plt.figure()
		fig.suptitle(OMEGA200_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA200_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(OMEGA200_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(OMEGA200_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, OMEGA200_test[:,:] - OMEGA200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + OMEGA200_units, fontsize = 8)
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


	if Z3850:
		#------------- Z3850 -------------#
		print("\nPlotting " + which + " Z3850 data...")
		fig = plt.figure()
		fig.suptitle(Z3850_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Z3850_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(Z3850_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Z3850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(Z3850_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Z3850_test[:,:] - Z3850_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + Z3850_units, fontsize = 8)
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

	if Z3500:
		#------------- Z3500 -------------#
		print("\nPlotting " + which + " Z3500 data...")
		fig = plt.figure()
		fig.suptitle(Z3500_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Z3500_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(Z3500_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Z3500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(Z3500_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Z3500_test[:,:] - Z3500_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + Z3500_units, fontsize = 8)
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

	if Z3200:
		#------------- Z3200 -------------#
		print("\nPlotting " + which + " Z3200 data...")
		fig = plt.figure()
		fig.suptitle(Z3200_long_name, fontweight = 'bold')
		# test data
		plt.subplot(3,1,1)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Z3200_test[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(Z3200_units, fontsize = 8)
		test_title = user_test_plot_title + ", " + which
		if user_test_plot_title is None:
			test_title = "test: " + testfn
		plt.title(test_title, fontsize = 8)
		# control data
		plt.subplot(3,1,2)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Z3200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(Z3200_units, fontsize = 8)
		control_title = user_control_plot_title + ", " + which
		if user_control_plot_title is None:
			control_title = "control: " + cntrlfn
		plt.title(control_title, fontsize = 8)
		# difference data 
		plt.subplot(3,1,3)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, Z3200_test[:,:] - Z3200_control[:,:], 18, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label("difference in " + Z3200_units, fontsize = 8)
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
	



# The following writes the directory that the figures are stored in to a hidden text file on the user's home directory
# This can be used in, for example, a shell script, to read where the figures are located without having to hard code it in
# E.g.: image_dir='cat ~/.tempfile.txt' in bash will store the image directory in the env variable called $image_dir
home = os.path.expanduser("~")
figdirtemp = open(home+"/.figdirtemp.txt", "w")
figdirtemp.write(figdir)

regiontemp = open(home+"/.regiontemp.txt", "w")
regiontemp.write(region_name)


