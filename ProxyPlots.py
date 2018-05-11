# This code creates a map of MH-PI PRECT_d18O data from the F.C5 amip runs.
# Proxy data from the P2C2 year one annual review are plotted on top of the PRECT_d18O contours.
# These proxy data are from Bronwen and are stored in a file called Holocene ITCZ Records

# This code is not robust and not made for easy switching in and out of variables.
# The proxy data sheet isn't event read in...it is all hardcoded into a pandas DataFrame

import os, glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import shiftgrid 


# # Set the directory that contains the ANN climatology file
root = '/Users/kyleniezgoda/LocalDocuments/Yellowstone_Runs/model_runs/'
cntrldir = root+"F.C5.2deg.wiso.defaultSSTICE_kn002/"
testdir = root+"F.C5.2deg.wiso.obs6kSST_kn003/"

#COAS Climate machine
# root = '/home/server/student/homes/kniezgod/model_runs/' 
# cntrldir = root+"F.C5.2deg.wiso.defaultSSTICE_kn002/"
# testdir = root+"F.C5.2deg.wiso.obs6kSST_kn003/"

# Location where file will be printed to
figdir = "."

# COAS Climate machine
# figdir = "/home/server/student/homes/kniezgod/python/figures/F.C5.2deg.wiso.obs6kSST_kn003-F.C5.2deg.wiso.defaultSSTICE_kn002/"

ftype = 'ps'

#------------- Set map coordinates
# Default
southern_lat = -50
northern_lat = 50
left_lon = 0
right_lon = 355

#------------- Show figures?
# This will take a while and it will show literally every plot along the way.
show_figure = True

#------------- Custom subplot title names 
# these are the names of the subplots, not the main figure --- the main figure always has the same name (the variable long name from the netcdf file)
# Set to None (no quotes) for default subplot names, the name of the file used
user_test_plot_title = "Test: Mid-Holocene (6kya)"
user_control_plot_title = "Control: Pre-industrial (1850)"

# Set the proxy data
d18O_proxy = pd.DataFrame({'Site' : ["Qunf", "Dongge", "LiangLuar", "GunungBuda", "KNI51", "RioGrandeDeNorte", "Botuvera"], \
	'lat' : [17.07, 25.17, -8.32, 4, -15.18, -5.36, -27.13], \
	'lon' : [53.41, 108.5, 120.26, 114, 128.37, -37.44, -49.09], \
	'6kya_value' : [-0.2, -8.51, -6.2, -9.65, -6.58, -6, -2.34], \
	'1850_value' : [-0.81, -7.74, -6.07, -9.16, -7.55, -2.5, -3.35]})
d18O_proxy['6kya-1850_difference'] = d18O_proxy['6kya_value'] - d18O_proxy['1850_value']
proxy_x, proxy_y = np.meshgrid(d18O_proxy['lon'], d18O_proxy['lat'])
num_proxies = d18O_proxy.shape[0]

which = "JJA"
# Find the file
print "Looking for control files in " + cntrldir
cntrldata = glob.glob(cntrldir+"*"+which+"*.nc")
print "Looking for test files in " + testdir
testdata = glob.glob(testdir+"*"+which+"*.nc")

if len(cntrldata) > 0:
	print "\nFound control " + which + " file: " + os.path.split(cntrldata[0])[1]
	cntrldata = cntrldata[0]
else:
	print "Could not find control " + which + " file..."
	print "No " + which + " plots will be created..."
if len(testdata) > 0:
	print "Found test " + which + " file: " + os.path.split(testdata[0])[1]
	testdata = testdata[0]
else:
	print "Could not find test " + which + " file..."
	print "No " + which + " plots will be created..."

# Assign file names for plotting later
cntrlfn = os.path.splitext(os.path.split(cntrldata)[1])[0]
testfn = os.path.splitext(os.path.split(testdata)[1])[0]

# Read the data
cntrldata = Dataset(cntrldata, mode='r')
testdata = Dataset(testdata, mode='r')

# Extract lat, lon, and lev data
lats = cntrldata.variables['lat'][:]
lons = cntrldata.variables['lon'][:]
bmlon, bmlat = np.meshgrid(lons, lats)
levs = cntrldata.variables['lev'][:]

# Extract d18O precip data
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

# Make a basemap of the d18O precip data
#------------- PRECT_d18O -------------#
print("\nPlotting " + which + " PRECT_d18O data...")
fig, ax = plt.subplots()
fig.suptitle(PRECT_d18O_long_name, fontweight = 'bold', fontsize = 14)

# Plot test data
plt.subplot(3,1,1)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
min_1850 = -22
max_1850 = 0
clevs = np.linspace(min_1850,max_1850,21)
cs = m.contourf(bmlon, bmlat, PRECT_d18O_test, clevs, shading = 'flat', latlon = True, cmap=plt.cm.coolwarm)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(PRECT_d18O_units, fontsize = 8)
proxy = m.scatter(proxy_x[0,:],proxy_y[:,0], c = d18O_proxy['6kya_value'], vmin = min_1850, vmax = max_1850, s= 100, cmap=plt.cm.coolwarm, latlon=True)
test_title = user_test_plot_title + ", " + which
if user_test_plot_title is None:
	test_title = "test: " + testfn
plt.title(test_title, fontsize = 8)

# Plot control data
plt.subplot(3,1,2)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary()
# m.fillcontinents(color='#cc9966',lake_color='#99ffff')
min_1850 = -22
max_1850 = 0
clevs = np.linspace(min_1850,max_1850,21)
cs = m.contourf(bmlon, bmlat, PRECT_d18O_control, clevs, shading = 'flat', latlon = True, cmap=plt.cm.coolwarm)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(PRECT_d18O_units, fontsize = 8)
proxy = m.scatter(proxy_x[0,:],proxy_y[:,0], c = d18O_proxy['1850_value'], vmin = min_1850, vmax = max_1850, s= 100, cmap=plt.cm.RdBu_r, latlon=True)
control_title = user_control_plot_title + ", " + which
if user_control_plot_title is None:
	control_title = "control: " + cntrlfn
plt.title(control_title, fontsize = 8)

# Plot difference data 
plt.subplot(3,1,3)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
min_diff = -4
max_diff = 4
clevs = np.linspace(min_diff,max_diff,21)
cs = m.contourf(bmlon, bmlat, PRECT_d18O_test-PRECT_d18O_control, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label("difference in " + PRECT_d18O_units, fontsize = 8)
proxy = m.scatter(proxy_x[0,:],proxy_y[:,0], c = d18O_proxy['6kya-1850_difference'], vmin = min_diff, vmax = max_diff, s= 100, cmap=plt.cm.RdBu_r, latlon=True)
plt.title("test - control difference", fontsize = 8)		


# Make things pretty
plt.subplots_adjust(hspace = .1)
plt.savefig(figdir + "/"  + "ProxyVsModel_PRECT_d18O" + "." + ftype, bbox_inches='tight', dpi = 500)
print("Created ProxyVsModel_PRECT_d18O" + "." + ftype + '\n')
if show_figure:
	plt.show()
plt.clf()
plt.cla()
plt.close(fig)
#---------------------------------#


