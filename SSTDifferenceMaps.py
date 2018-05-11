# Python code for differences plots
import os, glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import shiftgrid 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

######################################
####### User-entered data here #######
######################################

#------------- Do you want to see a bunch of text?
verbose = False

#------------- Set the input files
test_file = ""
control_file = ""

#------------- Define what to plot (a netcdf variable name, must match exactly)
# denom_var can be used to calculate ratios, e.g. PRECT_H218O/PRECT, the percent of rainfall that is 18O
var = "SST_cpl"
#------------- 

#------------- Show figure? Save figure?
showfig = True
savefig = False # savefig not working right now
#------------- 

#------------- If saving figure, what directory should the file be saved to? Default is current directory.
figdir = "./"
#------------- 

#------------- File type for saved image. Don't include any periods, just the suffix
ftype = 'ps'

#------------- Set map coordinates
southern_lat = -50
northern_lat = 50
left_lon = 0
right_lon = 355
#------------- 

#------------- Set the bounds for the color bar. A good idea is to set it as None first and then adjust the bounds as necessary
# bounds = [-.03, .03]
bounds = None
#------------- 


##############################################
###### Don't change anything below here ######
##############################################


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------- Beginning of code -------------#
# Assign file names for plotting later
cntrlfn = os.path.splitext(os.path.split(control_file)[1])[0]
testfn = os.path.splitext(os.path.split(test_file)[1])[0]

# Read the data
cntrldata = Dataset(control_file, mode='r')
testdata = Dataset(test_file, mode='r')

# Extract lat, lon, and lev data
lats = cntrldata.variables['lat'][:]
lons = cntrldata.variables['lon'][:]
bmlon, bmlat = np.meshgrid(lons, lats)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------- Extract data -------------#

testvar = testdata.variables[var][:].squeeze()
cntrlvar = cntrldata.variables[var][:].squeeze()
diff = testvar - cntrlvar

units = cntrldata.variables[var].units
long_name = cntrldata.variables[var].long_name


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------- Map creation -------------#
fig = plt.figure()
fig.suptitle(long_name)
# Plot test data
plt.subplot(1,3,1)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, testvar, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(units)
plt.title("Test: Mid-Holocene (6 kya)")

# Plot control data 
plt.subplot(1,3,2)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, cntrlvar, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(units)
plt.title("Control: Pre-industrial (1850)")

plt.subplot(1,3,3)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, plot, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(units)
plt.title("Test - control difference")

# Show and save
if savefig:
	plt.savefig(figdir + var + "_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
if showfig:
	plt.show()