#!/Users/kyleniezgoda/anaconda/bin/python
# Python code for differences plots
import os, glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import shiftgrid 
import sys
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

'''
Flags:
	-v: netcdf variable to plot
	-test: test netcdf file
	-control: control netcdf file
	-sl_test (optional): shift the longitude of the test file (deg)
	-sl_control (optional): shift the longitude of the control file
'''
######################################
####### User-entered data here #######
######################################
# Read in flags and arguments associated with each flag
flags = []
args = []
for n, arg in enumerate(sys.argv[1:]):
	if n % 2 == 0:
		flags.append(arg)
		args.append(sys.argv[n+2])
	else:
		continue

# Make sure flags exist that must exist
if ("-v" not in flags) | ("-test" not in flags) | ("-control" not in flags):
	print "A needed flag is missing. Required flags are:\n\t-v (variable to plot)\n\t-test (test file name)\n\t-control (control file name)"

# Assign the necessary args to variables for the code
var_ind = [n for n, f in enumerate(flags) if f == "-v"]
var = args[var_ind[0]]
test_file_ind = [n for n, f in enumerate(flags) if f == "-test"]
test_file = args[test_file_ind[0]]
control_file_ind = [n for n, f in enumerate(flags) if f == "-control"]
control_file = args[control_file_ind[0]]

#------------- Do you want to see a bunch of text?
verbose = False

#------------- Set the input files
#test_file = "ONLY_SST.MH6KAsstice.clim.obs.diddled.nc"
#control_file = "MH6KAsstice.clim.fromCCSM4b40.diddled.nc"

#------------- Define what to plot (a netcdf variable name, must match exactly)
# denom_var can be used to calculate ratios, e.g. PRECT_H218O/PRECT, the percent of rainfall that is 18O
#var = "SST_cpl"
#------------- 

#------------- Show figure? Save figure?
showfig = True
savefig = True 
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

# Shift the data if the flag exists
if ("-sl_test" in flags):
	sltest_ind = [n for n, f in enumerate(flags) if f == "-sl_test"]
	sltestdeg = int(args[sltest_ind[0]])
	print "Shifting the test longitude to make " + str(sltestdeg) + " the left-most edge of the figure."
	right = testvar[:,lons >= sltestdeg]
	testvar = np.hstack((testvar[:,lons >= sltestdeg], testvar[:,lons < sltestdeg]))


# Average over time
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
plt.subplot(311)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, testvar, 20, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(units)
plt.title('Test: ' + testfn)

# Plot control data 
plt.subplot(312)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, cntrlvar, 20,  shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(units)
plt.title("Control: " + cntrlfn)

plt.subplot(313)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, diff, 20, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(units)
plt.title("Test - control  difference")

# Show and save
if savefig:
	plt.savefig(figdir + var + "_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
	print "Saved figure " + figdir + var + "_" + testfn + "-" + cntrlfn + "." + ftype 
if showfig:
	plt.show()