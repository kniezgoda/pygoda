# Python code for differences plots
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
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

#------------- Set the input directory
# Kyle's personal laptop
directory = '/Users/kyleniezgoda/Library/Mobile Documents/com~apple~CloudDocs/Documents/Noone_Group/Yellowstone/model_runs/'
# Kyle's work computer
# directory = "/Users/kyleniezgoda/Yellowstone_Runs/model_runs/"
#------------- 

#------------- Set your control and test file here
cntrl = "F.C5.2deg.wiso.defaultSSTICE_kn002/F.C5.2deg.wiso.defaultSSTICE_kn002_JJA_climo.nc"
test = "F.C5.2deg.wiso.obs6kSST_kn003/F.C5.2deg.wiso.obs6kSST_kn003_JJA_climo.nc"
#------------- 

#------------- Define what to plot (a netcdf variable name, must match exactly)
# denom_var can be used to calculate ratios, e.g. PRECT_H218O/PRECT, the percent of rainfall that is 18O
var = "PRECT_HDO"
denom_var = "PRECT"
#------------- 

#------------- Are you plotting water vapor isotopes? Do you want them converted to delta values (most likely true if so!)? 
is_vapor_isotope = False
if is_vapor_isotope:
	print("\n!!!\n!!!\n!!!is_vapor_isotope IS SET TO TRUE!!!\n!!!\n!!!\n!!!")
#------------- 

#------------- Set model level (-1 = ground, -2 = next up. 0 = top of model, 1 = last level before top of model, etc...)
#------------- Only needed if the data is 3-dimensional in space. Otherwise, this variable won't be used
lev = -1
#------------- 

#------------- Show figure? Save figure?
showfig = True
savefig = False # savefig not working right now
#------------- 

#------------- If saving figure, what directory should the file be saved to? Default is current directory.
figdir = "/Users/kyleniezgoda/Library/Mobile Documents/com~apple~CloudDocs/Documents/Python/figures/"
#------------- 

#------------- File type for saved image. Don't include any periods, just the suffix
ftype = 'ps'

#------------- Set map coordinates
southern_lat = -50
northern_lat = 50
left_lon = 0
right_lon = 355
#------------- 

#------------- Add conversions here (e.g., to. convert temp data to celcius, set add = -273, multiply = 0)
# To convert from m/s to cm/hr, set add = 0, multiply = 100*60*60
add = 0
multiply = 0
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

def apply_algebra(x):
	if multiply != 0 and add == 0:
		return [x*multiply, units + " * " + str(multiply)]
	if multiply == 0 and add != 0:
		return [x + add, units + " + " + str(add)]
	if multiply != 0 and add !=0:
		return [x*multiply + add, units + " * " + str(multiply) + " + " + str(add)]
	if multiply == 0 and add == 0:
		return [x, units]

# removes the extension from the string, makes plotting and file names nicer later on
cntrlfn = os.path.splitext(os.path.split(cntrl)[1])[0]
testfn = os.path.splitext(os.path.split(test)[1])[0]

# Read the netcdf file
cntrldata = Dataset(directory+cntrl, mode='r')
testdata = Dataset(directory+test, mode='r')

# Extract lat, lon, and lev data
lats = cntrldata.variables['lat'][:]
lons = cntrldata.variables['lon'][:]
if verbose:
	print("Latitude data from netcdf:\n" + str(lats) + "\n")
	print("Longitude data from netcdf:\n" + str(lons) + "\n")
bmlon, bmlat = np.meshgrid(lons, lats)

levs = cntrldata.variables['lev'][:]

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------- Extract data -------------#

# Determine the number of dimensions
num_dims = len(cntrldata.variables[var].shape) - 1
# Extract the units - we will add the algebra conversion laters
units = str(cntrldata[var].units)

# Extract the variable data
# Take care of the possibility of 2-d and 3-d data
# Read in the denominator variable if needed
if num_dims == 2:
	cntrlvar = cntrldata.variables[var][0,:,:].squeeze()
	testvar = testdata.variables[var][0,:,:].squeeze()
	if denom_var is not None:
		denomcntrlvar = cntrldata.variables[denom_var][0,:,:].squeeze()
		denomtestvar = testdata.variables[denom_var][0,:,:].squeeze()
if num_dims == 3:
	cntrlvar = cntrldata.variables[var][0,lev,:,:].squeeze()
	testvar = testdata.variables[var][0,lev,:,:].squeeze()
	if denom_var is not None:
		denomcntrlvar = cntrldata.variables[denom_var][0,lev,:,:].squeeze()
		denomtestvar = testdata.variables[denom_var][0,lev,:,:].squeeze()

# Convert the isotope data to delta values if we are working with vapor isotopes
if is_vapor_isotope:
	cntrlH2OV = cntrldata.variables['H2OV'][0,lev,:,:].squeeze()
	testH2OV = testdata.variables['H2OV'][0,lev,:,:].squeeze()
	cntrlvar = (cntrlvar/cntrlH2OV - 1) * 1000
	testvar = (testvar/testH2OV - 1) * 1000

# Convert the data as per the user instructions
testvar = apply_algebra(testvar)[0]
cntrlvar = apply_algebra(cntrlvar)[0]

# This adds multiplication and addition strings to the units string for plotting
units = apply_algebra(testvar)[1]

# Compute differences and store in the "plot" variable
differences = testvar - cntrlvar

# Handle for the case when denom_var is set
if denom_var is not None:
	denomtestvar = apply_algebra(denomtestvar)[0]
	denomcntrlvar = apply_algebra(denomcntrlvar)[0]
	testvar = testvar/denomtestvar
	cntrlvar = cntrlvar/denomcntrlvar
	differences = testvar - cntrlvar
	units = "percent"

if verbose:
	print("Test data post algebra:\n" + str(testvar) + "\n")
	print("Control data post algebra:\n" + str(cntrlvar) + "\n")
	print("Difference data:\n" + str(differences) + '\n')

if is_vapor_isotope:
	units = "permil"


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------- Map creation -------------#

# Create the basemap
fig = plt.figure()
# Add a somewhat meaningful title
if denom_var is not None:
	fig.suptitle(str(var) + ' / ' + str(denom_var), fontweight = 'bold', fontsize = 14)
else:
	fig.suptitle(str(var), fontweight = 'bold', fontsize = 14)

if is_vapor_isotope:
	fig.suptitle(str(var) + " at level: " + str(lev))


# Plot test data
ax = fig.add_subplot(311)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, testvar, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(units)
ax.set_title(testfn, fontsize = 8)

# Plot control data 
ax = fig.add_subplot(312)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, cntrlvar, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label(units)
ax.set_title(cntrlfn, fontsize = 8)

ax = fig.add_subplot(313)
ax.set_title(testfn + " - " + cntrlfn, fontsize = 8)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')

# Describe the color levels if set by the user
# Also hack the data so we don't have grey areas in the map
if bounds is not None:
	if len(bounds) == 2:
		# Set any values in plot that are greater than bounds[1] to equal bounds[1]
		# Set and values in plot that are less than bounds[0] to equal bounds[0]
		# We do this because otherwise the plot shows a bunch of grey area for values of plot that are beyond the bounds
		# This is a bit hacky
		for i in range(differences.shape[0]):
			for j in range(differences.shape[1]):
				if differences[i,j] < bounds[0]:
					differences[i,j] = bounds[0]
				if differences[i,j] > bounds[1]:
					differences[i,j] = bounds[1]

		# Plot with clevs
		clevs = np.linspace(bounds[0], bounds[1], 21)
		cs = m.contourf(bmlon, bmlat, differences, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
else:
	# Plot without clevs
	cs = m.contourf(bmlon, bmlat, differences, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)

# Draw a color bar
cbar = m.colorbar(cs, location='right', pad="5%")
# Put units on the color bar - the algebra will be included here
cbar.set_label("difference in " + units)


# Show and save
if savefig:
	# plt.savefig(var + "_" + figdir + testfn + "-" + cntrlfn + ".pdf", bbox_inches='tight', dpi = 500)
	plt.savefig(figdir + var + "_" + testfn + "-" + cntrlfn + "." + ftype, bbox_inches='tight', dpi = 500)
if showfig:
	plt.show()
