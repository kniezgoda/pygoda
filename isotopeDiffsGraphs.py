# Python code for isotope differences plots
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

# Set the input directory
directory = "/Users/kyleniezgoda/Yellowstone_Runs/model_runs/"

# Set your control and test file here
cntrl = "F2deg.wiso.defaultSSTICE_kn001/F2deg.wiso.defaultSSTICE_kn001_DJF_climo.nc"
test = "isoCAM5CLM4_2deg_wdiddledobs6kSSTs/isoCAM5CLM4_2deg_wdiddledobs6kSSTs_DJF_climo.nc"

# Set level ("ground", "850hPa", "500hPa", 200hPa)
# These actually aren't the true pressures since CAM is terrain following.
# Ground should be the ground, pressure will be different for each location based on topology.
# 200 is safe, B gets really small that high up
# 850 and 500 are probably more effected by this, but not too much. 
level = "ground"

# Define what to plot ("d18OV", "dDV", "temp")
plotwhat = "d18OV"

showfig = True
savefig = False

##############################################
###### Don't change anything below here ######
##############################################


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-------------Beginning of code -------------#

# removes the extension from the string, makes plotting and file names nicer later on
cntrlfn = os.path.splitext(cntrl)[0]
testfn = os.path.splitext(test)[0]

# Read the netcdf file
cntrldata = Dataset(directory+cntrl, mode='r')
testdata = Dataset(directory+test, mode='r')

# Extract lat, lon, and lev data
lats = cntrldata.variables['lat'][:]
lons = cntrldata.variables['lon'][:]
bmlon, bmlat = np.meshgrid(lons, lats)

levs = cntrldata.variables['lev'][:]
lev850 = [n for n,i in enumerate([abs(l - 850) for l in levs]) if i == min([abs(l - 850) for l in levs])]
lev500 = [n for n,i in enumerate([abs(l - 500) for l in levs]) if i == min([abs(l - 500) for l in levs])]
levground = -1

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# User set level
if level == "ground":
	lev = levground
elif level == "850hPa":
	lev = lev850
elif level == "500hPa":
	lev = lev500
elif level == "200hPa":
	lev = lev200
else:
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "LEVEL NOT SET"
	print "Ground level being used, title of plot may be misleading!"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	lev = levground

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Extract data. Variable names are typical CAM5 output with isotopes enabled
cntrltemp = cntrldata.variables['T'][0,lev,:,:].squeeze()
cntrlHDOV = cntrldata.variables['HDOV'][0,lev,:,:].squeeze()
cntrlH218OV = cntrldata.variables['H218OV'][0,lev,:,:].squeeze()
cntrlH2OV = cntrldata.variables['H2OV'][0,lev,:,:].squeeze()
cntrl_PRECRC_H218Or = cntrldata.variables['PRECRC_H218Or'][0,:,:] * 100 * 60 * 60 * 24 # convert to cm / day


testtemp = testdata.variables['T'][0,lev,:,:].squeeze()
testHDOV = testdata.variables['HDOV'][0,lev,:,:].squeeze()
testH218OV = testdata.variables['H218OV'][0,lev,:,:].squeeze()
testH2OV = testdata.variables['H2OV'][0,lev,:,:].squeeze()
test_PRECRC_H218Or = testdata.variables['PRECRC_H218Or'][0,:,:] * 100 * 60 * 60 * 24 # convert to cm / day


# Compute d18O and dD in delta values
cntrldDV = (cntrlHDOV/cntrlH2OV - 1) * 1000
cntrld18OV = (cntrlH218OV/cntrlH2OV - 1) * 1000

# Compute d18O and dD in delta values
testdDV = (testHDOV/testH2OV - 1) * 1000
testd18OV = (testH218OV/testH2OV - 1) * 1000

# Convert temp data to celcius
cntrltemp = cntrltemp - 273
testtemp = testtemp - 273

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Subtract test - control for all the variables

d18OV = testd18OV - cntrld18OV
dDV = testdDV - cntrldDV
temp = testtemp - cntrltemp
PRECRC_H218Or = test_PRECRC_H218Or - cntrl_PRECRC_H218Or
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Creation of maps

# Regional box edge definitions
# box = pd.read_excel("/Users/kyleniezgoda/Documents/Noone_Group/Protocols/RegionalBoxes.xlsx")

# Set variable to plot here
if plotwhat == "d18OV":
	plot = d18OV
elif plotwhat == "dDV":
	plot = dDV
elif plotwhat == "temp":
	plot = temp
else:
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "PLOTWHAT NOT SET"
	print "d18OV will be plotted, title of plot may be misleading!"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	plot = d18OV


# # Plot precrc_h218or
# m = bm(projection = 'cea', llcrnrlat=-40,urcrnrlat=40, llcrnrlon=-55,urcrnrlon=300,resolution='c')
# m.drawcoastlines()
# m.drawmapboundary(fill_color='0.3')
# clevs = np.linspace(-.5, .5, 21)
# cs = m.contourf(bmlon, bmlat, PRECRC_H218Or, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
# cbar = m.colorbar(cs, location='bottom', pad="5%")
# plt.show()


# Plot plotwhat code:
m = bm(projection = 'cea', llcrnrlat=-40,urcrnrlat=40, llcrnrlon=-55,urcrnrlon=300,resolution='c')
m.drawcoastlines()
# m.drawparallels(np.arange(-90.,120.,15.))
# m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='0.3')
clevs = np.linspace(-3, 3, 21)
cs = m.contourf(bmlon, bmlat, plot, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='bottom', pad="5%")
cbar.set_label("change in d18OV, permil")

plt.title(plotwhat + " diffs at level " + level + ' (test-control) \n' + 'test file: '+ testfn + '\n' + 'control file: ' + cntrlfn)
if showfig:
	plt.show()
if savefig:
	plt.savefig("/Users/kyleniezgoda/Desktop/" + testfn + "-" + cntrlfn + ".ps", bbox_inches='tight', dpi = 500)
