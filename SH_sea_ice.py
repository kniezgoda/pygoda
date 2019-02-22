#Python code for basic isotope graphs
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import shiftgrid 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Set the input file
directory = "/data/kniezgod/IceData/"
fCCSM = "ice_cov.CCSM4b1850-Had1850.nc"
fEE = "ice_cov.ECEARTH-Had.nc"

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Read the netcdf file
CCSMdata = Dataset(directory+fCCSM, mode='r')
EEdata = Dataset(directory+fEE, mode='r')

#Extract lat, lon, and lev data
lats = CCSMdata.variables['lat'][:]
lons = CCSMdata.variables['lon'][:]
bmlon, bmlat = np.meshgrid(lons, lats)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Extract data. Variable names are typical CAM5 output with isotopes enabled
ice_cov1 = CCSMdata.variables['ice_cov'][0,:,:]
ice_cov2 = EEdata.variables['ice_cov'][0,:,:]

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Creation of maps

fig = plt.figure()

ax = fig.add_subplot(211)
ax.set_title("CCSM4 - Hadley")
m = bm(projection='spstere',boundinglat=-40,lon_0=180,resolution='l')
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,15.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, ice_cov1, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
m.colorbar(cs, location='bottom', pad="5%")

ax = fig.add_subplot(212)
ax.set_title("EC-EARTH - Hadley")
m = bm(projection='spstere',boundinglat=-40,lon_0=180,resolution='l')
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,15.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, ice_cov2, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
m.colorbar(cs, location='bottom', pad="5%")


fig.savefig("test.png")
