#!/Users/kyleniezgoda/anaconda/bin/python

from UsefulNetcdfFunctions import CorrelationArray as CA
from UsefulNetcdfFunctions import ncdump
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as bm
import numpy as np

files = "*cam.h1*"
field_var = "PRECT"
loc_var = "PRECT"
lat = 1.3
lon = 103

corr, lats, lons = CA(files, field_var, loc_var, lat, lon)
bmlon, bmlat = np.meshgrid(lons, lats)

# print ncdump(Dataset("F.C5.2deg.wiso.obs6kSST_kn003.cam.h1.0011-06-29-00000.nc", mode = 'r'))


m = bm(projection = 'cea', llcrnrlat=-90,urcrnrlat=90, llcrnrlon=0,urcrnrlon=360,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
clevs = np.linspace(-np.amax(abs(corr)), np.amax(abs(corr)), 19)
cs = m.contourf(bmlon, bmlat, corr, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label("correlation-coefficient", fontsize = 8)

plt.title("Correlation between\n" + loc_var + " at lat=" + str(lat) + " and lon=" + str(lon) + "\nand " + field_var, fontsize = 8)

plt.show()
#---------------------------------#