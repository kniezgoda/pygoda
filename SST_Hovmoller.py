#! /Users/kyleniezgoda/anaconda/bin/python
import numpy as np
from netCDF4 import Dataset
from UsefulNetcdfFunctions import find_indices
from pygoda import findClimoFile
import matplotlib.pyplot as plt

north = 5
south = -5 
east = 270
west = 180

cdir = "PI"
tdir = "MH"

# c = findClimoFile(cdir + "/tos_Omon_CCSM4_piControl_r1i1p1_080001-130012.nc")[0]
c = findClimoFile(cdir + "/allMonths*")[0]
control = Dataset(c, mode = 'r')
lat = control.variables['lat'][:,1]
lon = control.variables['lon'][1,:]

# t = findClimoFile(tdir + "/tos_Omon_CCSM4_midHolocene_r1i1p1_100001-130012.nc")[0]
t = findClimoFile(tdir + "/allMonths*")[0]
test = Dataset(t, mode = 'r')

box = find_indices([south, north, west, east], lat, lon)
c_sst = control.variables["tos"][:,box[0]:box[1], box[2]:box[3]].squeeze()
c_sst_latavg = np.mean(c_sst, axis = 1)

t_sst = test.variables["tos"][:,box[0]:box[1], box[2]:box[3]].squeeze()
t_sst_latavg = np.mean(t_sst, axis = 1)

plt.subplot(1,2,1)

plt.imshow(c_sst_latavg,interpolation='nearest',aspect='auto',cmap=plt.cm.get_cmap('coolwarm'))
plt.ylabel('month')
plt.xlabel('lon')
numlons = len(range(box[2], box[3]))
numtime = c_sst.shape[0]
loc = range(0, numlons, 5)
lons = np.linspace(west, east, len(loc))
lons = [int(round(l)) for l in lons]
plt.xticks(loc, lons)
plt.title("Pre-industrial")

plt.subplot(1,2,2)

plt.imshow(t_sst_latavg,interpolation='nearest',aspect='auto',cmap=plt.cm.get_cmap('coolwarm'))
plt.ylabel('month')
plt.xlabel('lon')
numlons = len(range(box[2], box[3]))
numtime = t_sst.shape[0]
loc = range(0, numlons, 5)
lons = np.linspace(west, east, len(loc))
lons = [int(round(l)) for l in lons]
plt.xticks(loc, lons)
plt.title("mid-Holocene")
plt.show()

