#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python
from pygoda import camgoda, findClimoFile
import numpy as np
#from scipy.stats import pearsonr as regress
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as bm
import argparse


'''
Author: Kyle Niezgoda, April 28, 2019

This code produces maps of correlation arrays between the local variable (-lv) at a point (-latlon) 
and the field variable (-fv) over the globe for a specified time (-dates)
'''
#########################
# Read in the arguments #
#########################

parser = argparse.ArgumentParser()
parser.add_argument('-cdir', '--control_directory', dest = 'cdir', default = 'fc5.2deg.wiso.piControl_kn028/atm/hist')
parser.add_argument('-tdir', '--test_directory', dest = 'tdir', default = 'fc5.2deg.wiso.mh6ka_kn032/atm/hist')
parser.add_argument('-r', '--region', dest = 'region', default = None)
parser.add_argument('-grep', dest = 'grep', default = 'cam.h0')
parser.add_argument('-lats', dest = 'lats', nargs = 2, default = [-50,50])
parser.add_argument('-lons', dest = 'lons', nargs = 2, default = [335,333])
parser.add_argument('-levs', dest = 'levs', type = float, nargs = 6, default = None)
parser.add_argument('-v', '--variables', dest = 'variables', nargs= 2, default = ["PRECT", "PRECT_d18O"])
parser.add_argument('-alpha', '--alpha', dest = 'alpha', default = .5)
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
cdir = ARGS.cdir
tdir = ARGS.tdir
levs = ARGS.levs
color_var, line_var = [v for v in ARGS.variables]
grep = ARGS.grep
region = ARGS.region
savefig = ARGS.savefig
showfig = ARGS.showfig
alpha = float(ARGS.alpha)
if ARGS.developer_mode:
    print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
    savfig = False
    showfig = True


#############
# Main code #
#############

# Set the lat bounds
# Default global tropics
region_name = "Box"
southern_lat, northern_lat = [int(l) for l in ARGS.lats]
left_lon, right_lon = [int(l) for l in ARGS.lons]

# Global
if (region == "GT") | (region == "GlobalTropics"):
	region_name = "GlobalTropics"
	southern_lat = -50
	northern_lat = 50
	left_lon = 0
	right_lon = 360
# Indian monsoon
if (region == "IM") | (region == "IndianMonsoon"):
	region_name = "IndianMonsoon"
	southern_lat = -5
	northern_lat = 45
	left_lon = 40
	right_lon = 110

# Greater maritime continent
if (region == "MC") | (region == "MaritimeContinent"):
	region_name = "MaritimeContinent"
	southern_lat = -20
	northern_lat = 20
	left_lon = 80
	right_lon = 160

# North Africa
if (region == "NA") | (region == "NorthAfrica"):
	region_name = "NorthAfrica"
	southern_lat = -20
	northern_lat = 45
	left_lon = -30
	right_lon = 70

#Central Equatorial Pacific to Western eq. Atlantic
if (region == "EP") | (region == "TropicalOceans"):
	region_name = "TropicalOceans"
	southern_lat = -25
	northern_lat = 25
	left_lon = 180
	right_lon = 355

#South America and the Carribean 
if (region == "SA") | (region == "SouthAmerica"):
	region_name = "SouthAmerica"
	southern_lat = -30
	northern_lat = 30
	left_lon = -120
	right_lon = -20

# Four proxy region in the warm pool
if (region == "4P") | (region == "FourProxies"):
	region_name = "FourProxies"
	southern_lat = -25
	northern_lat = 30
	left_lon = 85
	right_lon = 160

box = (southern_lat, northern_lat, left_lon, right_lon)

cpath, cfilename = findClimoFile('*' + grep + "*", directory = cdir)
tpath, tfilename = findClimoFile('*' + grep + "*", directory = tdir)
cnc = camgoda(cpath)
tnc = camgoda(tpath)

cnc.setBox(box)
data = np.zeros(shape = (len(cnc.boxlat), len(cnc.boxlon), 2, 2))
data[:,:,:,0] = cnc.ExtractData(color_var+','+line_var, box, returnData = True)
data[:,:,:,1] = tnc.ExtractData(color_var+','+line_var, box, returnData = True)

lats = cnc.boxlat
lons = cnc.boxlon
llcrnlat, urcrnlat, llcrnlon, urcrnrlon = [lats[0], lats[-1], lons[0], lons[-1]]
if 0 in lons[1:-2]: # if we cross the gml
	llcrnlon = lons[0]-360

color_lev = np.linspace(-5,5,11)
line_lev = np.linspace(-5,5,11)
if levs is not None:
	color_lev = np.linspace(levs[0], levs[1], levs[2])
	line_lev = np.linspace(levs[3], levs[4], levs[5])

bmlon, bmlat = np.meshgrid(lons, lats)
m = bm(projection = 'cea', llcrnrlat=llcrnlat,urcrnrlat=urcrnlat, llcrnrlon=llcrnlon,urcrnrlon=urcrnrlon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, data[:,:,0,1] - data[:,:,0,0], color_lev, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cont = m.contour(bmlon, bmlat, data[:,:,1,1] - data[:,:,1,0], line_lev, shading = 'flat', latlon = True, colors = 'k', alpha = alpha)
cbar = m.colorbar(cs, location='right', pad="5%")
clabel = plt.clabel(cont, fmt='%2.1f', colors = 'k', fontsize = 10)

plt.show()
