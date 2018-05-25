#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python

from pygoda import camgoda, camdates, findClimoFile, eof
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as bm 
import os, sys, glob, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-cdir', '--control_directory', dest = 'cdir', default = 'fc5.2deg.wiso.piControl_kn028/atm/hist')
parser.add_argument('-tdir', '--test_directory', dest = 'tdir', default = 'fc5.2deg.wiso.mh6ka_kn032/atm/hist')
parser.add_argument('-v', '--variable', dest = 'variable', default = None)
parser.add_argument('-years', dest = 'years', nargs = 2, default = [None, None])
parser.add_argument('-months', dest = 'months', nargs = '*', default = [1,2,3,4,5,6,7,8,9,10,11,12])
parser.add_argument('-box', dest = 'box', nargs = 4, default = [-50, 50, 0, 360])
parser.add_argument('-n', dest = 'num_eofs', default = 1)
parser.add_argument('-grep_pre', dest = 'grep_pre', default = '')
parser.add_argument('-grep_post', dest = 'grep_post', default = '')
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

# Read the ARGS
ARGS = parser.parse_args()
cdir = ARGS.cdir
tdir = ARGS.tdir
variable = ARGS.variable
box = [int(b) for b in ARGS.box]
num_eofs = int(ARGS.num_eofs)
start, end = ARGS.years
if start is not None and end is not None:
	start = int(start)
	end = int(end)
months = [int(m) for m in ARGS.months]
grep_pre = ARGS.grep_pre
grep_post = ARGS.grep_post
savefig = ARGS.savefig
showfig = ARGS.showfig
if ARGS.developer_mode:
    print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
    savfig = False
    showfig = True

# Extract the dates
dates = camdates(start, end)

for n, date in enumerate(dates):
	# Find the file for this date
	cfull_path, cfname = findClimoFile("*"+grep_pre+"*"+date+"*"+grep_post+"*", cdir)
	tfull_path, tfname = findClimoFile("*"+grep_pre+"*"+date+"*"+grep_post+"*", tdir)
	if cfname != 0:
		print cfname
	if tfname != 0:
		print tfname
	# Open the file
	cnc = camgoda(cfull_path)
	tnc = camgoda(tfull_path)
	is3d, var, vname = cnc.ExtractData(variable, box)
	is3d, var, vname = tnc.ExtractData(variable, box)
	if n == 0:
		nlats, nlons = cnc.data.shape
		boxlat = cnc.boxlat
		boxlon = cnc.boxlon
		d = np.zeros(shape = (len(dates), nlats*nlons))
	d[n,:] = np.ndarray.flatten(tnc.data - cnc.data)

# Compute the amplitude timeseries and EOF spatial distributions of the data array
a, F = eof(d)

# Reshape F into a spatial grid
F_grid = np.reshape(F, (F.shape[0], nlats, nlons))

# Make the maps 
southern_lat, northern_lat, left_lon, right_lon = box
bmlon, bmlat = np.meshgrid(boxlon, boxlat)
if 0 in boxlon[1:-2]: # if we cross the gml
	left_lon = boxlon[0]-360

num_subplots = 2*num_eofs
for subplot in range(num_subplots):
	if subplot % 2 == 0: #even, the F maps
		plt.subplot(2,num_eofs,subplot)
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, F_grid[subplots/2,:,:], shading = 'flat', latlon = True)
		cbar = m.colorbar(cs, location='right', pad="5%")
	else: # odd, the amplitude time series
		atx = np.arange(0, len(dates), 3)
		labx = np.array(dates)[atx]
		plt.subplot(2,num_eofs,subplot)
		plt.plot(a[:,(subplot-1)/2])
		plt.xticks(atx, labx, rotation=45)

plt.show()