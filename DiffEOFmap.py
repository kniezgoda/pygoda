#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python

from pygoda import camgoda, camdates, findClimoFile
from eofs.standard import Eof
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
parser.add_argument('-days', '--days', dest = 'days', action = 'store_true')
parser.add_argument('-box', dest = 'box', nargs = 4, default = [-50, 50, 0, 360])
parser.add_argument('-n', dest = 'num_eofs', default = 1)
parser.add_argument('-grep_pre', dest = 'grep_pre', default = '')
parser.add_argument('-grep_post', dest = 'grep_post', default = '')
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-raw', dest = 'raw', action = 'store_false')
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
days = ARGS.days
grep_pre = ARGS.grep_pre
grep_post = ARGS.grep_post
savefig = ARGS.savefig
showfig = ARGS.showfig
raw = ARGS.raw
removeMeans = True
if raw:
	removeMeans = False
if ARGS.developer_mode:
    print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
    savfig = False
    showfig = True

# Extract the dates
dates = camdates(start, end, months, days)

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
print "Computing the EOF..."
EOF = Eof(d, center = removeMeans)
eof = EOF.eofs(neofs = num_eofs)
pca = EOF.pcs(npcs = num_eofs, pcscaling = 1)
varfrac = EOF.varianceFraction()
print "Finished!"

# Reshape F into a spatial grid
eof_grid = np.reshape(eof, (eof.shape[0], nlats, nlons))

# Make the maps 
bmlon, bmlat = np.meshgrid(boxlon, boxlat)
southern_lat = boxlat[0]
northern_lat = boxlat[-1]
left_lon = boxlon[0]
right_lon = boxlon[-1]
if 0 in boxlon[1:-2]: # if we cross the gml
	left_lon = boxlon[0]-360

num_subplots = 2*num_eofs
for subplot in range(num_subplots):
	if subplot <= num_eofs-1: #the F maps
		# Generate the clevs
		x = round(np.max(np.abs(eof_grid[subplot,:,:])), 2)
		clev = np.linspace(-x,x,21)
		plt.subplot(2,num_eofs,subplot+1)
		plt.title(round(varfrac[subplot], 3))
		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		cs = m.contourf(bmlon, bmlat, eof_grid[subplot,:,:], clev, shading = 'flat', latlon = True)
		cbar = m.colorbar(cs, location='right', pad="5%")
	else: # the amplitude time series
		atx = np.arange(0, len(dates), 3)
		labx = np.array(dates)[atx]
		plt.subplot(2,num_eofs,subplot+1)
		plt.plot(pca[:,subplot-num_eofs])
		plt.xticks(atx, labx, rotation=45)

plt.show()