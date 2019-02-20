#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python
from pygoda import camgoda, camdates, findClimoFile, corr_2d, boxOut
import numpy as np
from scipy.stats import pearsonr as regress
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
parser.add_argument('-dir', '--directory', dest = 'directory', default = 'fc5.2deg.wiso.piControl_kn028/atm/hist')
parser.add_argument('-loc_var', '--local_variable', dest = 'lv', default = None)
parser.add_argument('-field_var', '--field_variable', dest = 'fv', default = None)
parser.add_argument('-latlon', dest = 'latlon', nargs = 2, default = (0,0))
parser.add_argument('-del', dest = 'delta', default = 0)
parser.add_argument('-box', dest = 'box', nargs = 4, default = [-40,40,335,333])
parser.add_argument('-years', dest = 'years', nargs = 2, default = [10,30])
parser.add_argument('-months', dest = 'months', nargs = '*', default = [1,2,3,4,5,6,7,8,9,10,11,12])
parser.add_argument('-grep_pre', dest = 'grep_pre', default = '')
parser.add_argument('-grep_post', dest = 'grep_post', default = '')
parser.add_argument('-alpha', '--alpha', dest = 'alpha', default = .05)
parser.add_argument('-stiple', '--stiple', dest = 'stiple', action = "store_true")
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
directory = ARGS.directory
lv = ARGS.lv
fv = ARGS.fv
lat, lon = [int(x) for x in ARGS.latlon]
delta = int(ARGS.delta)
start, end = [int(x) for x in ARGS.years]
months = [int(m) for m in ARGS.months]
grep_pre = ARGS.grep_pre
grep_post = ARGS.grep_post
savefig = ARGS.savefig
showfig = ARGS.showfig
alpha = float(ARGS.alpha)
stiple = ARGS.stiple
if ARGS.developer_mode:
    print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
    savfig = False
    showfig = True


#############
# Main code #
#############

# Make the date array
dates = camdates(start, end, months)

# Read in each file, extract the local variable at the point and field variable on the globe, track them in a master array
# Correlate the local variable 1-d array to each grid cell of the field variable array
# Return the correlation coefficients 

bottom = lat-delta
top = lat+delta
left = lon-delta
if left < 0:
	left += 360
right = lon+delta
if right > 360:
	right -= 360
	
region = [bottom, top, left, right]
box = [int(l) for l in ARGS.box]
southern_lat, northern_lat, left_lon, right_lon = box

for date_idx, date in enumerate(dates):
	# Find the file
  path, filename = findClimoFile(grep_pre+'*'+date+'*'+grep_post,directory)
  print filename
	# Open the file
  nc = camgoda(path)
  if date_idx == 0:
    data = np.zeros(shape = (len(dates), len(nc.lat), len(nc.lon), 2))
  data[date_idx, ...] = nc.ExtractData(lv+","+fv, returnData = True)
	#   
# Extract the lv timeseries 
lv_ts = np.nanmean(boxOut(data[...,0], region), axis = (-2,-1))

# Box out the fv 
fv_boxed_ts = boxOut(data[...,1], box)

corr_array = corr_2d(np.expand_dims(np.expand_dims(lv_ts,-1),-1), fv_boxed_ts, axis = 0)

# Make the map of the corr array
fig = plt.figure()

bmlon, bmlat = np.meshgrid(lons, lats)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
clevs = np.linspace(-1, 1, 21)
cs = m.contourf(bmlon, bmlat, corr, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label("correlation-coefficient", fontsize = 8)
plt.title(lv + " vs global " + fv)
x,y = m(lon,lat)
m.plot(x, y, 'gx')

if savefig:
	plt.savefig("correlationMap.pdf")
if showfig:
	plt.show()

