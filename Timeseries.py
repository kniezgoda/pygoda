#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python

'''

This function is intended to be used at the command line

Required cmd line flags:
-v VAR1 , (opt) VAR2... : the variables to plot the timeseries of
	multiple variables willbe plotted on separate subplots in order

Optional cmd line flags:
-cen CENTER : the center lat lon
-del DELTA : the latlon delta to move around the center; default is 5 deg
-r BOTTOM_LAT TOP_LAT LEFT_LON RIGHT_LON : the lat lon bounds
	if -cen flag exists, this option is disabled
-dir DIRECTORY : the directories to look for data in. 
	If this flag exists, multiple timeseries will be shown on the same subplot
-run RUNNING_AVG : how many days to average the data over in a running average fashion

Control flags
-dev : prints all plots to screen, does not save any plots, no directories created
-nosave : do not save plots
-show : show the plots to the screen

'''

from pygoda import camdates, findClimoFile, camgoda, corr
import numpy as np
import matplotlib.pyplot as plt
import argparse

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N 

parser = argparse.ArgumentParser()
parser.add_argument('-years', '--years', dest = 'years', nargs = 2, default = None)
parser.add_argument('-months', '--months', dest = 'months', nargs = '*', default = [1,2,3,4,5,6,7,8,9,10,11,12])
parser.add_argument('-days', dest = 'h1', action = 'store_true')
parser.add_argument('-cen', '--center_latlon', dest = 'center_latlon', nargs = 2, default = None)
parser.add_argument('-del', '--delta_latlon', dest = 'delta_latlon', default = 5)
parser.add_argument('-box', dest = 'box', nargs = 4, default = None)
parser.add_argument('-grep', dest = 'grep', default = "")
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variables', nargs= "*")
parser.add_argument('-dir', '--directory', dest = 'directory', default = '.')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')
parser.add_argument('-run', '--running_mean', dest = 'running_mean', default = 1)

##########################
# Read command-line args #
##########################

ARGS = parser.parse_args()
run = int(ARGS.running_mean)
delta = int(ARGS.delta_latlon)
grep = ARGS.grep
if ARGS.center_latlon is not None:
	lat, lon = [int(c) for c in ARGS.center_latlon]
	box = [lat-delta, lat+delta, lon-delta, lon+delta]
elif ARGS.box is not None:
	box = [int(b) for b in ARGS.box]
else:
	print "No region set! Defaulting to global"
	box = [-90,90,0,360]

print "Box array:"
print box

# Find the date array
dates = []
h1 = ARGS.h1 # boolean
months = [int(m) for m in ARGS.months]
if ARGS.years is not None:
	y1, y2 = [int(y) for y in ARGS.years]
	dates = camdates(y1, y2, months, days = h1)
else:
	y1 = y2 = 0
	# Remove the years and leave only months and days if needed
	# Useful for *dailyClim* files
	dates = [x[5:] for x in camdates(y1, y2, months, days = h1)]

print "Date array:" 
print dates

savefig = ARGS.savefig
showfig = ARGS.showfig
variables = ARGS.variables
mkdir = True
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savfig = False
	showfig = True
	mkdir = False

##################
# Main algorithm #
##################

# Creates the master array of the correct shape
var_master = np.zeros(shape = (len(dates), len(variables)))
long_name = []
units = []
for n, d in enumerate(dates):
	# Find the file for this date
	full_path, fname = findClimoFile("*"+grep+"*"+d+"*", directory = ARGS.directory)
	if fname != 0:
		print fname
	# Open the file
	nc = camgoda(full_path)
	for m, v in enumerate(variables):
		# Read the data
		var_is_3d, var, pressure = nc.ExtractData(v, box)
		data = nc.data
		# Average the data
		data_avg = np.nanmean(data)
		var_master[n,m] = data_avg
		if n == 0:
			long_name.append(nc.long_name)
			units.append(nc.units)

ntime, nvar = var_master.shape
add = 0
if nvar > 1:
	add = 1
	
for i in range(nvar):
	plt.subplot(nvar+add,1,i+1)
	plt.plot(var_master[:,i])
	plt.title(long_name[i])
	plt.ylabel(units[i])
	atx = [int(round(DATE)) for DATE in np.linspace(0, len(dates)-1, num = 10)]
	labx = ['' for whatever in atx]
	if i == nvar-1:
		labx = np.array(dates)[np.array(atx)]
	plt.xticks(atx,labx,rotation=45)

if add:
	print "Correlating the first two variables..."
	varcorr = [corr(var_master[:,0], var_master[:,1], lag = x) for x in range(-100,101)]
	plt.subplot(nvar+1,1,nvar+1)
	plt.plot(range(-100,101), varcorr)

plt.tight_layout()

if showfig:
	plt.show()
