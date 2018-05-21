#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python

from pygoda import camdates, findClimoFile, camgoda, corr, runningMean
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
parser.add_argument('-grep', dest = 'grep', default = None)
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variables', nargs = 2, dest = 'variable')
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
variables = [v for v in ARGS.variable]
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
var_master = np.zeros(shape = (len(dates), 2))
long_name = []
units = []
for n, d in enumerate(dates):
	# Find the file for this date
	full_path, fname = findClimoFile("*"+grep+"*"+d+"*", directory = ARGS.directory)
	if fname != 0:
		print fname
	# Open the file
	nc = camgoda(full_path)
	for m, v in variable:
		# Read the data
		var_is_3d, var, pressure = nc.ExtractData(v, box)
		data = nc.data
		# Average the data
		data_avg = np.nanmean(data)
		var_master[n,m] = data_avg
		if n == 0:
			long_name.append(nc.long_name)
			units.append(nc.units)

for D in range(2):
	var_master[:,D] = runningMean(var_master[:,D], run)

# Plot the timeseries
plt.subplot(311)
plt.plot(var_master[:,0])
plt.title(long_name[0])
plt.ylabel(units[0])

plt.subplot(312)
pt.plot(var_master[:,1])
plt.title(long_name[1])
plt.ylabel(units[1])
atx = [int(round(DATE)) for DATE in np.linspace(0, len(dates)-1, num = 10)]
labx = np.array(dates)[np.array(atx)]
plt.xticks(atx,labx,rotation=45)

# Compute and plot the ACF
plt.subplot(313)
acf = []
for lag in range(-100,101):
	acf.append(corr(var_master[:,0], var_master[:,1], lag))
plt.plot(range(-100,101), acf)
plt.xlabel("lag")
plt.ylabel("cross-correlation coefficient")

plt.tight_layout()

if showfig:
	plt.show()
