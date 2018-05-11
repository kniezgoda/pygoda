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
from pygoda import h1ts
import numpy as np
import matplotlib.pyplot as plt
import argparse

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N 

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--region', dest = 'region', nargs = 4, default = None)
parser.add_argument('-cen', '--center_latlon', dest = 'center_latlon', nargs = 2, default = None)
parser.add_argument('-del', '--delta_latlon', dest = 'delta_latlon', nargs = 1, default = 5)
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variables', nargs= "*")
parser.add_argument('-dir', '--directory', dest = 'directory', nargs= "*", default = ['.'])
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')
parser.add_argument('-run', '--running_mean', dest = 'running_mean', default = 1)


ARGS = parser.parse_args()
run = int(ARGS.running_mean)
wd = ARGS.directory
region = ARGS.region
center = ARGS.center_latlon
delta = ARGS.delta_latlon
if center is not None:
	box = [int(center[0])-delta, int(center[0])+delta, int(center[1])-delta, int(center[1])+delta]
elif region is not None:
	box = [int(r) for r in region]
else:
	print "No region set! Defaulting to global"
	box = [-90,90,0,360]

print "\nRegion is " + str(box) 
savefig = ARGS.savefig
showfig = ARGS.showfig
variables = ARGS.variables
mkdir = True
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savfig = False
	showfig = True
	mkdir = False

fig = plt.figure()
for n, v in enumerate(variables):
	for w in wd:
		varhold, long_name, units, dates = h1ts(v, box, w)
		varhold = running_mean(varhold, run)

		plt.subplot(len(variables),1,n+1)
		if (w == '.') | (w == './'):
			experiment = ""
		else:
			experiment = w
		if len(wd) > 1:
			plt.plot(varhold, label = experiment)
			plt.legend()
		else:
			plt.plot(varhold)
		plt.title(long_name)
		plt.ylabel(units)
		atx = [int(round(DATE)) for DATE in np.linspace(0, len(dates)-1, num = 10)]
		labx = ['' for whatever in atx]
		if v == variables[-1]:
			labx = np.array(dates)[np.array(atx)]
		plt.xticks(atx,labx,rotation=45)

fig.suptitle("Region: " + str(box))

if showfig:
	plt.show()
if savefig:
	None
