#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python
from pygoda import h1ts
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--region', dest = 'region', nargs = 4, default = None)
parser.add_argument('-cen', '--center_latlon', dest = 'center_latlon', nargs = 2, default = None)
parser.add_argument('-del', '--delta_latlon', dest = 'delta_latlon', nargs = 1, default = 5)
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variables', nargs= "*")
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')


ARGS = parser.parse_args()

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

print "\n\nRegion is " + str(box) 
savefig = ARGS.savefig
showfig = ARGS.showfig
variables = ARGS.variables
mkdir = True
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savfig = False
	showfig = True
	mkdir = False

for n, v in enumerate(variables):
	
	varhold, long_name, units, dates = h1ts(v, box)

	plt.subplot(len(variables),1,n+1)
	plt.plot(varhold)
	plt.title(long_name)
	plt.ylabel(units)
	if v != variables[-1]:
		plt.tick_params(axis='x', which='both', top='off', labelbottom='off')
	if v == variables[-1]:
		atx = [int(round(DATE)) for DATE in np.linspace(0, len(dates)-1, num = 10)]
		labx = np.array(dates)[np.array(atx)]
		plt.xticks(atx,labx,rotation=45)

if showfig:
	plt.show()
if savefig:
	None
