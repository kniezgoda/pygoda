#!/home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python

from pygoda import camgoda, findClimoFile
import os, sys, glob, argparse
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('-cdir', '--control_directory', dest = 'cdir', default = 'fc5.2deg.wiso.piControl_kn028/AMWG.climoFiles/')
parser.add_argument('-tdir', '--test_directory', dest = 'tdir', default = 'fc5.2deg.wiso.mh6ka_kn032/AMWG.climoFiles/')
parser.add_argument('-grep', dest = 'grep', default = "ANN")
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variable', nargs= '*', default = None)
parser.add_argument('-t', '--test', dest = 'testdatafname', default = None)
parser.add_argument('-c', '--control', dest = 'controldatafname', default = None)
parser.add_argument('-box', dest = 'box', nargs = 4, default = [-90, 90, 0, 360])
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
cdir = ARGS.cdir
tdir = ARGS.tdir
grep = ARGS.grep
cfile = ARGS.controldatafname
tfile = ARGS.testdatafname
bottom_lat, top_lat, left_lon, right_lon = [int(l) for l in ARGS.box]
box = [bottom_lat, top_lat, left_lon, right_lon]
variables = [str(V) for V in ARGS.variable]
savefig = ARGS.savefig
showfig = ARGS.showfig

mkdir = True
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savefig = False
	showfig = True
	mkdir = False


# Find the file
if (cfile is None) or (tfile is None):
	print "\nLooking for control " + grep + " files in " + cdir + "..."
	cfile, controlfn = findClimoFile("*" + grep + "*", cdir)
	if not cfile:
		sys.exit()
	else:
		print "Found file " + controlfn
	print "\nLooking for test " + grep + " files in " + tdir + "..."
	tfile, testfn = findClimoFile("*" + grep + "*", tdir)
	if not tfile:
		sys.exit()
	else:
		print "Found file " + testfn

else:
	print "\nControl file is " + cfile
	print "\nTest file is " + tfile
	controlfn = os.path.splitext(os.path.split(cfile)[1])[0]
	testfn = os.path.splitext(os.path.split(tfile)[1])[0]

# Open the file
cnc = camgoda(cfile)
tnc = camgoda(tfile)

for v in variables:
	# Extract the variable
	cnc.ExtractData(v, box)
	tnc.ExtractData(v, box)

	# Save the data
	cvar = cnc.data
	tvar = tnc.data
	long_name = cnc.long_name
	units = cnc.units
	lat = cnc.boxlat

	# Zonally average the data
	cvar_zonalMean = np.nanmean(cvar, axis = 1)
	tvar_zonalMean = np.nanmean(tvar, axis = 1)

	# Compute the difference
	dvar_zonalMean = tvar_zonalMean - cvar_zonalMean

	# Plot the data
	fig = plt.figure()

	plt.subplot(211)
	plt.plot(lat, tvar_zonalMean, label = "test")
	plt.plot(lat, cvar_zonalMean, label = "control")
	plt.ylabel(units)
	plt.legend()

	plt.subplot(212)
	plt.plot(lat, dvar_zonalMean, color = 'k')
	plt.ylabel("Difference in " + units)

	fig.suptitle(long_name + "\nAveraged over longitudes " + str(left_lon) + " - " + str(right_lon))

	if savefig:
		plt.savefig(v + "." + str(left_lon) + "_" + str(right_lon) + ".ps")
	if showfig:
		plt.show()
	
