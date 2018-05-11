#!/home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python

from pygoda_2 import ExtractData, findClimoFile
import os, sys, glob, argparse
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('-cdir', '--control_directory', dest = 'cdir', default = 'fc5.2deg.wiso.piControl_kn028/AMWG.climoFiles/')
parser.add_argument('-tdir', '--test_directory', dest = 'tdir', default = 'fc5.2deg.wiso.mh6ka_kn032/AMWG.climoFiles/')
parser.add_argument('-grep', dest = 'grep', default = "ANN")
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variable', nargs= 1, default = None)
parser.add_argument('-t', '--test', dest = 'testdatafname', default = None)
parser.add_argument('-c', '--control', dest = 'controldatafname', default = None)
parser.add_argument('-lats', dest = 'lats', default = [-90, 90])
parser.add_argument('-lons', dest = 'lons', default = [0, 360])
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
cdir = ARGS.cdir
tdir = ARGS.tdir
grep = ARGS.grep
cfile = ARGS.controldatafname
tfile = ARGS.testdatafname
bottom_lat, top_lat = [int(l) for l in ARGS.lats]
left_lon, right_lon = [int(l) for l in ARGS.lons]
box = [bottom_lat, top_lat, left_lon, right_lon]
variable = ARGS.variable
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
	controldatafname, controlfn = findClimoFile("*" + grep + "*", cdir)
	if not controldatafname:
		sys.exit()
	else:
		print "Found file " + controlfn
	print "\nLooking for test " + grep + " files in " + tdir + "..."
	testdatafname, testfn = findClimoFile("*" + grep + "*", tdir)
	if not testdatafname:
		sys.exit()
	else:
		print "Found file " + testfn

else:
	print "\nControl file is " + cfile
	print "\nTest file is " + tfile
	controlfn = os.path.splitext(os.path.split(cfile)[1])[0]
	testfn = os.path.splitext(os.path.split(tfile)[1])[0]


# Extract the variable
cnc = ExtractData(controldatafname, variable, box)
tnc = ExtractData(testdatafname, variable, box)
