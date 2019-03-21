#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python

from pygoda import camgoda, findClimoFile, niceClev
import os, sys
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--region', dest = 'region', default = None)
parser.add_argument('-cdir', '--control_directory', dest = 'controldir', default = "F.C5.2deg.wiso.defaultSSTICE_kn002")
parser.add_argument('-tdir', '--test_directory', dest = 'testdir', default = "F.C5.2deg.wiso.obs6kSST_kn003")
parser.add_argument('-t', '--test', dest = 'testdatafname', default = None)
parser.add_argument('-c', '--control', dest = 'controldatafname', default = None)
parser.add_argument('-grep', dest = 'grep', default = None)
parser.add_argument('-box', dest = 'box', nargs = 4, default = [-90,90,0,360])
parser.add_argument('-clev', dest = 'clev', type = float, nargs = 3, default = None)
parser.add_argument('-diffclev', dest = 'diffclev',type = float, nargs = 3, default = None)
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variables', dest = 'variables', nargs= "*", default = None)
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
region = ARGS.region
season = '' # This needs to be set to comply with legacy coding schemes
testdatafname = ARGS.testdatafname
controldatafname = ARGS.controldatafname
findFile = True
if (testdatafname is not None) & (controldatafname is not None):
	findFile = False
grep = ARGS.grep
testdir = ARGS.testdir
controldir = ARGS.controldir
savefig = ARGS.savefig # default is True
showfig = ARGS.showfig # default is false
mkdir = True
if showfig and not savefig:
	mkdir = False
variables = ARGS.variables
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savefig = False
	showfig = True
	mkdir = False
clev = ARGS.clev
diffclev = ARGS.diffclev

# Set the lat bounds
# Default arguments
region_name = "Box"
southern_lat, northern_lat, left_lon, right_lon = [int(l) for l in ARGS.box]

# Global
if (region == "GT") | (region == "GlobalTropics"):
	region_name = "GlobalTropics"
	southern_lat = -50
	northern_lat = 50
	left_lon = 0
	right_lon = 360
# Arabian Sea
if (region == "AS") | (region == "ArabianSea"):
	region_name = "ArabianSea"
	southern_lat = 0
	northern_lat = 30
	left_lon = 50
	right_lon = 75

# Bay of Bengal
if (region == "BB") | (region == "BayOfBengal"):
	region_name = "BayOfBengal"
	southern_lat = 0
	northern_lat = 30
	left_lon = 75
	right_lon = 95

# Greater maritime continent
if (region == "MC") | (region == "MaritimeContinent"):
	region_name = "MaritimeContinent"
	southern_lat = -10
	northern_lat = 10
	left_lon = 80
	right_lon = 160

# North Africa
if (region == "NA") | (region == "NorthAfrica"):
	region_name = "NorthAfrica"
	southern_lat = 5
	northern_lat = 15
	left_lon = -30
	right_lon = 70

#Central Equatorial Pacific to Western eq. Atlantic
if (region == "EP") | (region == "TropicalOceans"):
	region_name = "TropicalOceans"
	southern_lat = -10
	northern_lat = 10
	left_lon = 180
	right_lon = 355

# Region containing 4 proxies near the warm pool
if (region == "4P") | (region == "FourProxies"):
	region_name = "FourProxies"
	southern_lat = -25
	northern_lat = 30
	left_lon = 85
	right_lon = 160

box = [southern_lat, northern_lat, left_lon, right_lon]
pressures = range(100000, 0, -5000)

if mkdir:
	# Create maps directory is it doesn't exist
	if not os.path.exists("PressureVsLong"):
		os.mkdir("PressureVsLong")
		print "Created directory " + "PressureVsLong"

	# Create the region directory if it doesn't already exist
	if not os.path.exists("PressureVsLong/" + region_name):
		os.mkdir("PressureVsLong/" + region_name)
		print "Created directory " + "PressureVsLong/" + region_name

	# Create grep directory inside region directory
	if not os.path.exists("PressureVsLong/" + region_name + "/" + grep):
		os.mkdir("PressureVsLong/" + region_name + "/" + grep)
		print "Created directory " + "PressureVsLong/" + region_name + "/" + grep

if findFile:
	# Look for the climo files in the root directory
	print "\nLooking for control " + grep + " files in " + controldir
	controldatafname, controlfn = findClimoFile("*" + grep + "*", controldir)
	if not controldatafname:
		sys.exit()
	print "Found control file: " + controlfn
	print "\nLooking for test " + grep + " files in " + testdir
	testdatafname, testfn = findClimoFile("*" + grep + "*", testdir)
	if not testdatafname:
		sys.exit()
	print "Found test file: " + testfn
else:
	print "\nControl file is " + controldatafname
	print "\nTest file is " + testdatafname
	controlfn = os.path.splitext(os.path.split(controldatafname)[1])[0]
	testfn = os.path.splitext(os.path.split(testdatafname)[1])[0]
	
# Read the data
control = camgoda(controldatafname)
test = camgoda(testdatafname)

if mkdir:
	os.chdir("PressureVsLong/" + region_name + "/" + season)

for var in variables:
	print "\nPlotting " + var + " data...\n"

	control.ExtractData(var, box)
	test.ExtractData(var, box)
	lon = control.boxlon
	nlon = len(lon)
	data = np.zeros(shape = (len(pressures), nlon, 2))
	
	for p_idx, p in enumerate(pressures):
		print p
		# Average down to 1 horizontal dimension
		data[p_idx,:,0] = np.nanmean(control.isobar(p, setData = False), axis = 0)
		data[p_idx,:,1] = np.nanmean(test.isobar(p, setData = False), axis = 0)

	aty = np.arange(len(pressures), step = 5) 
	laby = np.array(pressures)[aty]/100
	labx = np.array(lon)
	atx = np.linspace(0, len(labx)-1, num = 8)
	labx = labx[np.array([int(round(a)) for a in atx])] 

	fig = plt.figure()

	test.prep_map(season,region)
	control.prep_map(season, region)
	testclev = niceClev(data[...,1])
	controlclev = niceClev(data[...,0])
	dlev = niceClev(data[...,1]-data[...,0])

	if clev is not None:
		testclev = np.linspace(clev[0], clev[1], clev[2])
		controlclev = np.linspace(clev[0], clev[1], clev[2])	
		
	if diffclev is not None:
		dlev = np.linspace(diffclev[0], diffclev[1], diffclev[2]) 

	plt.subplot(3,1,1)
	tplot = plt.contourf(data[...,1],testclev, cmap = test.cmap)
	# plt.title("mh", fontsize = 8)
	plt.xticks(atx, labx)
	plt.yticks(aty, laby)
	cbar = plt.colorbar(tplot)
	cbar.set_label(control.units)

	plt.subplot(3,1,2)
	cplot = plt.contourf(data[...,0],controlclev, cmap = control.cmap)
	# plt.title("pi", fontsize = 8)
	plt.xticks(atx, labx)
	plt.yticks(aty, laby)
	cbar = plt.colorbar(cplot)
	cbar.set_label(control.units)

	plt.subplot(3,1,3)
	dplot = plt.contourf(data[...,1]-data[...,0], dlev, cmap = control.diffcmap)
	# plt.title("diff", fontsize = 8)
	plt.xticks(atx, labx)
	plt.yticks(aty, laby)
	cbar = plt.colorbar(dplot)
	cbar.set_label(control.units)

	fig.suptitle(control.long_name, fontweight = 'bold', fontsize = 14)
	plt.figtext(.25, .02,"Averaged over latitudes " + str(southern_lat) + " and " + str(northern_lat), fontsize = 8)
	plt.subplots_adjust(hspace = .3)
	if showfig:
		plt.show()
	if savefig:
		plt.savefig(var + "_" + testfn + "-" + controlfn + ".ps", bbox_inches='tight', dpi = 500)
		print("Created " + var + "_" + testfn + "-" + controlfn + ".ps" + '\n')
