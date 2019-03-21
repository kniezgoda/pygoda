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
parser.add_argument('-grep', dest = 'grep', default = None)
parser.add_argument('-box', dest = 'box', nargs = 4, default = [-90,90,0,360])
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variables', dest = 'variables', nargs= "*", default = None)
parser.add_argument('-t', '--test', dest = 'testdatafname', default = None)
parser.add_argument('-c', '--control', dest = 'controldatafname', default = None)
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
region = ARGS.region
season = '' # This needs to be set to comply with legacy coding schemes
grep = ARGS.grep
testdatafname = ARGS.testdatafname
controldatafname = ARGS.controldatafname
findFile = True
if (testdatafname is not None) & (controldatafname is not None):
	findFile = False
testdir = ARGS.testdir
controldir = ARGS.controldir
mkdir = True
savefig = ARGS.savefig # default is True
showfig = ARGS.showfig # default is false
variables = ARGS.variables
if showfig and (not savefig):
	mkdir = False
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savefig = False
	showfig = True
	mkdir = False
	
# Set the lat bounds
# Default arguments
region_name = "Box"
southern_lat, northern_lat ,left_lon, right_lon= [int(l) for l in ARGS.box]

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
	if not os.path.exists("PressureVsLat"):
		os.mkdir("PressureVsLat")
		print "Created directory " + "PressureVsLat"

	# Create the region directory if it doesn't already exist
	if not os.path.exists("PressureVsLat/" + region_name):
		os.mkdir("PressureVsLat/" + region_name)
		print "Created directory " + "PressureVsLat/" + region_name

	# Create grep directory inside region directory
	if not os.path.exists("PressureVsLat/" + region_name + "/" + grep):
		os.mkdir("PressureVsLat/" + region_name + "/" + grep)
		print "Created directory " + "PressureVsLat/" + region_name + "/" + grep

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
	os.chdir("PressureVsLat/" + region_name + "/" + grep)

for var in variables:
	print "\nPlotting " + var + " data...\n"
	data = np.zeros(shape = (len(pressures), 2))
	for p_idx, p in enumerate(pressures):
		print p
		extractVar = "3d_"+var+"_"+str(p/100)
		data[p_idx, 0] = np.nanmean(control.ExtractData(extractVar, box, returnData = True))
		data[p_idx, 1] = np.nanmean(test.ExtractData(extractVar, box, returnData = True))

	aty = np.arange(len(pressures), step = 5) 
	laby = np.array(pressures)[aty]/100
	# labx = np.array(lat)
	# atx = np.linspace(0, len(labx)-1, num = 8)
	# labx = [round(l) for l in labx[np.array([int(round(a)) for a in atx])]]

	fig = plt.figure()

	plt.subplot(2,1,1)
	plt.plot(data[:,0], np.array(pressures)/100, label = "control")
	plt.plot(data[:,1], np.array(pressures)/100, label = "test")
	#plt.yticks(aty, laby)
	plt.gca().invert_yaxis()
	plt.legend()

	plt.subplot(2,1,2)
	plt.plot(data[:,1]-data[:,0], np.array(pressures)/100)
	plt.gca().invert_yaxis()

	# fig.suptitle(control.long_name, fontweight = 'bold', fontsize = 14)
	# plt.figtext(.25, .02,"Averaged over longitudes " + str(left_lon) + " and " + str(right_lon), fontsize = 8)
	# plt.subplots_adjust(hspace = .3)
	if showfig:
		plt.show()
	if savefig:
		plt.savefig(var + "_" + testfn + "-" + controlfn + ".ps", bbox_inches='tight', dpi = 500)
		print("Created " + var + "_" + testfn + "-" + controlfn + ".ps" + '\n')
