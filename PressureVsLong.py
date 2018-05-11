#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python


from pygoda import ncgoda, findClimoFile, niceClev
import os, sys
import argparse
import matplotlib.pyplot as plt
import numpy as np

use = {           \
	"Q" : 1,      \
	"V" : 1,      \
	"VT" : 1,     \
	"VQ" : 1,     \
	"U" : 1,      \
	"UT" : 1,     \
	"UQ" : 1,     \
	"T" : 1,      \
	"OMEGA" : 1,  \
	"Z3" : 0,     \
	"dDV" : 1,    \
	"d18OV" : 1,  \
	"dxsV" : 1,   \
	"psi" : 1	  \
	}


'''
6 optional args:
-r (--region) REGION : sets the region, default to '' (global tropics)
-s (--season) SEASON : sets season, defaults to 'ANN'
-tdir (test_directory) TDIR : sets the directory to look for test files in, defaults to *_kn003 directory
-cdir (control_directory) CDIR : see above, for control, defaults to *_kn002 directory
-show : sets showfig to True, will print all plots to screen. Default is not to show figs
-nosave : sets savefig to False, will not save figures to current directory. Default is to save images to current directory
'''

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--region', dest = 'region', default = '')
parser.add_argument('-s', '--season', dest = 'season', default = 'ANN')
parser.add_argument('-cdir', '--control_directory', dest = 'controldir', default = "F.C5.2deg.wiso.defaultSSTICE_kn002")
parser.add_argument('-tdir', '--test_directory', dest = 'testdir', default = "F.C5.2deg.wiso.obs6kSST_kn003")
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variable', nargs= "*", default = None)
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
season = ARGS.season
print "Season is " + season
region = ARGS.region
print "Region is " + region
testdir = ARGS.testdir
controldir = ARGS.controldir
savefig = ARGS.savefig # default is True
showfig = ARGS.showfig # default is false
variable = ARGS.variable
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savfig = False
	showfig = True
	mkdir = False

# Set the new variable list if variable is not None
if variable is not None:
	use = {u : 0 for u in use}
	for v in variable:
		if v in use:
			use[v] = 1
		else:
			print "\nVariable " + v + " not in master variable list.\nWill not plot this variable."


# Set the lat bounds
# Default global tropics
region_name = "GlobalTropics"
southern_lat = -10
northern_lat = 10
left_lon = 0
right_lon = 359

# Indian monsoon
if (region == "IM") | (region == "IndianMonsoon"):
	region_name = "IndianMonsoon"
	southern_lat = 20
	northern_lat = 30
	left_lon = 40
	right_lon = 110

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

	# Create season directory inside region directory
	if not os.path.exists("PressureVsLong/" + region_name + "/" + season):
		os.mkdir("PressureVsLong/" + region_name + "/" + season)
		print "Created directory " + "PressureVsLong/" + region_name + "/" + season

# Look for the climo files in the root directory
print "\nLooking for control " + season + " files in " + controldir
controldatafname, controlfn = findClimoFile(season, controldir)
if not controldatafname:
	sys.exit()
print "Found control file: " + controlfn
print "\nLooking for test " + season + " files in " + testdir
testdatafname, testfn = findClimoFile(season, testdir)
if not testdatafname:
	sys.exit()
print "Found test file: " + testfn

# Read the data
control = ncgoda(controldatafname)
test = ncgoda(testdatafname)


# Initialize a variable extraction from ncgoda so that boxlat and boxlon are set to correct values
control.variable("T", box, setData = False)
lon = control.boxlon
numlons = len(lon)

class Niezgoda:
	import numpy as np
	def __init__(self, dimlength):
		self.control = np.zeros(dimlength)
		self.test = np.zeros(dimlength)
	#
	def stack(self, control, test):
		self.control = np.vstack((self.control, control))
		self.test = np.vstack((self.test, test))
	#
	def finish(self):
		self.control = np.delete(self.control, np.where(np.sum(self.control,1)==0)[0][0], 0)
		self.test = np.delete(self.test, np.where(np.sum(self.test,1)==0)[0][0], 0)

if mkdir:
	os.chdir("PressureVsLong/" + region_name + "/" + season)

for var in use:
	if not use.get(var):
		continue
	
	print "\nPlotting " + var + " data...\n"
	
	if var == "d18OV":
		control.d18OV(box)
		test.d18OV(box)
	elif var == "dDV":
		control.dDV(box)
		test.dDV(box)
	elif var == 'dxsV':
		control.dxsV(box)
		test.dxsV(box)
	elif var == "psi":
		control.psi(box)
		test.psi(box)
	else:
		try:
			control.variable(var, box)
			test.variable(var, box)
		except KeyError:
			print "Not able to plot variable " + var + "...\nSkipping this variable."
			continue
	
	master = Niezgoda(numlons)

	for p in pressures:
		# print test.isobar(p)

		print p
		# Average down to 1 horizontal dimension
		hold_c = np.nanmean(control.isobar(p), axis = 0)
		hold_t = np.nanmean(test.isobar(p), axis = 0)

		master.stack(hold_c, hold_t)

	master.finish()

	aty = np.arange(len(pressures), step = 5) 
	laby = np.array(pressures)[aty]/100
	labx = np.array(lon)
	atx = np.linspace(0, len(labx)-1, num = 8)
	labx = labx[np.array([int(round(a)) for a in atx])] 

	fig = plt.figure()

	test.prep_map(season,region)
	control.prep_map(season, region)
	testclev = niceClev(master.test)
	controlclev = niceClev(master.control)
	diffclev = niceClev(master.test - master.control)

	plt.subplot(3,1,1)
	tplot = plt.contourf(master.test,testclev, cmap = test.cmap)
	# plt.title("mh", fontsize = 8)
	plt.xticks(atx, labx)
	plt.yticks(aty, laby)
	cbar = plt.colorbar(tplot)
	cbar.set_label(control.units)

	plt.subplot(3,1,2)
	cplot = plt.contourf(master.control,controlclev, cmap = control.cmap)
	# plt.title("pi", fontsize = 8)
	plt.xticks(atx, labx)
	plt.yticks(aty, laby)
	cbar = plt.colorbar(cplot)
	cbar.set_label(control.units)

	plt.subplot(3,1,3)
	dplot = plt.contourf(master.test - master.control, diffclev, cmap = control.diffcmap)
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