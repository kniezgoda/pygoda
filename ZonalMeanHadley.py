#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python
from pygoda import ncgoda, findClimoFile, niceClev
import matplotlib.pyplot as plt
import numpy as np
import os, argparse
root = os.getcwd()

# Read in cmd line args
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--season', dest = 'season', default = 'ANN')
parser.add_argument('-cdir', '--control_directory', dest = 'controldir', default = "F.C5.2deg.wiso.defaultSSTICE_kn002")
parser.add_argument('-tdir', '--test_directory', dest = 'testdir', default = "F.C5.2deg.wiso.obs6kSST_kn003")
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')
parser.add_argument('-dots', '--dots_on_plot', dest = 'plotdots', action = 'store_true')
parser.add_argument('-lons', '--longitudes', dest = 'longitudes', nargs= 2, default = None)

ARGS = parser.parse_args()
season = ARGS.season
print "\nSeason is " + season
testdir = ARGS.testdir
controldir = ARGS.controldir
savefig = ARGS.savefig
showfig = ARGS.showfig
plotdots = ARGS.plotdots
if ARGS.longitudes is None:
	box = None
	region = "global"
else:
	l1 = int(ARGS.longitudes[0])
	l2 = int(ARGS.longitudes[1])
	box = [-90,90,l1,l2]
	if l1 % 360 > 180:
		r1 = str(360 - (l1 % 360)) + "W-"
	else:
		r1 = str(l1) + "E-"
	if l2 % 360 > 180:
		r2 = str(360 - (l2 % 360)) + "W"
	else:
		r2 = str(l2) + "E"
	region = r1+r2
print "Region is " + region
mkdir = True
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savefig = False
	showfig = True
	mkdir = False

# Look for the climo files in the root directory
print "\nLooking for control " + season + " files in " + controldir + "..."
control_fullpath, controlfn = findClimoFile(season, controldir)
if not control_fullpath:
	sys.exit()
else:
	print "Found file " + controlfn
print "\nLooking for test " + season + " files in " + testdir + "..."
test_fullpath, testfn = findClimoFile(season, testdir)
if not test_fullpath:
	sys.exit()
else:
	print "Found file " + testfn

if mkdir:
	# Create maps directory is it doesn't exist
	if not os.path.exists("ZonalMeanHadley"):
		os.mkdir("ZonalMeanHadley")
		print "Created directory " + "ZonalMeanHadley"

	# Create the region directory if it doesn't already exist
	if not os.path.exists("ZonalMeanHadley/" + season):
		os.mkdir("ZonalMeanHadley/" + season)
		print "Created directory " + "ZonalMeanHadley/" + season

# Read in the data and run HadleyCellInfo
c = ncgoda(control_fullpath)
t = ncgoda(test_fullpath)
cZonalPsi, cCellInfo = c.HadleyCellInfo(box)
tZonalPsi, tCellInfo = t.HadleyCellInfo(box)
diffCellInfo = {k : tCellInfo[k] - cCellInfo[k] for k in cCellInfo.keys()}

# Echo some text to the screen
if showfig:
	print 																	   			\
		'\nTest\n----\n' 																\
			'ITCZ Center : ' + str(tCellInfo['ITCZcenter']) + '\n'   			    	\
			'Northern Edge : ' + str(tCellInfo['HadleyCellNorthernEdge']) + '\n'     	\
			'Southern Edge : ' + str(tCellInfo['HadleyCellSouthernEdge']) + '\n\n'    	\
		'Control\n-------\n' 															\
			'ITCZ Center : ' + str(cCellInfo['ITCZcenter']) + '\n'   			    	\
			'Northern Edge : ' + str(cCellInfo['HadleyCellNorthernEdge']) + '\n'     	\
			'Southern Edge : ' + str(cCellInfo['HadleyCellSouthernEdge']) + '\n\n'    	\
		'Test-Control difference\n-----------------------\n' 				   			\
			'ITCZ Center : ' + str(diffCellInfo['ITCZcenter']) + '\n'   			   	\
			'Northern Edge : ' + str(diffCellInfo['HadleyCellNorthernEdge']) + '\n'   	\
			'Southern Edge : ' + str(diffCellInfo['HadleyCellSouthernEdge']) + '\n'  	

if savefig:
	import pandas as pd
	InfoDataFrame = pd.DataFrame(columns=('Case', 'ITCZcenter', 'HadleyCellNorthernEdge', 'HadleyCellSouthernEdge'))
	tDF = pd.DataFrame({'Case' : 'Test, MH', 'ITCZcenter' : tCellInfo['ITCZcenter'],  'HadleyCellNorthernEdge' : tCellInfo['HadleyCellNorthernEdge'], 'HadleyCellSouthernEdge' : tCellInfo['HadleyCellSouthernEdge']}, index = [0])
	InfoDataFrame = InfoDataFrame.append(tDF)
	cDF = pd.DataFrame({'Case' : 'Control, PI', 'ITCZcenter' : cCellInfo['ITCZcenter'],  'HadleyCellNorthernEdge' : cCellInfo['HadleyCellNorthernEdge'], 'HadleyCellSouthernEdge' : cCellInfo['HadleyCellSouthernEdge']}, index = [0])
	InfoDataFrame = InfoDataFrame.append(cDF)
	diffDF = pd.DataFrame({'Case' : 'Difference, MH-PI', 'ITCZcenter' : diffCellInfo['ITCZcenter'],  'HadleyCellNorthernEdge' : diffCellInfo['HadleyCellNorthernEdge'], 'HadleyCellSouthernEdge' : diffCellInfo['HadleyCellSouthernEdge']}, index = [0])
	InfoDataFrame = InfoDataFrame.append(diffDF)
	InfoDataFrame.to_csv("ZonalMeanHadley/" + season + "/" + region + "_" + season + "_HadleyData.csv")
	print("\nCreated " + "ZonalMeanHadley/" + season + "/" + region + "_" + season + "_HadleyData.csv")

# Print plots
# psi vs lat
cLine, = plt.plot(c.boxlat, cZonalPsi/10**9, color = 'red', label = 'Control')
tLine, = plt.plot(t.boxlat, tZonalPsi/10**9, color = 'blue', label = 'Test')

if plotdots:
	#control
	itczcenter, = plt.plot(cCellInfo['ITCZcenter'], 0, 'ro')
	plt.text(cCellInfo['ITCZcenter']+2, 5, "lat: " + str(round(cCellInfo['ITCZcenter'], 1)), color = 'red')
	southernedge, = plt.plot(cCellInfo['HadleyCellSouthernEdge'], 0, 'ro')
	plt.text(cCellInfo['HadleyCellSouthernEdge']+2, 5, "lat: " + str(round(cCellInfo['HadleyCellSouthernEdge'], 1)), color = 'red')
	northernedge, = plt.plot(cCellInfo['HadleyCellNorthernEdge'], 0, 'ro')
	plt.text(cCellInfo['HadleyCellNorthernEdge']+2, 5, "lat: " + str(round(cCellInfo['HadleyCellNorthernEdge'], 1)), color = 'red')
	#test
	itczcenter, = plt.plot(tCellInfo['ITCZcenter'], 0, 'bo')
	plt.text(tCellInfo['ITCZcenter']+2, -5, "lat: " + str(round(tCellInfo['ITCZcenter'], 1)), color = 'blue')
	southernedge, = plt.plot(tCellInfo['HadleyCellSouthernEdge'], 0, 'bo')
	plt.text(tCellInfo['HadleyCellSouthernEdge']+2, -5, "lat: " + str(round(tCellInfo['HadleyCellSouthernEdge'], 1)), color = 'blue')
	northernedge, = plt.plot(tCellInfo['HadleyCellNorthernEdge'], 0, 'bo')
	plt.text(tCellInfo['HadleyCellNorthernEdge']+2, -5, "lat: " + str(round(tCellInfo['HadleyCellNorthernEdge'], 1)), color = 'blue')

plt.legend(handles=[cLine, tLine], loc = 'lower right')
plt.title('Top-down integrated steam function at 500 mb, ' + season)
plt.xlabel('Latitude')
plt.ylabel('Tg/s')
if showfig:
	plt.show()
if savefig:
	plt.savefig("ZonalMeanHadley/" + season + "/" + region + "_" + season + "_PsiLines_" + testfn + "-" + controlfn + ".ps", bbox_inches='tight', dpi = 500)
	print("\nCreated " + "ZonalMeanHadley/" + season + "/" + region + "_" + season + "_PsiLines_" + testfn + "-" + controlfn + ".ps" + '\n')


# Pressure vs Lat
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

pressures = range(90000, 0, -5000)
master = Niezgoda(len(c.boxlat))

for p in pressures:
	print p
	# Average down to 1 horizontal dimension
	hold_c = np.nanmean(c.isobar(p), axis = 1)
	hold_t = np.nanmean(t.isobar(p), axis = 1)

	master.stack(hold_c, hold_t)

master.finish()

aty = np.arange(len(pressures), step = 5) 
laby = np.array(pressures)[aty]/100
labx = np.array(c.boxlat)
atx = np.linspace(0, len(labx)-1, num = 8)
labx = [round(l) for l in labx[np.array([int(round(a)) for a in atx])]]

fig = plt.figure()

plt.subplot(3,1,1)
plt.title("Test, " + season)
clev = niceClev(master.test/10**9)
tplot = plt.contourf(master.test/10**9, clev, cmap = plt.cm.RdBu_r)
plt.xticks(atx, labx)
plt.yticks(aty, laby)
cbar = plt.colorbar(tplot)
cbar.set_label("Tg/s")

plt.subplot(3,1,2)
plt.title("Control, " + season)
clev = niceClev(master.control/10**9)
cplot = plt.contourf(master.control/10**9, clev, cmap = plt.cm.RdBu_r)
plt.xticks(atx, labx)
plt.yticks(aty, laby)
cbar = plt.colorbar(cplot)
cbar.set_label("Tg/s")

plt.subplot(3,1,3)
plt.title("Test - control difference, " + season)
diff = master.test/10**9 - master.control/10**9
clev = niceClev(diff)
dplot = plt.contourf(diff, clev, cmap = plt.cm.RdBu_r)
plt.xticks(atx, labx)
plt.yticks(aty, laby)
cbar = plt.colorbar(dplot)
cbar.set_label("Tg/s")

fig.suptitle("Top-down integrated meridional stream function", fontweight = 'bold', fontsize = 14)
plt.subplots_adjust(hspace = .5)

if showfig:
	plt.show()
if savefig:
	plt.savefig("ZonalMeanHadley/" + season + "/" + region + "_" + season + "_PsiContour_" + testfn + "-" + controlfn + ".ps", bbox_inches='tight', dpi = 500)
	print("\nCreated " + "ZonalMeanHadley/" + season + "/" + region + "_" + season + "_PsiContour_" + testfn + "-" + controlfn + ".ps" + '\n')