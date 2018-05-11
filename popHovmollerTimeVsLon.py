#!/glade/u/apps/ch/opt/python/2.7.13/gnu/6.2.0/bin/python
import argparse
import numpy as np
from pygoda import findClimoFile, popgoda, h0dates, niceClev
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument('-dir', '--directory', dest = 'directory', default = '.')
parser.add_argument('-grep', '--regular_expression', dest = 'regular_expression', default = '*pop.h.')
parser.add_argument('-dates', '--date_bounds', dest = 'date_bounds', nargs = 2)
parser.add_argument('-lats', type = int, dest = 'lats', nargs = 2, default = (-5, 5))
parser.add_argument('-lons', type = int, dest = 'lons', nargs = 2, default = (0, 360))
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variable', nargs= "*", default = None)
parser.add_argument('-ft', '--file_type', dest = 'file_type', default = 'ps')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')
parser.add_argument('-lines', '--draw_lines', dest = 'draw_lines', action = 'store_true')

ARGS = parser.parse_args()
lats = ARGS.lats
lons = ARGS.lons
box = [lats[0], lats[1], lons[0], lons[1]]
dates = ARGS.date_bounds
start = dates[0]
end = dates[1]
lines = ARGS.draw_lines
directory = ARGS.directory
grep = ARGS.regular_expression
savefig = ARGS.savefig
showfig = ARGS.showfig
ftype = ARGS.file_type
mkdir = True
if ARGS.developer_mode:
	print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savefig = False
	showfig = True
	mkdir = False
date_list = h0dates(start, end)
print "Dates : "
print date_list
variable = ARGS.variable

for V in variable:
	print "Variable is : " + V

	var_is_3d = False
	if V[:3] == "3d_":
		var_is_3d = True

	if not var_is_3d:
		var = V
		vname = V
	else:
		var = V[3:-3]
		vname = V[3:]
		pressure = int(V[-3:]) * 100

	for n, date in enumerate(date_list):
		print date

		# Find the file and open it
		nc_path, nc_file = findClimoFile(grep + date + "*", directory)
		nc = popgoda(nc_path)
		
		try:
			nc.variable(var, box)
		except KeyError:
			print "Not able to plot variable " + var + "...\nSkipping this variable."
			continue

		# Average over latitudes
		hold = np.mean(nc.data, 0)
		# Stack things up
		if n == 0:
			units = nc.units
			master = hold
			nlons = len(nc.boxlon)
			nlats = len(nc.boxlat)
		else:
			master = np.vstack((master, hold))

	# fig = plt.figure(figsize=(20,10))

	# clev = niceClev(master, alpha = .9)
	# try:
	# 	vmin = clev[0]
	# 	vmax = clev[-1]
	# except TypeError:
	# 	vmin = None
	# 	vmax = None

	# master_plot = plt.imshow(master, cmap = plt.cm.RdBu_r, vmin = vmin, vmax = vmax, aspect = 'auto')
	master_plot = plt.imshow(master, cmap = plt.cm.RdBu_r, aspect = 'auto')

	plt.ylabel('months since ' + start)
	plt.xlabel('longitude')

	cbar = plt.colorbar(master_plot, orientation = 'horizontal')
	cbar.set_label("")

	numtime = master.shape[0]
	if lines:
		for i in range(numtime):
		        if ((i+1)%12 == 6):
		                plt.axhline(i, color = 'r')
		        if ((i+1)%12 == 8):
		                plt.axhline(i, color = 'k')

	loc = range(0, nlons, 30)
	lons = np.linspace(box[2], box[3], len(loc))
	lons = [int(round(l)) for l in lons]
	plt.xticks(loc, lons)

	lat1_sign = "N"
	lat2_sign = "N"
	if lats[0] < 0:
	        lat1_sign = "S"
	if lats[1] < 0:
	        lat2_sign = "S"

	plottitle = vname + "hovmollerh0" + "_" + str(abs(lats[0])) + lat1_sign + str(abs(lats[1])) + lat2_sign

	plt.title(plottitle)

	if showfig:
	        plt.show()
	if savefig:
	        plt.savefig(plottitle + ".ps")
	plt.clf()
	plt.cla()


