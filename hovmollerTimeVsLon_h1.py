#!/glade/u/apps/ch/opt/python/2.7.13/gnu/6.2.0/bin/python
import argparse
import numpy as np
from pygoda import findClimoFile, ncgoda, h1dates, niceClev
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument('-dir', '--directory', dest = 'directory', default = '.')
parser.add_argument('-grep', '--regular_expression', dest = 'regular_expression', default = '*cam.h1.')
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
date_list = h1dates(start, end)
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
		nc = ncgoda(nc_path)
		
		# Read in the variable data
		if var == "PRECT_d18O":
			nc.PRECT_d18O(box)
		elif var == "PRECT_dD":
			nc.PRECT_dD(box)		
		elif var == "PRECT_dxs":
			nc.PRECT_dxs(box)
		elif var == "QFLX_d18O":
			nc.QFLX_d18O(box)		
		elif var == "QFLX_dD":
			nc.QFLX_dD(box)		
		elif var == "fluxDelta":
			nc.fluxDelta(box)
		elif var == "Column_d18OV":
			nc.variable('H2OV', box)
			denom = nc.columnSum(box)
			nc.variable('H218OV', box)
			num = nc.columnSum(box)
			nc.data = (num/denom - 1) * 1000	
		elif var == "Column_dDV":
			nc.variable('H2OV', box)
			denom = nc.columnSum(box)
			nc.variable('HDOV', box)
			num = nc.columnSum(box)
			nc.data = (num/denom - 1) * 1000
		elif var == "P_E":
			nc.data = (nc.variable('PRECT', box, math = False)*1000 - nc.variable('QFLX', box, math = False)) * 60 * 60 * 24
			nc.units = "kg/m2/day"
			nc.long_name = "Advective moisture flux"
		elif var == "d18OV":
			nc.d18OV(box)
		elif var == "dDV":
			nc.dDV(box)
		elif var == "dxsV":
			nc.dxsV(box)
		elif var == "psi":
			nc.psi(box)
		elif var == "RH":
			nc.RH(box)
		elif var == "VQ_d18O":
			nc.VQ_d18O(box)
		elif var == "VQ_dD":
			nc.VQ_dD(box)
		elif var == "UQ_d18O":
			nc.UQ_d18O(box)
		elif var == "UQ_dD":
			nc.UQ_dD( box)
		elif var == "QFLX_d18O":
			nc.QFLX_d18O(box)
		
		# Regular variables inside the netcdf file
		else:
			try:
				nc.variable(var, box)
			except KeyError:
				print "Not able to plot variable " + var + "...\nSkipping this variable."
				print "Is this a 3-spatial-dimension variable? If so, append 3d_ to the beginning of the variable name."
				continue
		if var_is_3d:
			nc.data = nc.isobar(pressure)

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

	clev = niceClev(master, alpha = .9)
	try:
		vmin = clev[0]
		vmax = clev[-1]
	except TypeError:
		vmin = None
		vmax = None

	master_plot = plt.imshow(master, cmap = plt.cm.RdBu_r, vmin = vmin, vmax = vmax, aspect = 'auto')

	plt.ylabel('months since ' + start)
	plt.xlabel('longitude')

	cbar = plt.colorbar(master_plot, orientation = 'horizontal')
	cbar.set_label(units)

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

	plottitle = vname + "hovmollerh1" + "_" + str(abs(lats[0])) + lat1_sign + str(abs(lats[1])) + lat2_sign

	plt.title(plottitle)

	if showfig:
	        plt.show()
	if savefig:
	        plt.savefig(plottitle + ".ps")
	plt.clf()
	plt.cla()



