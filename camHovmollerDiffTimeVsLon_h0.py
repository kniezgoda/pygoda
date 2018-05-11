#!/Users/kyleniezgoda/anaconda/bin/python
import argparse
import numpy as np
from pygoda import findClimoFile, camgoda, h0dates, niceClev
import matplotlib.pyplot as plt
import sys

month_heating_on = 6
month_heating_off = 8
heating_left_lon = 0
heating_right_lon = 90

parser = argparse.ArgumentParser()
parser.add_argument('-cdir', '--control_directory', dest = 'control_directory', default = '.')
parser.add_argument('-tdir', '--test_directory', dest = 'test_directory', default = '.')
parser.add_argument('-grep', '--regular_expression', dest = 'regular_expression', default = '*cam.h0.')
parser.add_argument('-dates', '--date_bounds', dest = 'date_bounds', nargs = 2)
parser.add_argument('-lats', type = int, dest = 'lats', nargs = 2, default = (-5, 5))
parser.add_argument('-lons', type = int, dest = 'lons', nargs = 2, default = (0, 360))
parser.add_argument('-vmin', type = float, dest = 'vmin', nargs = 1, default = None)
parser.add_argument('-vmax', type = float, dest = 'vmax', nargs = 1, default = None)
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variable', nargs= "*", default = None)
parser.add_argument('-ft', '--file_type', dest = 'file_type', default = 'ps')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')
# parser.add_argument('-3d', '--var_is_3d', dest = 'var_is_3d', action = 'store_true')

# Handle the args
ARGS = parser.parse_args()
# var_is_3d = ARGS.var_is_3d
lats = ARGS.lats
lons = ARGS.lons
vmin = ARGS.vmin
vmax = ARGS.vmax
if vmin is not None:
	vmin = vmin[0]
if vmax is not None:
	vmax = vmax[0]
box = [lats[0], lats[1], lons[0], lons[1]]
dates = ARGS.date_bounds
start = dates[0]
end = dates[1]
control_directory = ARGS.control_directory
test_directory = ARGS.test_directory
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

	# Is the variable 2 or 3 dimensions
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
		c_nc_path, c_nc_file = findClimoFile(grep + date + "*", control_directory)
		t_nc_path, t_nc_file = findClimoFile(grep + date + "*", test_directory)
		#print c_nc_path
		#print t_nc_path
		c_nc = camgoda(c_nc_path)
		t_nc = camgoda(t_nc_path)
		
		# Extract the variable data
		# Special variables
		if var == "PRECT_d18O":
			t_nc.PRECT_d18O(box)
			c_nc.PRECT_d18O(box)
		
		elif var == "PRECT_dD":
			t_nc.PRECT_dD(box)
			c_nc.PRECT_dD(box)
		
		elif var == "PRECT_dxs":
			t_nc.PRECT_dxs(box)
			c_nc.PRECT_dxs(box)
		
		elif var == "QFLX_d18O":
			t_nc.QFLX_d18O(box)
			c_nc.QFLX_d18O(box)
		
		elif var == "QFLX_dD":
			t_nc.QFLX_dD(box)
			c_nc.QFLX_dD(box)
		
		elif var == "fluxDelta":
			t_nc.fluxDelta(box)
			c_nc.fluxDelta(box)
		
		elif var == "Column_d18OV":
			#test
			t_nc.variable('H2OV', box)
			denom = t_nc.columnSum(box)
			t_nc.variable('H218OV', box)
			num = t_nc.columnSum(box)
			t_nc.data = (num/denom - 1) * 1000
			#control
			c_nc.variable('H2OV', box)
			denom = c_nc.columnSum(box)
			c_nc.variable('H218OV', box)
			num = c_nc.columnSum(box)
			c_nc.data = (num/denom - 1) * 1000
		
		elif var == "Column_dDV":
			#test
			t_nc.variable('H2OV', box)
			denom = t_nc.columnSum(box)
			t_nc.variable('HDOV', box)
			num = t_nc.columnSum(box)
			t_nc.data = (num/denom - 1) * 1000
			#control
			c_nc.variable('H2OV', box)
			denom = c_nc.columnSum(box)
			c_nc.variable('HDOV', box)
			num = c_nc.columnSum(box)
			c_nc.data = (num/denom - 1) * 1000

		elif var == "P_E":
			#test
			t_nc.data = (t_nc.variable('PRECT', box, math = False)*1000 - t_nc.variable('QFLX', box, math = False)) * 60 * 60 * 24
			c_nc.data = (c_nc.variable('PRECT', box, math = False)*1000 - c_nc.variable('QFLX', box, math = False)) * 60 * 60 * 24
			t_nc.units = "kg/m2/day"
			t_nc.long_name = "Moisture flux"
			c_nc.units = "kg/m2/day"
			c_nc.long_name = "Advective moisture flux"

		elif var == "d18OV":
			t_nc.d18OV(box)
			c_nc.d18OV(box)
		elif var == "dDV":
			t_nc.dDV(box)
			c_nc.dDV(box)
		elif var == "dxsV":
			t_nc.dxsV(box)
			c_nc.dxsV(box)
		elif var == "psi":
			t_nc.psi(box)
			c_nc.psi(box)
		elif var == "RH":
			t_nc.RH(box)
			c_nc.RH(box)
		elif var == "VQ_d18O":
			t_nc.VQ_d18O(box)
			c_nc.VQ_d18O(box)
		elif var == "VQ_dD":
			t_nc.VQ_dD(box)
			c_nc.VQ_dD(box)
		elif var == "UQ_d18O":
			t_nc.UQ_d18O(box)
			c_nc.UQ_d18O(box)
		elif var == "UQ_dD":
			t_nc.UQ_dD(box)
			c_nc.UQ_dD( box)
		elif var == "QFLX_d18O":
			t_nc.QFLX_d18O(box)
			c_nc.QFLX_d18O(box)
		
		# Regular variables inside the netcdf file
		else:
			try:
				t_nc.variable(var, box, verb = True)
				c_nc.variable(var, box)
			except KeyError:
				print "Not able to plot variable " + var + "...\nSkipping this variable."
				print "Is this a 3-spatial-dimension variable? If so, append 3d_ to the beginning of the variable name."
				continue
		if var_is_3d:
			t_nc.data = t_nc.isobar(pressure)
			c_nc.data = c_nc.isobar(pressure)
		
		# Average over the latitudes
		c_hold = np.mean(c_nc.data, 0)
		t_hold = np.mean(t_nc.data, 0)
		
		# Take the diff
		hold = t_hold - c_hold
		
		# Set some extraneous variables if it's the first instance of the for loop
		if n == 0:
			units = c_nc.units
			master = hold
			nlons = len(c_nc.boxlon)
			boxlons = c_nc.boxlon
			nlats = len(c_nc.boxlat)
		# Otherwise stack the data together
		else:
			master = np.vstack((master, hold))

	clev = niceClev(master)
	if vmin is None:
		vmin = clev[0]
	if vmax is None:
		vmax = clev[-1]

	master_plot = plt.imshow(master, cmap = plt.cm.RdBu_r, vmin = vmin, vmax = vmax, aspect = 'auto', interpolation = 'none')

	plt.ylabel('months since ' + start)
	plt.xlabel('longitude')

	cbar = plt.colorbar(master_plot, orientation = 'horizontal')
	cbar.set_label(units)

	numtime = master.shape[0] 

	# for i in range(numtime):
	# 	if ((i+1)%12 == 6):
	# 		plt.axhline(i, color = 'r')
	# 	if ((i+1)%12 == 8):
	# 		plt.axhline(i, color = 'k')	

	loc = range(0, nlons, 20)
	loc = loc[0:]
	lons = [boxlons[b] for b in loc]
	plt.xticks(loc, lons)

	loc = range(0, numtime)
	months = [str(x) for x in range(1,numtime/30)]
	#plt.yticks(loc, months)

	lat1_sign = "N"
	lat2_sign = "N"
	if lats[0] < 0:
		lat1_sign = "S"
	if lats[1] < 0:
		lat2_sign = "S"
	plottitle = vname + "hovmollerDiffh0" + "_" + str(abs(lats[0])) + lat1_sign + str((lats[1])) + lat2_sign 

	plt.title(plottitle)

	if showfig:
		plt.show()
	if savefig:
		plt.savefig(plottitle + ".ps")
	plt.clf()
	plt.cla()
