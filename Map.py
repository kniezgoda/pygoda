#!/glade/u/apps/ch/opt/python/2.7.13/gnu/6.2.0/bin/python
'''
Creates one map of a variable

Author: Kyle Niezgoda
'''

import os, sys, glob, argparse
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap as bm 
from pygoda import camgoda, findClimoFile, niceClev, RegularClev
root = os.getcwd()
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
Probably will want to put the following variables in a text file
or some file structure and then have the user supply the file in a 
command line flag or argument.
'''

#------------- Custom subplot title names 
# these are the names of the subplots, not the main figure --- the main figure always has the same name (the variable long name from the netcdf file)
# Set to None (no quotes) for default subplot names, the name of the file used
# Default 
user_control_plot_title = None
# Regular 
#user_control_plot_title = "Control: Pre-industrial (1850)"
#user_control_plot_title = "Control: Simulated MH"


##############################################
###### Don't change anything below here ######
##############################################


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#################################################
####### Parse arguments from command line #######
#################################################

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
parser.add_argument('-f', '--files', dest = 'files', nargs = "*")
parser.add_argument('-r', '--region', dest = 'region', default = 'GlobalTropics')
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-v', '--variable', dest = 'variable', nargs= "*", default = None)
parser.add_argument('-ft', '--file_type', dest = 'file_type', default = 'ps')
parser.add_argument('-clev', dest = 'clev', type = float, nargs = 3, default = None)
parser.add_argument('-barbs', '--wind_barb_pressure', dest = 'wind_barb_pressure', nargs = 1, type = float, default = None)
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')
parser.add_argument('-verb', '--verbose', dest = 'verbose', action = 'store_true')

ARGS = parser.parse_args()
clev = ARGS.clev
region = ARGS.region
verbose = ARGS.verbose

barbs = ARGS.wind_barb_pressure
show_barbs = False
if barbs is not None:
	barb_pressure = ARGS.wind_barb_pressure[0] * 100
	show_barbs = True

if verbose:
	print "Region is " + region
files = ARGS.files
savefig = ARGS.savefig
showfig = ARGS.showfig
variable = ARGS.variable
ftype = ARGS.file_type
mkdir = True
if ARGS.developer_mode:
	if verbose:
		print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
	savefig = False
	showfig = True
	mkdir = False


# Set the lat bounds
# Default global tropics
region_name = "GlobalTropics"
southern_lat = -85
northern_lat = 85
left_lon = 0
right_lon = 355

# Indian monsoon
if (region == "IM") | (region == "IndianMonsoon"):
	region_name = "IndianMonsoon"
	southern_lat = -5
	northern_lat = 45
	left_lon = 40
	right_lon = 110

# Greater maritime continent
if (region == "MC") | (region == "MaritimeContinent"):
	region_name = "MaritimeContinent"
	southern_lat = -20
	northern_lat = 20
	left_lon = 80
	right_lon = 160

# North Africa
if (region == "NA") | (region == "NorthAfrica"):
	region_name = "NorthAfrica"
	southern_lat = -20
	northern_lat = 45
	left_lon = -30
	right_lon = 70

# Rodwell and Hoskins region
if (region == "RH") | (region == "RodwellHoskins"):
	region_name = "RodwellHoskins"
	southern_lat = 15
	northern_lat = 55
	left_lon = -10
	right_lon = 95

#Central Equatorial Pacific to Western eq. Atlantic
if (region == "EP") | (region == "TropicalOceans"):
	region_name = "TropicalOceans"
	southern_lat = -25
	northern_lat = 25
	left_lon = 180
	right_lon = 355

#South America and the Carribean 
if (region == "SA") | (region == "SouthAmerica"):
	region_name = "SouthAmerica"
	southern_lat = -30
	northern_lat = 30
	left_lon = -120
	right_lon = -20

# Four proxy region in the warm pool
if (region == "4P") | (region == "FourProxies"):
	region_name = "FourProxies"
	southern_lat = -25
	northern_lat = 30
	left_lon = 85
	right_lon = 160

if mkdir:
	# Create Map directory is it doesn't exist
	if not os.path.exists("Map"):
		os.mkdir("Map")
		if verbose:
			print "Created directory " + "Map"

	# Create the region directory if it doesn't already exist
	if not os.path.exists("Map/" + region_name):
		os.mkdir("Map/" + region_name)
		if verbose:
			print "Created directory " + "Map/" + region_name

for filein in files:
	print "Current files is : " + filein + '\n'
	# Read the data
	ncdata = camgoda(filein)

	# Set the boxlat and boxlon
	box = (southern_lat, northern_lat, left_lon, right_lon)
	ncdata.variable(ncdata.vars[0], box, setData = False) # this sets self.boxlon and self.boxlat

	# Create bm coords from region bounds
	# bm lonitude coords need to be 0 < coord < 360 
	bmlon, bmlat = np.meshgrid(ncdata.boxlon, ncdata.boxlat)

	# Reset the lat and lon bounds so that Map don't show grey areas 
	# southern_lat, northern_lat = np.array(ncdata.boxlat)[[0,-1]]
	# Change lons to be negative is 180 < lon < 360 because that's how bm works for 'cea' projection
	# left_lon, right_lon = np.array(ncdata.boxlon)[[0,-1]]
	if 0 in ncdata.boxlon[1:-2]: # if we cross the gml
		left_lon = ncdata.boxlon[0]-360


	g = -9.8 # gravitational constant

	# Change into figure directory (root/region/season/) for image creation
	if mkdir:
		os.chdir("Map/" + region_name)

	#------------------#
	# Read in the data #
	#------------------#

	for V in variable:
		# Extract data for wind barbs if needed
		if show_barbs:
			nc_v = ncdata.variable("V", box)
			nc_v = ncdata.isobar(barb_pressure)
			nc_u = ncdata.variable("U", box)
			nc_u = ncdata.isobar(barb_pressure)
			
			density = 5
			colsby4 = np.arange(0, bmlon.shape[1], density)
			rowsby4 = np.arange(0, bmlat.shape[0], density)
			barb_lon = bmlon[:,colsby4]
			barb_lon = barb_lon[rowsby4,:]
			barb_lat = bmlat[:,colsby4]
			barb_lat = barb_lat[rowsby4,:]
				
			nc_v = nc_v[rowsby4,:]	
			nc_v = nc_v[:,colsby4]
			nc_u = nc_u[rowsby4,:]	
			nc_u = nc_u[:,colsby4]

		# Extract variable info (sets var, vname, and pressure)
		var_is_3d = False
		if V[:3] == '3d_':
			var_is_3d = True
			
		if not var_is_3d:
			var = V
			vname = V
			pressure = None
		else:
			var = V[3:-3]
			vname = V[3:]
			pressure = int(V[-3:]) * 100

		print("\nPlotting " + vname  +  " data...")

		# Extract the variable data
		# Special variables
		if var == "PRECT_d18O":
			ncdata.PRECT_d18O(box)
		
		elif var == "PRECT_dD":
			ncdata.PRECT_dD(box)
		
		elif var == "PRECT_dxs":
			ncdata.PRECT_dxs(box)
		
		elif var == "QFLX_d18O":
			ncdata.QFLX_d18O(box)
		
		elif var == "QFLX_dD":
			ncdata.QFLX_dD(box)
		
		elif var == "fluxDelta":
			ncdata.fluxDelta(box)
		
		elif var == "Column_d18OV":
			ncdata.variable('H2OV', box)
			denom = ncdata.columnSum(box)
			ncdata.variable('H218OV', box)
			num = ncdata.columnSum(box)
			ncdata.data = (num/denom - 1) * 1000
		
		elif var == "Column_dDV":
			ncdata.variable('H2OV', box)
			denom = ncdata.columnSum(box)
			ncdata.variable('HDOV', box)
			num = ncdata.columnSum(box)
			ncdata.data = (num/denom - 1) * 1000

		elif var == "P_E":
			ncdata.data = (ncdata.variable('PRECT', box, math = False)*1000 - ncdata.variable('QFLX', box, math = False)) * 60 * 60 * 24
			ncdata.units = "kg/m2/day"
			ncdata.long_name = "Advective moisture flux"

		elif var == "d18OV":
			ncdata.d18OV(box)
		elif var == "dDV":
			ncdata.dDV(box)
		elif var == "dxsV":
			ncdata.dxsV(box)
		elif var == "psi":
			ncdata.psi(box)
		elif var == "RH":
			ncdata.RH(box)
		elif var == "VQ_d18O":
			ncdata.VQ_d18O(box)
		elif var == "VQ_dD":
			ncdata.VQ_dD(box)
		elif var == "UQ_d18O":
			ncdata.UQ_d18O(box)
		elif var == "UQ_dD":
			ncdata.UQ_dD( box)
		elif var == "QFLX_d18O":
			ncdata.QFLX_d18O(box)
		
		# Regular variables inside the netcdf file
		else:
			try:
				ncdata.variable(var, box)
			except KeyError:
				print "Not able to plot variable " + var + "...\nSkipping this variable."
				print "Is this a 3-spatial-dimension variable? If so, append 3d_ to the beginning of the variable name."
				continue
		if var_is_3d:
			ncdata.data = ncdata.isobar(pressure)

		#----------------#
		# Create the map #
		#----------------#

		fig = plt.figure()

		# Set up some map properties
		ncdata.clevs = RegularClev(ncdata.data)
		ncdata.prep_map('ANN', region) # This sets the cmap 

		if clev is not None:
			ncdata.clevs = np.linspace(clev[0], clev[1], clev[2])

		m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
		m.drawcoastlines()
		m.drawmapboundary(fill_color='0.3')
		print ncdata.data
		cs = m.contourf(bmlon, bmlat, ncdata.data, ncdata.clevs, shading = 'flat', latlon = True, cmap=ncdata.cmap)
		cbar = m.colorbar(cs, location='right', pad="5%")
		cbar.set_label(ncdata.units, fontsize = 8)
		if user_control_plot_title is None:
			control_title = vname + " from " + filein
		else:
			control_title = user_control_plot_title 
		plt.title(control_title, fontsize = 8)
			
		# Make things pretty
		if savefig:
			plt.savefig(vname + "_" + filein + "." + ftype, bbox_inches='tight', dpi = 500)
			if verbose:
				print("Created " + vname + "-" + filein + "." + ftype + '\n')
		if showfig:
			plt.show()

		plt.clf()
		plt.cla()
		plt.close(fig)

		os.chdir(root)

sys.exit()
