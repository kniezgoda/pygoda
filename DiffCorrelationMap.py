#! /home/server/student/homes/kniezgod/.conda/envs/condagoda/bin/python
from pygoda import camgoda, camdates, findClimoFile
import numpy as np
from scipy.stats import pearsonr as regress
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as bm
import argparse


'''
Author: Kyle Niezgoda, April 28, 2019

This code produces maps of correlation arrays between the local variable (-lv) at a point (-latlon) 
and the field variable (-fv) over the globe for a specified time (-dates)
'''
#########################
# Read in the arguments #
#########################

parser = argparse.ArgumentParser()
parser.add_argument('-cdir', '--control_directory', dest = 'cdir', nargs = 1, default = 'fc5.2deg.wiso.piControl_kn028/atm/hist')
parser.add_argument('-tdir', '--test_directory', dest = 'tdir', nargs = 1, default = 'fc5.2deg.wiso.mh6ka_kn032/atm/hist')
parser.add_argument('-loc_var', '--local_variable', dest = 'lv', nargs = 1, default = None)
parser.add_argument('-field_var', '--field_variable', dest = 'fv', nargs = 1, default = None)
parser.add_argument('-latlon', dest = 'latlon', nargs = 2, default = (0,0))
parser.add_argument('-del', dest = 'delta', nargs = 1, default = [0])
parser.add_argument('-years', dest = 'years', nargs = 2, default = [10,30])
parser.add_argument('-months', dest = 'months', nargs = '*', default = [1,2,3,4,5,6,7,8,9,10,11,12])
parser.add_argument('-grep_pre', dest = 'grep_pre', nargs = 1, default = [''])
parser.add_argument('-grep_post', dest = 'grep_post', nargs = 1, default = [''])
parser.add_argument('-alpha', '--alpha', dest = 'alpha', nargs = 1, default = .05)
parser.add_argument('-stiple', '--stiple', dest = 'stiple', action = "store_true")
parser.add_argument('-nosave', '--dont_save_figure', dest = 'savefig', action = 'store_false')
parser.add_argument('-show', '--showfig', dest = 'showfig', action = 'store_true')
parser.add_argument('-dev', '--developer_mode', dest = 'developer_mode', action = 'store_true')

ARGS = parser.parse_args()
cdir = ARGS.cdir
tdir = ARGS.tdir
lv = ARGS.lv[0]
fv = ARGS.fv[0]
lat, lon = [int(x) for x in ARGS.latlon]
delta = int(ARGS.delta[0])
start, end = [int(x) for x in ARGS.years]
months = [int(m) for m in ARGS.months]
grep_pre = ARGS.grep_pre[0]
grep_post = ARGS.grep_post[0]
savefig = ARGS.savefig
showfig = ARGS.showfig
alpha = float(ARGS.alpha)
stiple = ARGS.stiple
if ARGS.developer_mode:
    print "\nRunning in dev mode. No files will be saved, no directories will be created, and all plots will be printed to the screen."
    savfig = False
    showfig = True


#############
# Main code #
#############

# Make the date array
dates = camdates(start, end, months)

# Read in each file, extract the local variable at the point and field variable on the globe, track them in a master array
# Correlate the local variable 1-d array to each grid cell of the field variable array
# Return the correlation coefficients 

bottom = lat-delta
top = lat+delta
left = lon-delta
if left < 0:
	left += 360
right = lon+delta
if right > 360:
	right -= 360
	
region = [bottom, top, left, right]
tloc_var_master = []
cloc_var_master = []
for n, d in enumerate(dates):
	# Find the file
	cpath, cfilename = findClimoFile(cdir + "/" + grep_pre+'*'+d+'*'+grep_post)
	tpath, tfilename = findClimoFile(tdir + "/" + grep_pre+'*'+d+'*'+grep_post)
	print tfilename
	print cfilename
	
	# Open the file
	cnc = camgoda(cpath)
	tnc = camgoda(tpath)

	# Read the local and field variables in
	# each = 0 ---> lv
	# each = 1 ---> fc
	for each in range(2):
		if each == 0:
			V = lv
			box = region
		else:
			V = fv
			box = None

		tnc.ExtractData(V, box)
		cnc.ExtractData(V, box)
		
		if each == 0:
			tloc_var = tnc.data
			cloc_var = cnc.data
		else:
			tfield_var = tnc.data
			cfield_var = cnc.data

	tloc_var = np.mean(tloc_var)
	# Add lv to the master list
	tloc_var_master.append(tloc_var)
	
	cloc_var = np.mean(cloc_var)
	cloc_var_master.append(cloc_var)
	    
	if n == 0:
		# Create the master array for the field variable
		tfield_var_master = np.zeros((len(dates),) + tfield_var.shape)
		cfield_var_master = np.zeros((len(dates),) + cfield_var.shape)
		lats = cnc.lat
		lons = cnc.lon
	# Add fv to the master array
	tfield_var_master[n,:,:] = tfield_var
	cfield_var_master[n,:,:] = cfield_var


# Subtract the test from control for loc var and field var to get the diff arrays
loc_var_master = np.subtract(tloc_var_master, cloc_var_master)
field_var_master = np.subtract(tfield_var_master, cfield_var_master)

# We now have the 3-d field variable array and the 1-d local variable array
# Now correlate the local variable array to each grid cell of the field variable (fielf_var_master[:,i,j])
# and save the correlation coefficient to a 2-d array
dim1 = field_var_master.shape[1]
dim2 = field_var_master.shape[2]
corr_array = np.zeros((dim1, dim2))
for i in range(dim1):
	for j in range(dim2):
		r, pval = regress(loc_var_master, field_var_master[:,i,j])
		if stiple:
			if pval > alpha:
				corr_array[i,j] = np.nan
			else:
				corr_array[i,j] = r
		else:
			corr_array[i,j] = r

# Make the map of the corr array
bmlon, bmlat = np.meshgrid(lons, lats)
m = bm(projection = 'cea', llcrnrlat=-50,urcrnrlat=50, llcrnrlon=0,urcrnrlon=360,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
clevs = np.linspace(-1, 1, 21)
cs = m.contourf(bmlon, bmlat, corr_array, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
cbar.set_label("correlation-coefficient", fontsize = 8)

# Put a point on the map to show the location
x,y = m(lon,lat)
m.plot(x, y, 'gx')

# Put some words on it
if bottom == top:
	if left == right:
		title = "Correlation between\n" + lv + " at lat=" + str(bottom) + ", lon=" + str(right) + "\nand global " + fv
	else:
		title = "Correlation between\n" + lv + " at lats=" + str(top) + ", lons=" + str(left) + "-" + str(right) + "\nand global " + fv
else:
	if left == right:
		title = "Correlation between\n" + lv + " at lats=" + str(bottom) + "-" + str(top) + ", lon=" + str(right) + "\nand global " + fv
	else:
		title = "Correlation between\n" + lv + " at lats=" + str(bottom) + "-" + str(top) + ", lons=" + str(left) + "-" + str(right) + "\nand global " + fv
		
plt.title(title)
if savefig:
	plt.savefig("correlationMap.pdf")
if showfig:
	plt.show()

