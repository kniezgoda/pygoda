from pygoda import ncgoda
from glob import glob
import numpy as np
import os

# Set the variables you want to use
mh_path = '/home/server/student/homes/kniezgod/model_runs/F.C5.2deg.wiso.obs6kSST_kn003'
pi_path = '/home/server/student/homes/kniezgod/model_runs/F.C5.2deg.wiso.defaultSSTICE_kn002'
regexp = "*cam.h0.00*"
var = "PRECT"
region = (10, 30, 0, 90)
# Set the way the dates are spelled in the netcdf file name
# Only need to worry about parts of the netcdf file name with dashes
# datetype = "yyyy-mm-dd-ssssss"
datetype = "yyyy-mm"


#------------------#
# Code starts here #
 #------------------#

# Creates a dict for the months to add data to
mh_var_dict = {"01" : [], "02" : [], "03" : [], "04" : [], "05" : [], "06" : [], \
		       "07" : [], "08" : [], "09" : [], "10" : [], "11" : [], "12" : [], }

pi_var_dict = {"01" : [], "02" : [], "03" : [], "04" : [], "05" : [], "06" : [], \
		       "07" : [], "08" : [], "09" : [], "10" : [], "11" : [], "12" : [], }

# Calculate for MH
os.chdir(mh_path)
counter = 0
h1files = glob(regexp)
for f in h1files:
	print counter
	print f
	# Find the month
	if datetype == "yyyy-mm-dd-ssssss":
		mon = f.split("-")[-3]
	#
	if datetype == "yyyy-mm":
		mon = f.split("-")[-1][:2]
	#
	print "Month is " + mon
	# Load the file
	d = ncgoda(f)
	# Extract the variable data
	data = d.variable(var, box = region)
	# Add the average to the correct list in the master dict
	mh_var_dict[mon].append(np.mean(data))
	counter += 1

# Calculate for PI
os.chdir(pi_path)
counter = 0
h1files = glob(regexp)
for f in h1files:
	print counter
	print f
	# Find the month
	if datetype == "yyyy-mm-dd-ssssss":
		mon = f.split("-")[-3]
	#
	if datetype == "yyyy-mm":
		mon = f.split("-")[-1][:2]
	#
	print "Month is " + mon
	# Load the file
	d = ncgoda(f)
	# Extract the variable data
	data = d.variable(var, box = region)
	# Add the average to the correct list in the master dict
	pi_var_dict[mon].append(np.mean(data))
	counter += 1

diff_var_dict = {"01" : [], "02" : [], "03" : [], "04" : [], "05" : [], "06" : [], \
		         "07" : [], "08" : [], "09" : [], "10" : [], "11" : [], "12" : [], }

for m in diff_var_dict.keys():
	diff_var_dict[m] = np.mean(mh_var_dict[m]) - np.mean(pi_var_dict[m])