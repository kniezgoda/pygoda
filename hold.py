# A collection of functions useful for manipulating and working with netcdf files in python


# =========================================================================================== #


####################
### find_indices ### 
####################

def find_indices(box, lats, lons):
 	'''
	Author : Kyle Niezgoda
	
	This function is useful for finding the indices of lat and lon lists
	that correspond to the latitude and longitudes in a box.
	
	For instance, maybe you want to extract data for India only, or within
	the lat range of 5N to 35N and lon range of 65E to 90E. To do this, you need
	to know the indices in the lat and lon list that are closest to the India 
	corrdinates. 
	
	To find these indices, you would do the following:
	
	lats = f.variables['lat'][:]
	lons = f.variables['lons'][:]
	region = [5, 35, 65, 90]
	bottom_index, top_index, left_index, right_index = find_indices(region, lats, lons)
	
	Arguments
	---------
	box : list, [bottom_edge, top_edge, left_edge, right_edge] all in degrees
	lats : list of latitudes from the netcdf file
	lons : list of longitudes from the netcdf file
	
	Returns
	-------
	list, [bottom_index, top_index, left_index, right_index]
	'''
	import numpy as np

	if box is None:
		return [np.arange(len(lats)), np.arange(len(lons))]

	box = list(box)

	if box[2] < 0:
		box[2] = 360 + box[2]

	if box[3] < 0:
		box[3] = 360 + box[3]	

	def cat_idx():
		if l > r:
			return np.concatenate((np.arange(l, len(lons)), np.arange(0,r)))
		elif l == r:
			return [l]
		else:
			return np.arange(l,r)

	b = abs(lats - box[0]).tolist().index(min(abs(lats-box[0])))
	t = abs(lats - box[1]).tolist().index(min(abs(lats-box[1])))
	l = abs(lons-(box[2])).tolist().index(min(abs(lons-(box[2]))))
	r = abs(lons-(box[3])).tolist().index(min(abs(lons-(box[3]))))

	if b == t:
		return [[b], cat_idx()]

	return [np.arange(b,t), cat_idx()]



# =========================================================================================== #


##############
### ncdump ###
##############

def ncdump(nc_Dataset, var='total', verb=True):
	'''
	Copied and modified from http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html

	ncdump outputs dimensions, variables and their attribute information.
	The information is similar to that of NCAR's ncdump utility.
	ncdump requires a valid instance of Dataset.

	Arguments
	---------
	nc_Dataset : netCDF4.Dataset
	    A netCDF4 dateset object
	var : Name of variable (str)
		Default is 'total' for all the variable info
	verb : Verbose (boolean)

	Returns
	-------
	nc_dims : list
	    A Python list of the NetCDF file dimensions
	nc_vars : list
	    A Python list of the NetCDF file variables
	'''
	from netCDF4 import Dataset
	def print_ncattr(key):
	    """
	    Prints the NetCDF file attributes for a given key

	    Parameters
	    ----------
	    key : unicode
	        a valid netCDF4.Dataset.variables key
	    """
	    try:
	        print "\t\ttype:", repr(nc_Dataset.variables[key].dtype)
	        for ncattr in nc_Dataset.variables[key].ncattrs():
	            print '\t\t%s:' % ncattr,\
	                  repr(nc_Dataset.variables[key].getncattr(ncattr))
	    except KeyError:
	        print "\t\tWARNING: %s does not contain variable attributes" % key

	nc_dims = [dim for dim in nc_Dataset.dimensions]  # list of nc dimensions  
	# Dimension information.
	if verb:
	    print "NetCDF dimension information:"
	    for dim in nc_dims:
	        print "\tName:", dim, "\tsize:", len(nc_Dataset.dimensions[dim])

	print "\n"

	# Variable information.
	if var == 'total':
	    nc_vars = [v for v in nc_Dataset.variables]  # list of all nc variables
	else:
	    if type(var) is list: #list of variables
	        nc_vars = var
	    else: 
	        nc_vars = var.split("!@#$%^&*") #single variable
	if verb:
	    print "NetCDF variable information:"
	    for v in nc_vars:
	        if v not in nc_dims:
	            print '\tName:', v
	            print "\t\tdimensions:", nc_Dataset.variables[v].dimensions
	            print "\t\tsize:", nc_Dataset.variables[v].size
	            print_ncattr(v)

	return nc_dims, nc_vars


# =========================================================================================== #


########################
### CorrelationArray ###
########################

def CorrelationArray(files, field_var, loc_var, lat, lon, loc_var_lev = None, field_var_lev = None):
	'''
	Author : Kyle Niezgoda
	Date : June 2, 2017

	This function creates an array of correlation coefficients
	The coefficients are calculated between a variable at a 
	specific location and a different variable over the entire 
	spatial field. For instance, to correlate the d18O of precip
	in Singapore to the OLR at each lat/lon point in the world, 
	you would execute:
		CorrelationArray(files, "OLR", "PRECT_d18O", 1.3, 103)
	where files are the datafiles for the timeseries over which
	you wish to calculate the correlation.

	Arguments
	---------
	files : list or str
		-if str, use wildcards to search for multiple files 
		e.g.: *cam.h1* will find all cam h1 files
		-if list, enter each file name as an index of the list
	lat : numerical within the range [-90, 90]
	lon : numerical within the range [0, 360] w/ 0 = GMT, prime meridian 
	field_var : the variable to extract at every lat/lon (str)
	loc_var : the variable to extract only at the lat and lon location (str)
	loc_var_lev (optional) : the level of the variable, if the variable is 3d
	field_var_lev (optional) : same as loc_var_lev

	Returns
	-------
	array : correlation array (2d)
	'''
	import os, glob
	from netCDF4 import Dataset
	from UsefulNetcdfFunctions import find_indices
	import numpy as np
	region = [lat, lat, lon, lon]

	# Find all the cam.h1 files
	flist = files
	if type(files) is str:
		flist = glob.glob(files)
	print "Files found: \n" + str(flist)	
	loc_master = []
	field_master = []
	loc_coords = [lat, lat, lon, lon]
	for n, f in enumerate(flist):
		data = Dataset(f, mode = 'r')

		# Do the location stuff first
		# Find the location coordinate indices
		lats = data.variables['lat'][:]
		lons = data.variables['lon'][:]
		loc_lat = find_indices(region, lats, lons)[0]
		loc_lon = find_indices(region, lats, lons)[2]

		# Extract the location variable data and add it to the master list
		if loc_var_lev is not None:
			loc_master.append(data.variables[loc_var][:, lov_var_lev, loc_lat, loc_lon].squeeze().tolist())
		if loc_var_lev is None: 
			loc_master.append(data.variables[loc_var][:, loc_lat, loc_lon].squeeze().tolist())
		
		# Extract all the field data and keep track of the data
		if n == 0: # Initialize
			if field_var_lev is not None:
				field_master = data.variables[field_var][:, field_var_lev, :, :].squeeze()
			if field_var_lev is None:
				field_master = data.variables[field_var][:].squeeze()
		else: # Stack
			if field_var_lev is not None:
				field_master = np.dstack((field_master, data.variables[field_var][:, field_var_lev, :, :].squeeze()))
			if field_var_lev is None:
				field_master = np.dstack((field_master, data.variables[field_var][:].squeeze()))
			
	
	# Now we have a the 1-d list of the loc_var and 3-d list of that field_var
	# The lengths of both of these lists should correspond to the number of files we read in
	
	# Loop through lats and lons and run the correlation
	num_ind1 = field_master.shape[0]
	num_ind2 = field_master.shape[1]
	corr_array = np.empty([num_ind1, num_ind2])
	for i in range(num_ind1):
		for j in range(num_ind2):
			corr_array[i,j] = np.corrcoef(loc_master, field_master[i,j,:])[0,1]

	return [corr_array, lats, lons]


# =========================================================================================== #


############
### base ###
############

def base(Dataset, idxbox = None):
	'''
	returns lat, lon, A, B, and PS from a netCDF4 Dataset
	'''
	lat = Dataset.variables['lat'][:].squeeze()
	lon = Dataset.variables['lon'][:].squeeze()
	PS = Dataset.variables['PS'][:].squeeze()
	A = Dataset.variables['hyam'][:].squeeze()
	B = Dataset.variables['hybm'][:].squeeze()
	if idxbox is not None:
		lat = lat[idxbox[0]:idxbox[1]]
		lon = lon[idxbox[2]:idxbox[3]]
		PS = PS[idxbox[0]:idxbox[1],idxbox[2]:idxbox[3]]
	return lat, lon, A, B, PS


# =========================================================================================== #


####################
### PressureCalc ###
####################

def PressureCalc(A, B, PS):
	'''
	Author : Kyle Niezgoda
	Date : June 9, 2017

	This function creates a 3-d array of pressure values.
	For CESM, the 3-d pressure variable created and printed 
	out in the history files is not the true pressure. You 
	have to apply the coefficients A and B to the values, 
	along with the surface pressure, to get the real pressure
	since CESM has terrain following pycnals (or something
	like that.) The equation is P = A + B * PS w/ PS = surface
	pressure.

	Arguments
	---------
	A : list of length n where n is the number of vertical levels
	B : same as A
	PS : 2-d array of surface pressure values

	Returns
	-------
	array : 3-d array of pressure values with dimensions =
	[n, PS.shape[0], PS.shape[1]] (lev, lat, lon)
	'''
	import numpy as np
	pressure = []
	for n in range(len(B)):
		pressure.append(np.array(A[n] * 100000 + B[n] * PS))
	return np.array(pressure)


# =========================================================================================== #


#############
### dR_dt ###
#############

def dR_dt(f, isotope, lev = -1, time = 0):
	'''
	Author : Kyle Niezgoda
	Date : June 9, 2017

	This function calculates dR/dt due to horizontal advection 
	where R is the ratio of isotope to normal water in vapor form.

	Arguments
	---------
	file : string, the file name to calculate dR/dt from
	isotope : string, which isotope to calculate, either HDO or H218O
	lev : int, the level to caulate dR/dt at
		default : -1 (surface)
	time : int, if the file has more than one time file, which one to use?
		default : 0

	Returns
	-------
	array : 2-d array of dR/dt values in units of R per second
	'''
	from netCDF4 import Dataset
	import numpy as np

	Rad = 6371000 # meters, radius of Earth

	# Read the file
	nc = Dataset(f, mode = 'r')

	# Import and calculate some spatial variables
	lat = nc.variables['lat'][:]
	# Only works if longitudes and latitudes are equally spaced on a finite volume grid!
	dy = 2 * np.pi * Rad * np.diff(lat)[0] / 360
	lon = nc.variables['lon'][:]
	# Only works if longitudes and latitudes are equally spaced on a finite volume grid!
	delta_lon = np.diff(lon)[0]
	dx = 2 * np.pi / 360 * Rad * delta_lon * np.cos(np.radians(lat))

	# Figure out what the user meant with the isotope argument
	if (isotope == "hdo") | (isotope == "HDO"):
		which = "HDO"
		#print "Extracting data for " + which
	elif (isotope == "h218o") | (isotope == "H218O") | (isotope == "H218o") | (isotope == "h2180") | (isotope == "H2180"):
		which = "H218O"
		#print "Extracting data for " + which
	else:
		print "Error from dR_dt.py..."
		print "Do not understand the isotope argument or no isotope argument entered!"
		print "Fix the isotope argument, make sure it is a string, a try again."
		print "Exiting program..."
		return

	# Read in the important stuff
	UQ = nc.variables['UQ'][time,lev,:,:].squeeze()
	UQ_HDO = nc.variables['UQ_'+which][time,lev,:,:].squeeze()
	HDOV = nc.variables[which+'V'][time,lev,:,:].squeeze()
	Q = nc.variables['Q'][time,lev,:,:].squeeze()
	R = HDOV/Q

	# Initialize the spatial gradient master arrays
	dUQ_dx = np.zeros(UQ.shape)
	dUQ_dy = np.zeros(UQ.shape)
	dUQHDO_dx = np.zeros(UQ_HDO.shape)
	dUQHDO_dy = np.zeros(UQ_HDO.shape)

	######################
	### Plan of attack ###
	######################
	'''
	We need to calculate d/dt(UQ) and d/dt(UQ_HDO)
	To do so we'll need to find spatial gradients for
	UQ_HDO and UQ in both x and y. We will ignore 
	gradients in the z. 

	We will set d/dy = 0 for the top and bottom lat
	as a boundary condition for both variables.

	Algorithm:
	# 1) loop through rows and columns
	# 2) calculate the dx one by one inside each column loop
	# 3) calculate all the dy at once inside the first column loop (j==0)
	# 4) dR/dx = -d/dx(UQ_HDO) + R*d/dx(UQ), same for dy
	# 5) dR/dt = (dR/dx + dR/dy) / q
	'''

	# Loop through rows
	for i in range(UQ.shape[0]):
		# Extract data for this iteration
		uqrow = UQ[i,:]
		uqhdorow = UQ_HDO[i,:]
		
		# Initialize dx carriers
		duq_dx = []
		duqhdo_dx = []

		# Loop through columns
		for j in range(UQ.shape[1]):
			# Solve dx stuff one at a time
			duq_dx.append((uqrow[(j+1)%UQ.shape[1]] - uqrow[(j-1)])/dx[i])
			duqhdo_dx.append((uqhdorow[(j+1)%UQ.shape[1]] - uqhdorow[(j-1)])/dx[i])
			
			################
			### dy start ###
			################
			# Solve all the dy stuff here
			# Only have to do it once inside the i loop
			if i == 0:
				# Initialize dy carriers as 0 for the top lat
				duq_dy = [0]
				duqhdo_dy = [0]

				# Extract the column data for this iteration
				uqcol = UQ[:,j]
				uqhdocol = UQ_HDO[:,j]

				# Loop through the extraacted column and compute the d/dy
				for n in range(len(uqcol)):
					# This is for top and bottom lat, we hard code in 0's already
					if (n == 0) | (n == len(uqcol)-1):
						continue
					else:
						duq_dy.append((uqcol[n+1] - uqcol[n-1])/dy)
						duqhdo_dy.append((uqhdocol[n+1] - uqhdocol[n-1])/dy)

				# Add another zero for the -90 lat
				duq_dy.append(0)
				duqhdo_dy.append(0)

				# Keep track with master array
				dUQ_dy[:,j] = duq_dy
				dUQHDO_dy[:,j] = duqhdo_dy
			###############
			### end dy ####
			###############

		# Keep track with master array
		dUQ_dx[i,:] = duq_dx
		dUQHDO_dx[i,:] = duqhdo_dx

	dR_dy = -dUQHDO_dy + np.multiply(R, dUQ_dy)
	dR_dx = -dUQHDO_dx + np.multiply(R, dUQ_dx)
	dR_dt = np.divide(dR_dx + dR_dy, Q)
	return dR_dt


# =========================================================================================== #


############################
### ExtractVarAtPressure ###
############################

def extract_isobar(Dataset, var, pressure=50000, box = None):
	'''
	Author : Kyle Niezgoda
	Date : June 15, 2017

	This function extracts 2d data from a netCDF4 Dataset object at a specified pressure level.
	In CESM, model levels are surface following, so the level does not reflect the true pressure.
	Instead, P = A + B*PS. To get a specific pressure, you have to interpolate between levels. 
	Assuming the variable interpolates the same way (linearly), you can extract the variable value
	at the specified pressure level.

	Arguments
	---------

	Returns
	-------
	array : 2-d array of var at specified pressure
	'''
	import numpy as np
	from pygoda import find_indices

	# Extract necessary data
	lats = Dataset.variables['lat'][:].squeeze()
	lons = Dataset.variables['lon'][:].squeeze()
	PS = Dataset.variables['PS'][:].squeeze()
	A = Dataset.variables['hyam'][:].squeeze()
	B = Dataset.variables['hybm'][:].squeeze()
	var = Dataset.variables[var][:].squeeze()

	if box is not None:
		idx = find_indices(box, lats, lons)
		lats = lats[idx[0] : idx[1]]
		lons = lons[idx[2] : idx[3]]
		PS = PS[idx[0] : idx[1], idx[2] : idx[3]]
		var = var[:, idx[0] : idx[1], idx[2] : idx[3]]

	# Create 3d pressure field
	P = []
	for n in range(len(B)):
		P.append(np.array(A[n] * 100000 + B[n] * PS))
	P = np.array(P)

	# Initialize a mster array
	VAR = np.zeros(len(lons))

	for i in range(len(lats)):
		# Initialize empty array to hold var at this latitude
		var_array_by_latzone_hold = []

		for j in range(len(lons)):
			# Extract pressure and var arrays for this lat/lon coordinate
			pressure_hold = P[:,i,j]
			var_hold = var[:,i,j]

			# Here, "below" and "above" mean in space, thus 60000 is "below" 50000 even though 6 > 5
			pressures_above = pressure_hold - pressure < 0
			pressures_below = pressure_hold - pressure > 0
			
			# Check for the case when all pressures are above our value for p
			# This can happen at the poles or over mountain ranges where PS values are < our p value
			# Might want to set this to just return the value at the ground instead of NaN
			if sum(pressures_below) == 0:
				var_array_by_latzone_hold.append(np.nan)

			else:
				specific_pressure_above = pressure_hold[pressures_above][-1]
				specific_pressure_below = pressure_hold[pressures_below][0]

				# Computes the fraction of the pressurelevel that 500 falls in
				# We will integrate from top down until we hit the last pressure level that is less than pressure
				# At that point, we use the fraction of the next pressure jump that will correspond to pressure
				# E.g: if specific_pressure_above = 450mb and specific_pressure_below = 550mb, then f = .5
				# If specific_pressure_above = 499mb and specific_pressure_below = 599mb, then f = .01 
				f = (pressure - specific_pressure_above) / abs(specific_pressure_below - specific_pressure_above)

				specific_var_above = var_hold[pressures_above][-1]
				specific_var_below = var_hold[pressures_below][0]

				# Now we know f, want to solve for var at this i,j coordinate
				# This assumes the variable changes linearly in the z direction
				# Might want to code in some other ways to interpolate var values
				# e.g., wind speeds increase exponentially with height, not linearly
				# Well-mixed substances (temp, water vapor) should change linearly 
				var_ij = f * abs(specific_var_below - specific_var_above) + specific_var_above
				
				# Add to the latitude band hold array
				var_array_by_latzone_hold.append(var_ij)

		# Add to the master array
		VAR = np.vstack([VAR, var_array_by_latzone_hold])

	# Get rid of the initialzation row (full of zeros)
	VAR = VAR[1:,:]
	return VAR


# =========================================================================================== #


############################
### ExtractVarAtPressure ###
############################

def isobar(data, pressure, PS, hyam, hybm, lats, lons):
	'''
	Author : Kyle Niezgoda
	Date : June 15, 2017

	This function extracts 2d data from a netCDF4 Dataset object at a specified pressure level.
	In CESM, model levels are surface following, so the level does not reflect the true pressure.
	Instead, P = A + B*PS. To get a specific pressure, you have to interpolate between levels. 
	Assuming the variable interpolates the same way (linearly), you can extract the variable value
	at the specified pressure level.

	Arguments
	---------
	data : a 3d numpy array of data
	pressure : int of pressure in Pa
	PS : PS from cam history file
	hyam : hyam from cam history file
	hybm : see above
	lats : lats to be used, must be same len as PS.shape[1]
	lons : lons to be used, must be same len as PS.shape[2]

	Returns
	-------
	array : 2-d array of var at specified pressure
	'''
	import numpy as np
	from pygoda import PressureCalc

	# Create 3d pressure field
	P = PressureCalc(hyam,hybm,PS)

	# Check for problems
	if P.shape[1] != len(lats):
		print "Input PS array does not have the same latitudes at the lats argument! Exiting..."
		return None

	if P.shape[2] != len(lons):
		print "Input PS array does not have the same longitudes at the lons argument! Exiting..."
		return None

	# Initialize a master array
	VAR = np.zeros(len(lons))

	for i in range(len(lats)):
		# Initialize empty array to hold var at this latitude
		var_array_by_latzone_hold = []

		for j in range(len(lons)):
			# Extract pressure and var arrays for this lat/lon coordinate
			pressure_hold = P[:,i,j]
			var_hold = data[:,i,j]

			# Here, "below" and "above" mean in space, thus 60000 is "below" 50000 even though 6 > 5
			pressures_above = pressure_hold - pressure < 0
			pressures_below = pressure_hold - pressure > 0
			
			# Check for the case when all pressures are above our value for p
			# This can happen at the poles or over mountain ranges where PS values are < our p value
			# Might want to set this to just return the value at the ground instead of NaN
			if sum(pressures_below) == 0:
				var_array_by_latzone_hold.append(np.nan)

			else:
				specific_pressure_above = pressure_hold[pressures_above][-1]
				specific_pressure_below = pressure_hold[pressures_below][0]

				# Computes the fraction of the pressurelevel that 500 falls in
				# We will integrate from top down until we hit the last pressure level that is less than pressure
				# At that point, we use the fraction of the next pressure jump that will correspond to pressure
				# E.g: if specific_pressure_above = 450mb and specific_pressure_below = 550mb, then f = .5
				# If specific_pressure_above = 499mb and specific_pressure_below = 599mb, then f = .01 
				f = (pressure - specific_pressure_above) / abs(specific_pressure_below - specific_pressure_above)

				specific_var_above = var_hold[pressures_above][-1]
				specific_var_below = var_hold[pressures_below][0]

				# Now we know f, want to solve for var at this i,j coordinate
				# This assumes the variable changes linearly in the z direction
				# Might want to code in some other ways to interpolate var values
				# e.g., wind speeds increase exponentially with height, not linearly
				# Well-mixed substances (temp, water vapor) should change linearly 
				var_ij = f * abs(specific_var_below - specific_var_above) + specific_var_above
				
				# Add to the latitude band hold array
				var_array_by_latzone_hold.append(var_ij)

		# Add to the master array
		VAR = np.vstack([VAR, var_array_by_latzone_hold])

	# Get rid of the initialzation row (full of zeros)
	VAR = VAR[1:,:]
	return VAR


###############
### h1dates ###
###############

def h1dates(sy, ey, sm = 1, em = 12):
	'''
	Author : Kyle Niezgoda
	Date : July 27, 2017

	This function creates date string for use in glob to search for the right
	h1 cam files. A list of dates of the form yyyy-mm-dd that starts at 
	sy-sm-01 and ends at ey-em-31 is created, which can be looped over to 
	search for h1 files you needs to use.

	Arguments
	---------
	sy : start year, an integer
	ey : end year, int
	sm : start month, int (jan = 1, dec = 12)
		default : 1
	em : end month, int
		default : 12

	Returns
	-------
	list : list of 'yyyy-mm-dd' strings, same form as cam h1 file names.
	'''
	start_year = sy
	end_year = ey
	start_month = sm
	end_month = em

	# Contruct dates to look for
	# All months will have 31 days, but that is ok because glob will ignore any files it doesn't find
	years = []
	if end_year < 10:
		years = years + ["000" + str(y) for y in range(start_year,(end_year+1))] 
	elif end_year >= 10:
		if start_year < 10:
			years = years + ["000" + str(y) for y in range(start_year,10)]
			years = years + ["00" + str(y) for y in range(10,(end_year+1))]
		if start_year > 10:
			years = years + ["00" + str(y) for y in range(start_year,(end_year+1))]


	months = []
	if end_month < 10 :
		months = months + ["0" + str(y) for y in range(start_year,(end_month+1))] 
	elif end_month >= 10:
		if start_month < 10:
			months = months + ["0" + str(y) for y in range(start_month,10)]
			months = months + [str(y) for y in range(10,(end_month+1))]
		if start_month > 10:
			months = months + [str(y) for y in range(start_month,(end_month+1))]


	dates = []
	for y in years:
		for m in months:
			dates = dates + [y + "-" + m + "-0" + str(d) for d in range(1,10)]
			dates = dates + [y + "-" + m + "-" + str(d) for d in range(10,32)]

	return dates


# =========================================================================================== #


############
### h1ts ###
############

def h1ts(var, box, math, sy, ey, sm = 1, em = 12):
	'''
	Author : Kyle Niezgoda
	Date : July 28, 2017

	This function calculates area averages or sums of "var" and reports them
	as timeseries data. The dates are all the dates between sy-sm-01 and ey-em-31 
	(sy = "start year", etc.) Data are extracted from cam.h1 files,
	so the function must be executed in a directory with files that will 
	match the wildcard expression *cam.h1.yyyy-mm-dd*.nc where yyyy, mm, and
	dd are strings created inside this function from the input arguments, 
	sy, ey, sm, and em.

	The data are either averaged or summed according to the math argument.

	Arguments
	---------
	var : string, the name of a 2-dimensional variable in the h1 files, all sensitive (spelling, case)
		***This function supports calculation of delta - var can be prect_d18O, prect_dD, d18OV, and dDV
		for vapor, only the lowest model level can be reported (have to create new h1 files with 3d data
		if you want to report other levels for vapor data)***
	box : list, region over which to extract data, [bottom, top, left, right]
	math : string, currently only supports "sum" or "mean"/"average"
	sy : int, start year
	ey : int, end year
	sm : int, start month (jan = 1, dec = 12)
		default : 1
	em : int, end month
		default : 12

	Returns
	-------
	[0] : var area-averaged or -summed data over the timeseries
	[1] : list of dates used
	[2] : list of files used
	'''
	import os, glob
	from pygoda import h1dates, find_indices
	from netCDF4 import Dataset
	import pandas as pd
	import numpy as np

	math_lower = math.lower()
	var_lower = var.lower()

	# special_vars = ['prect_d18o', 'prect_dd', 'd18ov', 'ddv']

	dates = h1dates(sy, ey, sm, em)

	southern_lat = box[0]
	northern_lat = box[1] 
	left_lon = box[2] 
	right_lon = box[3]
	
	VAR = []
	files = []
	used_dates = []
	first = True
	for n, d in enumerate(dates):
		# Find the file for this date
		f = glob.glob("*cam.h1."+d+"*.nc")
		
		if len(f) > 0: # If a file was found
			# Open the netcdf file for reading
			nc = Dataset(f[0], 'r')
			
			while first:
				print "Extracting lat and lon data, will only do this once"
				# Only have to do this once
				lat = nc.variables['lat'][:]
				lon = nc.variables['lon'][:]
				bounds = find_indices([southern_lat, northern_lat, left_lon, right_lon], lat, lon)
				first = False

			# Extract the data, checking for special cases
			if var_lower == 'prect_d18o':
				hold = (nc.variables['PRECT_H218O'][:,bounds[0]:bounds[1],bounds[2]:bounds[3]]/nc.variables['PRECT_H2O'][:,bounds[0]:bounds[1],bounds[2]:bounds[3]] - 1) * 1000
			elif var_lower == 'prect_dd':
				hold = (nc.variables['PRECT_HDO'][:,bounds[0]:bounds[1],bounds[2]:bounds[3]]/nc.variables['PRECT_H2O'][:,bounds[0]:bounds[1],bounds[2]:bounds[3]] - 1) * 1000
			elif var_lower == 'd18ov':
				hold = (nc.variables['H218OVBT'][:,bounds[0]:bounds[1],bounds[2]:bounds[3]]/nc.variables['H2OVBT'][:,bounds[0]:bounds[1],bounds[2]:bounds[3]] - 1) * 1000
			elif var_lower == 'ddv':
				hold = (nc.variables['HDOVBT'][:,bounds[0]:bounds[1],bounds[2]:bounds[3]]/nc.variables['H2OVBT'][:,bounds[0]:bounds[1],bounds[2]:bounds[3]] - 1) * 1000
			else:	
				hold = nc.variables[var][:,bounds[0]:bounds[1],bounds[2]:bounds[3]]

			# Do arithmetic on the data
			if math_lower == "sum":
				VAR.append(np.nansum(hold))
			elif (math_lower == "mean") | (math_lower == "average"):
				VAR.append(np.nanmean(hold))

			#Progress tracker
			if n % (len(dates)/10) == 0:
				print str(round(float(n)/float(len(dates)) * 100)) + "% done extracting ts data..."
		
			# Track the data
			files.append(f[0])
			used_dates.append(d)
	return [np.array(VAR), used_dates, files]


# =========================================================================================== #


#####################
### isnot_outlier ###
#####################

def isnot_outlier(points, thresh=3.5):
    """
    Copied and slightly modified from: 
    https://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
    
    Returns a boolean array with False if points are outliers and True 
    otherwise.

    Parameters:
    -----------
    points : An numobservations by numdimensions array of observations
    thresh : The modified z-score to use as a threshold. Observations with
        a modified z-score (based on the median absolute deviation) greater
        than this value will be classified as outliers.

    Returns:
    --------
    mask : A numobservations-length boolean array.

    References:
    ----------
    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
    Handle Outliers", The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    import numpy as np
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.nanmedian(points, axis=0)
    diff = np.nansum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score > thresh


# =========================================================================================== #


#####################
### findClimoFile ###
#####################

def findClimoFile(season, directory = '.'):
	"""
	Finds a climo cam h0 file from a directory using glob
	"""
	import os, glob, sys
	fname = glob.glob(directory+"/*"+season+"*.nc")
	try:
		path = os.path.abspath(fname[0])
		name = os.path.splitext(os.path.split(fname[0])[1])[0]
		return [path, name]
	except IndexError:
		print "No file found!"
		return [0,0]


# =========================================================================================== #


######################
### zeroCenterClev ###
######################

def zeroCenterClev(d, alpha = .99, verb = False):
	import numpy as np
	d = d[~np.isnan(d)]
	numcells = d.size
	MAX = np.max(d)
	MIN = np.min(d)
	noZero = False
	if (MAX < 0) | (MIN > 0):
		noZero = True
	
	if noZero:
		if verb:
			print "No way to center color bar on zero"
		return 19
	else:
		delta = MAX - MIN
		test = min(abs(MIN), abs(MAX))
		a = len(np.where(np.logical_and(d>-test, d<test))[0])
		frac = float(a)/float(numcells)
		trace = 0
		while frac < alpha:
			test = test + delta/20
			a = len(np.where(np.logical_and(d>-test, d<test))[0])
			frac = float(a)/float(numcells)
			trace += 1
			if trace > 1000:
				if verb:
					print "Infinite loop in the diffColor function.\nExiting function with return value of 19"
				return 19
	
	return np.linspace(-test, test, 21)


# =========================================================================================== #


def niceClev(data, minnstep = 8, maxnstep = 22, alpha = .99, verb = False):
	from pygoda import zeroCenterClev, ncgoda
	import numpy as np
	import math
	#
	not_pretty_clevs = zeroCenterClev(data, alpha, verb)
	if type(not_pretty_clevs) == int:
		return not_pretty_clevs
	else:
		maxmin = abs(not_pretty_clevs[0])
	delta = maxmin*2
	step = (1, 2, 3, 4, 5)
	power_ten = (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5)
	#
	found = False
	for p in power_ten:
		if found:
			break
		#
		for s in step:
			if found:
				break
			#
			nstep = math.ceil(delta/(s * (10**p)))
			if verb:
				print nstep
			#
			if (minnstep <= nstep) & (nstep <= maxnstep) & (nstep % 2 == 0):
				if verb:
					print "Found a combo that works"
				#
				found = True
				STEP = s * (10**p)
				NSTEP = nstep
			#
	#
	if not found:
		if verb: 
			print "cant do it"
		#
		return not_pretty_clevs
	#
	else:
		# if NSTEP % 2 == 0:
		return np.linspace(-float(NSTEP)/2*STEP, float(NSTEP)/2*STEP, num = NSTEP+1)
		# if NSTEP % 2 == 1:
		# 	return np.linspace(-float(NSTEP+1)/2, num = NSTEP+1)



##############
### ncgoda ###
##############
	
class ncgoda:
	"""
*** Class ncgoda ***

Class useful for working with iCAM5 history files. Instances are initiated with .dataset, .lat, and .lon
from the netcdf dataset. 

====================
!!!Need to update!!!
====================
1) Right now, class is not safely capable of handling time dimensions with length > 1.
That means ncgoda should be used to open .nc files with only 1 time slice.

Not the best way to do it, and should incorporate time dimensions into isobar and variable methods for completion.

Right now, best way to deal with this is to just loop over files with 1 time slice and create a new
ncgoda instance each time. 

---

2) Incorporate way to slice model level rather than isobar.

Along with this, might want to make a way to show all the 3d data rather than a slice.

---
=============
How it works:
=============

Initialize:
-----------
Arugments : [0] dataset (netCDF4 Dataset)

Methods:
--------
variable : returns 2 or 3d numpy array data. 
	Sets : .var, .box, .data, .boxlat, .boxlon, .units, and .long_name
	Arugments : [0] string variable name; [1] region box (oprtional) in form of [bottom, top, left, right]

isobar : returns 2d np array along an isobar
	Sets : .pressure, .data
	Arguments : [0] string variable name; [1] integer pressure in Pa; [2] region box (optional)

PRECT_d18O and PRECT_dD : returns 2d numpy array data. 
	Sets : .h2o, .h218o/.hdo, .var, .data, and .long_name
	Arguments : [0] region box (optional)

d18OV and dDV : returns 2d numpy array data. 
	Sets : .h2ov_p, .h218ov_p/.hdov_p, .var, .data, and .long_name
	Arguments : [0] integer pressure in Pa; [1] region box (optional)

	"""
	def __init__(self, dataset):
		from netCDF4 import Dataset as ds
		from pygoda import ncdump
		self.dataset = dataset
		if type(dataset) is str:
			self.dataset = ds(dataset, "r")
		self.dims = [d for d in self.dataset.dimensions]
		self.vars = [v for v in self.dataset.variables]
		self.dimlen = [len(self.dataset.dimensions[d]) for d in self.dims]
		self.lat = self.dataset.variables['lat'][:]
		self.lon = self.dataset.variables['lon'][:]
		self.isTime = any(d == "time" for d in self.dims)
		if self.isTime:
			self.ntime = self.dimlen[self.dims.index("time")]
			if self.ntime == 1:
				self.isTime = False

	def variable(self, var, box = None, verb = False):
		'''
		Main function for extracting data from netcdf file
		Performs conversions on data if there is documentation for it in algebra.py
		'''
		from pygoda import find_indices
		import xarray as xr
		import numpy as np

		# Read the variable 
		data = self.dataset.variables[var]

		# Set some variables for the whole class to see
		self.vartype = "2d"
		self.var = var
		self.box = box

		# Read in numbers and units for the data
		self.units = data.units
		self.data = data[:].squeeze()

		ndims = len(self.data.shape)
		
		# Strip the numbers down to the region from the box arg
		idxbox = find_indices(box, self.lat, self.lon)

		# No lev, no time
		if ndims == 2:
			self.data = np.array(xr.DataArray(self.data)[idxbox[0], idxbox[1]])
		# Either lev or time but not both
		if ndims == 3:
			self.data = np.array(xr.DataArray(self.data)[:, idxbox[0], idxbox[1]])
		if ndims == 4:
			self.data = np.array(xr.DataArray(self.data)[:, :, idxbox[0], idxbox[1]])

		# Do conversion (default of no conversion if there isn't documentation for a conversion in algebra.py)
		try:
			from algebra import algebra
			# Convert the data using algebra.py file
			# Set new units to correspond to the algebra changes
			mult, add, units = algebra.get(self.var, [1, 0, self.units]) # default to no change if var is not found
			self.data = (self.data * mult) + add
			self.units = units
			if verb:
				if (mult != 1) | (add != 0):
					print "From pygoda.variable(): Converted data as follows:\n\tmultiplied by " + str(mult) + "\n\tadded " + str(add)
		except ImportError:
			print "Could not find algebra.py, not able to do conversion!"

		self.boxlat = self.lat[idxbox[0]]
		self.boxlon = self.lon[idxbox[1]]
		self.long_name = self.dataset.variables[var].long_name
		return self.data

	def isobar(self, var, pressure, box = None, verb = False):
		from pygoda import PressureCalc
		import numpy as np
		self.pressure = pressure
		PS = self.variable("PS", box)
		data = self.variable(var, box, verb)
		if pressure <= 30:
			'''
			Add in feature that can ready pressure args <= 30 
			and print model level rather than actuall isobar
			'''
			
			# print "Extracting along model level rather than isobar"
			# print "Will return model level " + str(pressure)
			None
		A = self.dataset.variables["hyam"][:]
		B = self.dataset.variables["hybm"][:]

		# Create 3d pressure field
		P = PressureCalc(A,B,PS)
		ndims = len(P.shape)

		# Initialize a master array
		VAR = np.zeros(len(self.boxlon))

		for i in range(len(self.boxlat)):
			# Initialize empty array to hold var at this latitude
			var_array_by_latzone_hold = []

			for j in range(len(self.boxlon)):
				# Extract pressure and var arrays for this lat/lon coordinate
				pressure_hold = P[:,i,j]
				var_hold = data[:,i,j]

				# Here, "below" and "above" mean in space, thus 60000 is "below" 50000 even though 6 > 5
				pressures_above = pressure_hold - pressure < 0
				pressures_below = pressure_hold - pressure > 0
				
				# Check for the case when all pressures are above our value for p
				# This can happen at the poles or over mountain ranges where PS values are < our p value
				# Might want to set this to just return the value at the ground instead of NaN
				if sum(pressures_below) == 0:
					var_array_by_latzone_hold.append(np.nan)

				else:
					specific_pressure_above = pressure_hold[pressures_above][-1]
					specific_pressure_below = pressure_hold[pressures_below][0]

					# Computes the fraction of the pressurelevel that 500 falls in
					# We will integrate from top down until we hit the last pressure level that is less than pressure
					# At that point, we use the fraction of the next pressure jump that will correspond to pressure
					# E.g: if specific_pressure_above = 450mb and specific_pressure_below = 550mb, then f = .5
					# If specific_pressure_above = 499mb and specific_pressure_below = 599mb, then f = .01 
					# f = (pressure - specific_pressure_above) / abs(specific_pressure_below - specific_pressure_above)
					# New way
					f = (pressure - specific_pressure_above) / (specific_pressure_below - specific_pressure_above)

					specific_var_above = var_hold[pressures_above][-1]
					specific_var_below = var_hold[pressures_below][0]

					# Now we know f, want to solve for var at this i,j coordinate
					# This assumes the variable changes linearly in the z direction
					# Might want to code in some other ways to interpolate var values
					# e.g., wind speeds increase exponentially with height, not linearly
					# Well-mixed substances (temp, water vapor) should change linearly 
					# var_ij = f * abs(specific_var_below - specific_var_above) + specific_var_above
					# New way
					var_ij = f * (specific_var_below - specific_var_above) + specific_var_above
					
					# Add to the latitude band hold array
					var_array_by_latzone_hold.append(var_ij)

			# Add to the master array
			VAR = np.vstack([VAR, var_array_by_latzone_hold])

		# Get rid of the initialzation row (full of zeros)
		VAR = VAR[1:,:]
		self.data = VAR
		self.vartype = "3d"
		return self.data

	def prep_map(self, season, region):
		try:
			from clevs import getClev, getCmap
			import matplotlib.pyplot as plt
			pressure = ""
			if self.vartype == "3d":
				pressure = self.pressure
			self.clevs = getClev(self.var+str(pressure), season+"_"+region, self.data)
			self.cmap = getCmap(self.var, 'cmap')
			self.diffcmap = getCmap(self.var, 'diffcmap')
		except ImportError:
			return "Could not find clevs.py, not able to set clevs!"


	def PRECT_d18O(self, box = None):
		h2o = self.variable("PRECT_H2O", box)
		h218o = self.variable("PRECT_H218O", box)
		self.var = "PRECT_d18O"
		self.long_name = "d18O of PRECT"
		self.units = "permil"
		self.data = (h218o / h2o - 1) * 1000
		return self.data

	def PRECT_dD(self, box = None):
		h2o = self.variable("PRECT_H2O", box)
		hdo = self.variable("PRECT_HDO", box)
		self.var = "PRECT_dD"
		self.long_name = "dD of PRECT"
		self.units = "permil"
		self.data = (hdo / h2o - 1) * 1000
		return self.data

	def PRECT_dxs(self, box = None):
		dd = self.PRECT_dD(box)
		d18o = self.PRECT_d18O(box)
		self.var = "PRECT_dxs"
		self.long_name = "d-excess of PRECT"
		self.units = "permil"
		self.data = dd - 8 * d18o
		return self.data

	def d18OV(self, pressure, box = None):
		h218ov_p = self.isobar("H218OV", pressure, box) 
		h2ov_p = self.isobar("H2OV", pressure, box) 
		self.var = "d18OV"
		self.long_name = "d18O of Vapor at "
		self.units = "permil"
		self.data = (h218ov_p / h2ov_p - 1) * 1000
		return self.data

	def dDV(self, pressure, box = None):
		hdov_p = self.isobar("HDOV", pressure, box) # only talk once
		h2ov_p = self.isobar("H2OV", pressure, box) 
		self.var = "dDV"
		self.long_name = "dD of Vapor at "
		self.units = "permil"
		self.data = (hdov_p / h2ov_p - 1) * 1000
		return self.data

	def dxsV(self, pressure, box = None):
		ddv_p = self.dDV(pressure, box)
		d18ov_p = self.d18OV(pressure, box)
		self.var = "dxsV"
		self.long_name = "d-excess of VAPOR at "
		self.units = "permil"
		self.data = ddv_p - 8 * d18ov_p
		return self.data