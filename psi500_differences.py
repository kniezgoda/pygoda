'''
Code to calculate Hadley Cell location data using the Kang and Polvani metrics
Kyle Niezgoda, Spring 2017
niezgodk@oregonstate.edu

We calculate the zonal-mean meridional mass stream function (psi) at every vertical level by

(1) psi = V * P * cos(lat) / g

where V is meridional wind speed in m/s, P is pressure in Pa, and g is the acceleration due to gravity, 9.8 m/s2

We then calculate psi500 by summing psi from the top of the atmosphere to the 500mb (50000 Pa) pressure level.
Model pressure coords do not land exactly at 500mb, so we sum all psi values above the last pressure value before 500 
and then add a fraction of the next pressure level which corresponds to the position of 500mb relative to the pressure gap it falls in:

Example:
---------------483mb pressure level--------------

(500mb falls in the middle of this gap)

---------------529mb pressure level--------------

To calculate f using the example above, we do abs(500-438)/abs(529-438). 
The final calulation of psi then becomes:

(2) sum(psi values above 500mb) + f * sum(psi value of the pressure gap that 500mb falls in)

We then calculate the latitude at which psi500 changes sign going poleward from the extremum, as in Kang and Polvani.
This latitude represents the northern and southern edge of the Hadley Cell.
The zero-crossing between these two points is the center of the ITCZ, although Kang and Polvani do not report it as such.
'''

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
np.set_printoptions(threshold=np.inf)
def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Set the input file
# Multiple files will result in maps displayed one after another of each file
directory = "/home/server/student/homes/kniezgod/model_runs"
Control = "F.C5.2deg.wiso.defaultSSTICE_kn002/F.C5.2deg.wiso.defaultSSTICE_kn002_ANN_climo.nc"
Test = "F.C5.2deg.wiso.obs6kSST_kn003/F.C5.2deg.wiso.obs6kSST_kn003_kn003_ANN_climo.nc"

show_data = True
write_csv = False
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

g = 9.8 #m/s2

# Initialize empty df and start for loop through files
# InfoDataFrame = pd.DataFrame(columns=('File', 'Time', 'ITCZcenter', 'HadleyCellNorthernEdge', 'HadleyCellSouthernEdge', "NorthernCellWidth", 'SouthernCellWidth'))



# file name for plotting
controlf_name = os.path.splitext(os.path.split(Control)[1])[0]

#Read the netcdf file
Control_data = Dataset(directory+Control, mode='r')

#Extract data
time = data.variables['time'][:]
# We want to cycle through times if there is more than one time in the netcdf file
for timenum, t in enumerate(time):
	
	#----------Begin for loop through time----------#

	if len(time) == 1:
		month = "Climo" # This is for dataframe creation at the end
	else:
		month = "Month " + str(timenum+1) # Also only for data frame creation at the end
	
	lats = data.variables['lat'][:]
	lons = data.variables['lon'][:]
	levs = data.variables['lev'][:] * 100
	V = data.variables['V'][timenum,:,:,:]
	A = data.variables['hyam'][:]
	B = data.variables['hybm'][:]
	PS = data.variables['PS'][timenum,:,:]
	
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	# We now want to compute psi at 500 mbar
	# We let dpsi/dz = rho*v 
	# We know rho is p/gh from hydrostatic equation
	# So, d-psi/dz = pv/gh and, in discrete terms, delta-psi = pv/g
	# Thus we are left to calculate p and v

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	# Pressure calcs
	# Pressure is not as simple as looking at levs because CAM has terrain following pressure scheme
	# Pressure = A*P0 + B*PS
	# P0 = 100000, PS = surface pressure from cam

	#creates 3-d pressure array 
	pressure = []
	for n in range(len(B)):
		pressure.append(np.array(A[n] * 100000 + B[n] * PS))

	#print some stuff
	pressure = np.array(pressure)
	# print np.shape(pressure) # 30, 96, 144

	#average over longitudes to compute zonal average
	zonal_pressure = np.mean(pressure,2) # dimensions (lev, lat)

	# print lats[15]
	# print zonal_pressure[:,15] #vertical pressure field at latitude lats[15]

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	# Wind speed calcs 
	zonal_V = np.mean(V, 2) # dimensions (lev, lat)

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	# Psi calcs
	# First calculate the lat by lev psi grid
	# Have to multiply psi by cos(latitude) to get the area-normalized mass 
	# Then, calculate psi500 = integrate from top down until reaching 500mb pressure level

	psi = []
	for n, l in enumerate(lats):
		psi.append(zonal_V[:,n] * zonal_pressure[:,n] * np.cos(l*np.pi/180) / g)

	psi = np.array(psi)
	# print np.shape(psi) # (96, 30)

	psi500 = []
	for n, l in enumerate(lats):
		pressure_hold = zonal_pressure[:,n]
		#
		pressures_above500mb = pressure_hold - 50000 < 0
		pressures_below500mb = pressure_hold - 50000 > 0
		#
		last_pressure_before500mb = pressure_hold[pressures_above500mb][-1]
		first_pressure_after500mb = pressure_hold[pressures_below500mb][0]
		fractional_level_index = np.argmax(pressure_hold>50000)
		#
		# Computes the fraction of the pressurelevel that 500 falls in
		# We will integrate from top down until we hit the last pressure level that is less than 500mb
		# At that point, we use the fraction of the next pressure jump that will correspond to 500mb
		# E.g: if last = 450mb and first = 550mb, then f = .5
		# If last = 499mb and first = 599mb, then f = .01 
		f = (50000 - last_pressure_before500mb) / abs(first_pressure_after500mb - last_pressure_before500mb)
		#
		#compute the psi fractional addition ---> that pressure jump which 500mb falls in
		fractional_addition = (f * psi[n, fractional_level_index])
		#
		#compute psi500 ---> sum(psi above 500mb) + fractional component of pressure jump that 500mb falls in
		psi500.append(np.sum(psi[n,pressures_above500mb]) + fractional_addition)


	psi_df = pd.DataFrame({'latitude' : lats, 'psi500' : psi500})

	# plt.plot(psi_df['latitude'], psi_df['psi500'])
	# plt.show()

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	# Now, if we assume that changes in psi are linear, which they pretty much are (look at graph to be sure),
	# we can calculate the "red-x" value from Kang and Polvani
	# which is phi-psi-500, or the latitude where psi500 switches signs going poleward from the maxima. 

	#find the absolute value maximum 
	absmax = max(abs(psi_df['psi500']))
	absmax_index = [n for n, m in enumerate(psi_df['psi500']) if abs(m) == absmax][0]

	#find the value at the absmax index
	absmax_psival = psi_df['psi500'][absmax_index]

	################################################################################################
	# We have a situation where the absmax can be either positive or negative. If the value is 
	# positive, we can go north to the first zero crossing and locatate the NH edge. We can then go 
	# south to the first zero crossing to find the ITCZ, and again to the second zero-crossing to
	# find the SH edge. If absmax is negative, we go south once (SH edge) and north twice (ITCZ and
	# NH edge). 
	################################################################################################
	if absmax_psival > 0:
		#find the first negative value after the psimax and identify it's index
		#found argmax in a google search, not sure what it really does or why it works
		first_neg_index = np.argmax(psi_df['psi500'][absmax_index:] < 0)
		
		#calculate the value at and above that index
		first_neg_value = psi_df['psi500'][first_neg_index]
		last_pos_value = psi_df['psi500'][first_neg_index-1]
		
		#find the latitudes at those two points
		first_neg_lat = psi_df['latitude'][first_neg_index]
		last_pos_lat = psi_df['latitude'][first_neg_index-1]
		
		#zero should fall in between first_neg_value and last_pos_value
		#calculate the fraction of the psi gap between first_neg and last_pos that 0 falls in 
		f = abs(last_pos_value) / abs(last_pos_value - first_neg_value)
		
		#calculate the latitude using f and assuming lat changes are linear in this lat gap
		zerolat = last_pos_lat + f * abs(first_neg_lat - last_pos_lat)
		
		#this latitude is the northern edge of the hadley cell
		HadleyCellNorthernEdge = zerolat
		
		################################################################################################
		# now go backwards from the absmax_index and find the next zero-crossing to find the ITCZ center
		# following the same basic algorithm as above
		################################################################################################
		first_neg_index = np.argmax(psi_df['psi500'][list(reversed(range(0,absmax_index)))] < 0)
		
		first_neg_value = psi_df['psi500'][first_neg_index]
		last_pos_value = psi_df['psi500'][first_neg_index+1]
		first_neg_lat = psi_df['latitude'][first_neg_index]
		last_pos_lat = psi_df['latitude'][first_neg_index+1]
		
		f = abs(last_pos_value) / abs(last_pos_value - first_neg_value)
		
		#going up the data frame (decreasing latitude), so we subtract the fractional latitude change from the last positive latitude
		zerolat = last_pos_lat - f * abs(first_neg_lat - last_pos_lat)
		#this latitude is the center of the ITCZ
		ITCZcenter = zerolat
		
		#################################################################################################
		# now go backwards from the ITCZ center and find the next zero-crossing to find the SH Hadley 
		# edge. Same basic algorithm. 
		#################################################################################################
		first_pos_index = np.argmax(psi_df['psi500'][list(reversed(range(0,first_neg_index)))] > 0)
		
		first_pos_value = psi_df['psi500'][first_pos_index]
		last_neg_value = psi_df['psi500'][first_pos_index+1]
		first_pos_lat = psi_df['latitude'][first_pos_index]
		last_neg_lat = psi_df['latitude'][first_pos_index+1]
		
		f = abs(last_neg_value) / abs(last_neg_value - first_pos_value)
		
		#going up the data frame so subtract the fractional latitude change from the last positive latitude
		zerolat = last_neg_lat - f * abs(first_neg_lat - last_pos_lat)
		HadleyCellSouthernEdge = zerolat 


	##############################################
	# Do practically the same thing for absmax < 0
	# Just kind of flip the operation around
	##############################################
	if absmax_psival < 0:
		first_pos_index = np.argmax(psi_df['psi500'][list(reversed(range(0,absmax_index)))] > 0)
		
		first_pos_value = psi_df['psi500'][first_pos_index]
		last_neg_value = psi_df['psi500'][first_pos_index+1]
		first_pos_lat = psi_df['latitude'][first_pos_index]
		last_neg_lat = psi_df['latitude'][first_pos_index+1]
		
		f = abs(last_neg_value) / abs(last_neg_value - first_pos_value)
		
		#going up the data frame so subtract the fractional latitude change from the last positive latitude
		zerolat = last_neg_lat - f * abs(last_neg_lat - first_pos_lat)
		HadleyCellSouthernEdge = zerolat 
		
		################################################################################################
		# now go forwards from the absmax_index and find the next zero-crossing to find the ITCZ center
		################################################################################################
		first_pos_index = np.argmax(psi_df['psi500'][absmax_index:] > 0)
		
		first_pos_value = psi_df['psi500'][first_pos_index]
		last_neg_value = psi_df['psi500'][first_pos_index-1]
		first_pos_lat = psi_df['latitude'][first_pos_index]
		last_neg_lat = psi_df['latitude'][first_pos_index-1]
		
		f = abs(last_neg_value) / abs(last_neg_value - first_pos_value)
		
		zerolat = last_neg_lat + f * abs(last_neg_lat - first_pos_lat)
		ITCZcenter = zerolat

		################################################################################################
		# now go forwards from the ITCZ center and find the next zero-crossing to find the NH Hadley edge
		################################################################################################
		first_neg_index = np.argmax(psi_df['psi500'][first_pos_index:] < 0)
		
		first_neg_value = psi_df['psi500'][first_neg_index]
		last_pos_value = psi_df['psi500'][first_neg_index-1]
		first_neg_lat = psi_df['latitude'][first_neg_index]
		last_pos_lat = psi_df['latitude'][first_neg_index-1]
		
		f = abs(first_neg_value) / abs(first_neg_value - last_pos_value)
		
		zerolat = last_pos_lat + f * abs(last_pos_lat - first_neg_lat)
		HadleyCellNorthernEdge = zerolat

	# Compute some other variables like cell width
	WidthNorthernCell = abs(HadleyCellNorthernEdge - ITCZcenter)
	WidthSouthernCell = abs(HadleyCellSouthernEdge - ITCZcenter)

	# Temporary variable for this instance of the for loop
	hold = pd.DataFrame({'Time' : [month], 'ITCZcenter' : [ITCZcenter], 'HadleyCellNorthernEdge' : [HadleyCellNorthernEdge], 'HadleyCellSouthernEdge' : [HadleyCellSouthernEdge], 'NorthernCellWidth' : [WidthNorthernCell], 'SouthernCellWidth' : [WidthSouthernCell]})
	
	# Log data from temp variable into master file
	InfoDataFrame = InfoDataFrame.append(hold)

	plt.plot(psi_df['latitude'], psi_df['psi500'])
	
	itczcenter, = plt.plot([ITCZcenter], 0, 'ro')
	plt.text(ITCZcenter+2, 0, "lat: " + str(round(ITCZcenter, 1)), color = 'red')
	southernedge, = plt.plot([HadleyCellSouthernEdge], 0, 'ko')
	plt.text(HadleyCellSouthernEdge+2, 0, "lat: " + str(round(HadleyCellSouthernEdge, 1)), color = 'black')
	northernedge, = plt.plot([HadleyCellNorthernEdge], 0, 'go')
	plt.text(HadleyCellNorthernEdge+2, 0, "lat: " + str(round(HadleyCellNorthernEdge, 1)), color = 'green')
	# plt.legend([itczcenter, southernedge, northernedge], ['ITCZ Center', 'Hadley cell southern edge', 'Hadley cell northern edge'])
	
	plt.xlabel("Latitude")
	plt.ylabel("zonal mean psi-500mb (Kg/m/s)")
	plt.title(f_name)
	plt.show()

		#----------End for loop through time----------#

	#----------End for loop through files----------#

InfoDataFrame = InfoDataFrame.reset_index(drop=True)
if show_data:
	print InfoDataFrame
if write_csv:
	InfoDataFrame.to_csv("HadleyData.csv")




