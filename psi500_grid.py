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


ERRORS TO FIX:
Need to multiply psi by the area of the zonal band, not the cos of the latitude. This is because we are still in 
two dimensions here, rather than one (as in the psi500.py code which averages over longitudes)
'''

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap as bm
np.set_printoptions(threshold=np.inf)
def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Set the input file
directory = "/Users/kyleniezgoda/Yellowstone_Runs/model_runs/F.C5.2deg.wiso.defaultSSTICE_kn002/"
# Files = ["F2deg.wiso.defaultSSTICE_kn001_MAM_climo.nc", "isoCAM5CLM4_2deg_wdiddledobs6kSSTs_DJF_climo.nc", "isoCAM5CLM4_2deg_wdiddledobs6kSSTs_JJA_climo.nc", "F2deg.wiso.defaultSSTICE_kn001_ts.nc"]
Files = ["F.C5.2deg.wiso.defaultSSTICE_kn002_ANN_climo.nc", "F.C5.2deg.wiso.defaultSSTICE_kn002_DJF_climo.nc"]

g = 9.8 #m/s2

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Initialize empty df and start for loop through files
InfoDataFrame = pd.DataFrame(columns=('File', 'Time', 'ITCZcenter', 'HadleyCellNorthernEdge', 'HadleyCellSouthernEdge', "NorthernCellWidth", 'SouthernCellWidth'))

for N, fi in enumerate(Files):

	#----------Begin for loop through files----------#

	#Read the netcdf file
	data = Dataset(directory+fi, mode='r')

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
		bmlon, bmlat = np.meshgrid(lons, lats)
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
		
		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		
		# 

		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		

		# Psi calcs
		# First calculate the 3-d psi grid
		# Have to multiply psi by cos(latitude) to get the area-normalized mass 
		# Then, calculate psi500 = integrate from top down until reaching 500mb pressure level

		psi_before_latnormalized = np.multiply(pressure, V) / g

		# Pretty sure we don't need to normalize by lat because we are not doing zonal means
		# Instead, we need to normalize by grid size or something like that
		psi = []
		# Normalize by latitude
		for n, lat in enumerate(lats):
			psi.append(psi_before_latnormalized[:,n,:] * np.cos(lat*np.pi/180))

		# Combine into 3-d array and move axes around so the data is in a normal format (levs, lats, lons)
		psi = np.swapaxes(np.array(psi), 0, 1)
		# print(np.shape(psi)) # (30, 96, 144)

		def psi500calc(psi, P):
			# psi and P must be same length
			# psi is a 1-d vector of psi values, where the length of the vector is the number of vertical levs
			# P is the 1-d vector of pressure values in Pa, where the length of the vector is the number of vertical levs

			pressures_above500mb = P - 50000 < 0
			pressures_below500mb = P - 50000 > 0
			
			last_pressure_before500mb = P[pressures_above500mb][-1]
			first_pressure_after500mb = P[pressures_below500mb][0]
			fractional_level_index = np.argmax(P>50000)
			
			# Computes the fraction of the pressurelevel that 500 falls in
			# We will integrate from top down until we hit the last pressure level that is less than 500mb
			# At that point, we use the fraction of the next pressure jump that will correspond to 500mb
			# E.g: if last = 450mb and first = 550mb, then f = .5
			# If last = 499mb and first = 599mb, then f = .01 
			f = (50000 - last_pressure_before500mb) / abs(first_pressure_after500mb - last_pressure_before500mb)
			#
			#compute the psi fractional addition ---> that pressure jump which 500mb falls in
			fractional_addition = (f * psi[fractional_level_index])
			#
			#compute psi500 ---> sum(psi above 500mb) + fractional component of pressure jump that 500mb falls in
			return np.sum(psi[pressures_above500mb]) + fractional_addition

		psi500 = []
		for i in range(len(lats)):
			for j in range(len(lons)):
				if j == 0:
					hold = []
				hold.append(psi500calc(psi[:,i,j], pressure[:,i,j]))
			psi500.append(hold)
			hold = []

		psi500 = np.array(psi500)
		# print np.shape(psi500) # (96, 144)

		# cs = plt.contour(psi500, 18)
		# plt.clabel(cs, fmt = '%.0f', inline = True)
		# plt.yticks(np.arange(6,96,10),lats[np.arange(6,96,10)])
		# plt.xticks(np.arange(0,144,10), lons[np.arange(0,144,10)])
		# plt.show()

		m = bm(projection = 'cea', llcrnrlat=-90,urcrnrlat=90, llcrnrlon=-60,urcrnrlon=300,resolution='c')
		m.drawcoastlines()
		m.drawparallels(np.arange(-90.,120.,30.))
		m.drawmeridians(np.arange(0.,360.,120.))
		m.drawmapboundary(fill_color='0.3')


		cs = m.contourf(bmlon, bmlat, psi500, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
		m.colorbar(cs, location='bottom', pad="5%")
		plt.show()
'''

		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		# Now we want to calculate the phi-psi-500 value for each longitude
		# This is the latitude at which phi500 changes sign
		# There will be three main instance of this sign change. These three instance represent the ITCZ and the NH and SH Hadley cell edge

		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		for lon_index, lon in enumerate(lons):

			print "Longitude is " + str(lon)

			psihold = psi500[:,lon_index]

			#find the absolute value maximum 
			absmax = max(abs(psihold))
			absmax_index = [n for n, m in enumerate(psihold) if abs(m) == absmax][0]

			#find the value at the absmax index
			absmax_psival = psihold[absmax_index]

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
				first_neg_index = np.argmax(psihold[absmax_index:] < 0)
				
				#calculate the value at and above that index
				first_neg_value = psihold[first_neg_index]
				last_pos_value = psihold[first_neg_index-1]
				
				#find the latitudes at those two points
				first_neg_lat = lats[first_neg_index]
				last_pos_lat = lats[first_neg_index-1]
				
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
				first_neg_index = np.argmax(psihold[list(reversed(range(0,absmax_index)))] < 0)
				
				first_neg_value = psihold[first_neg_index]
				last_pos_value = psihold[first_neg_index+1]
				first_neg_lat = lats[first_neg_index]
				last_pos_lat = lats[first_neg_index+1]
				
				f = abs(last_pos_value) / abs(last_pos_value - first_neg_value)
				
				#going up the data frame (decreasing latitude), so we subtract the fractional latitude change from the last positive latitude
				zerolat = last_pos_lat - f * abs(first_neg_lat - last_pos_lat)
				#this latitude is the center of the ITCZ
				ITCZcenter = zerolat
				
				#################################################################################################
				# now go backwards from the ITCZ center and find the next zero-crossing to find the SH Hadley 
				# edge. Same basic algorithm. 
				#################################################################################################
				first_pos_index = np.argmax(psihold[list(reversed(range(0,first_neg_index)))] > 0)
				
				first_pos_value = psihold[first_pos_index]
				last_neg_value = psihold[first_pos_index+1]
				first_pos_lat = lats[first_pos_index]
				last_neg_lat = lats[first_pos_index+1]
				
				f = abs(last_neg_value) / abs(last_neg_value - first_pos_value)
				
				#going up the data frame so subtract the fractional latitude change from the last positive latitude
				zerolat = last_neg_lat - f * abs(first_neg_lat - last_pos_lat)
				HadleyCellSouthernEdge = zerolat 


			##############################################
			# Do practically the same thing for absmax < 0
			# Just kind of flip the operation around
			##############################################
			if absmax_psival < 0:
				first_pos_index = np.argmax(psihold[list(reversed(range(0,absmax_index)))] > 0)
				
				first_pos_value = psihold[first_pos_index]
				last_neg_value = psihold[first_pos_index+1]
				first_pos_lat = lats[first_pos_index]
				last_neg_lat = lats[first_pos_index+1]
				
				f = abs(last_neg_value) / abs(last_neg_value - first_pos_value)
				
				#going up the data frame so subtract the fractional latitude change from the last positive latitude
				zerolat = last_neg_lat - f * abs(last_neg_lat - first_pos_lat)
				HadleyCellSouthernEdge = zerolat 
				
				################################################################################################
				# now go forwards from the absmax_index and find the next zero-crossing to find the ITCZ center
				################################################################################################
				first_pos_index = np.argmax(psihold[absmax_index:] > 0)
				
				first_pos_value = psihold[first_pos_index]
				last_neg_value = psihold[first_pos_index-1]
				first_pos_lat = lats[first_pos_index]
				last_neg_lat = lats[first_pos_index-1]
				
				f = abs(last_neg_value) / abs(last_neg_value - first_pos_value)
				
				zerolat = last_neg_lat + f * abs(last_neg_lat - first_pos_lat)
				ITCZcenter = zerolat

				################################################################################################
				# now go forwards from the ITCZ center and find the next zero-crossing to find the NH Hadley edge
				################################################################################################
				first_neg_index = np.argmax(psihold[first_pos_index:] < 0)
				
				first_neg_value = psihold[first_neg_index]
				last_pos_value = psihold[first_neg_index-1]
				first_neg_lat = lats[first_neg_index]
				last_pos_lat = lats[first_neg_index-1]
				
				f = abs(first_neg_value) / abs(first_neg_value - last_pos_value)
				
				zerolat = last_pos_lat + f * abs(last_pos_lat - first_neg_lat)
				HadleyCellNorthernEdge = zerolat

			print "Hadley cell northern edge is " + str(HadleyCellNorthernEdge)
			print "Hadley cell southern edge is " + str(HadleyCellSouthernEdge)
			print "\n"

		# Compute some other variables like cell width
		# WidthNorthernCell = abs(HadleyCellNorthernEdge - ITCZcenter)
		# WidthSouthernCell = abs(HadleyCellSouthernEdge - ITCZcenter)

		# # Temporary variable for this instance of the for loop
		# hold = pd.DataFrame({'File' : [fi], 'Time' : [month], 'ITCZcenter' : [ITCZcenter], 'HadleyCellNorthernEdge' : [HadleyCellNorthernEdge], 'HadleyCellSouthernEdge' : [HadleyCellSouthernEdge], 'NorthernCellWidth' : [WidthNorthernCell], 'SouthernCellWidth' : [WidthSouthernCell]})
		
		# # Log data from temp variable into master file
		# InfoDataFrame = InfoDataFrame.append(hold)

		#----------End for loop through time----------#

	#----------End for loop through files----------#

# InfoDataFrame = InfoDataFrame.reset_index(drop=True)
# InfoDataFrame.to_csv("HadleyData.csv")

'''


