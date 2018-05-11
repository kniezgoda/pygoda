#!/Users/kyleniezgoda/anaconda/bin/python

'''
This code calculates omega (vertical pressure velocity) at a specific pressure level. 
A basemap is created and colored contours of omega are plotted.
Can also plot wind barbs if requested.

Author: Kyle Niezgoda
Date: May 2, 2017
'''

import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import shiftgrid 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- User-entered data here --- #
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set the input directory
# directory = "/Users/kyleniezgoda/Yellowstone_Runs/model_runs/F2deg.wiso.defaultSSTICE_kn001/"
directory = "/Users/kyleniezgoda/Yellowstone_Runs/model_runs/F.C5.2deg.wiso.obs6kSST_kn003/"

# Set your input file here
# Mid-holocene data
fi = "F.C5.2deg.wiso.obs6kSST_kn003_JJA_climo.nc"

# Pre-industrial era data
# fi = "F2deg.wiso.defaultSSTICE_kn001_JJA_climo.nc"
# fi = "F2deg.wiso.defaultSSTICE_kn001_DJF_climo.nc"
# fi = "F2deg.wiso.defaultSSTICE_kn001_ANN_climo.nc"

# Set the pressure at which you wish to calculate vertical velocity 
p = 50000 # Pascals

# Do you wish to plot wind barbs with the omega contours? (boolean)
plot_wind_barbs = True
wind_barb_density = 3 # int, larger numbers indicate lower density (reasonable numbers are 1 - 8)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Main code block starts here --- #
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Read the netcdf file and extract data
data = Dataset(directory+fi, mode='r')

lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
levs = data.variables['lev'][:] * 100
OMEGA = data.variables['OMEGA'][0,:,:,:]
A = data.variables['hyam'][:]
B = data.variables['hybm'][:]
PS = data.variables['PS'][0,:,:]
U = data.variables['U'][0,:,:,:]
V = data.variables['V'][0,:,:,:]

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Pressure calcs --- #
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Pressure is not as simple as looking at levs because CAM has terrain following pressure scheme
# Pressure = A*P0 + B*PS
# P0 = 100000, PS = surface pressure from cam

#creates 3-d pressure array 
pressure = []
for n in range(len(B)):
	pressure.append(np.array(A[n] * 100000 + B[n] * PS))
pressure = np.array(pressure)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Algorithm for finding 500mb vertical wind speed --- #
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 1) loop through lats and lons in the pressure array
# 2) find the first pressure level above and below 500mb
# 3) calculate the fraction, f, corresponding to where 500mb falls within the range [p_above, p_below]
# 4) apply that fraction to the omega value at that lat/lon column

# Initialize master array
# We will stack each latitude band onto the bottom of this
# Will remove the top row (all zeros) at the end
Omega = np.zeros(len(lons))
u = np.zeros(len(lons))
v = np.zeros(len(lons))

for i in range(len(lats)):
	# Initialize empty array to hold omega, u, and v values for this latitude band
	omega_array_by_latzone_hold = []
	u_array_by_latzone_hold = []
	v_array_by_latzone_hold = []

	for j in range(len(lons)):
		# Extract pressure and omega arrays for this lat/lon coordinate
		pressure_hold = pressure[:,i,j]
		omega_hold = OMEGA[:,i,j]
		u_hold = U[:,i,j]
		v_hold = V[:,i,j]

		# Here, "below" and "above" mean in space, thus 60000 is "below" 50000 even though 6 > 5
		pressures_above500mb = pressure_hold - p < 0
		pressures_below500mb = pressure_hold - p > 0
		
		# Check for the case when all pressures are above our value for p
		# This can happen at the poles or over mountain ranges where PS values are < our p value
		# Might want to set this to just return the value at the ground instead of NaN
		if sum(pressures_below500mb) == 0:
			omega_array_by_latzone_hold.append(np.nan)
			u_array_by_latzone_hold.append(np.nan)
			v_array_by_latzone_hold.append(np.nan)

		else:
			specific_pressure_above = pressure_hold[pressures_above500mb][-1]
			specific_pressure_below = pressure_hold[pressures_below500mb][0]

			# Computes the fraction of the pressurelevel that 500 falls in
			# We will integrate from top down until we hit the last pressure level that is less than 500mb
			# At that point, we use the fraction of the next pressure jump that will correspond to 500mb
			# E.g: if specific_pressure_above = 450mb and specific_pressure_below = 550mb, then f = .5
			# If specific_pressure_above = 499mb and specific_pressure_below = 599mb, then f = .01 
			f = (p - specific_pressure_above) / abs(specific_pressure_below - specific_pressure_above)

			specific_omega_above = omega_hold[pressures_above500mb][-1]
			specific_omega_below = omega_hold[pressures_below500mb][0]

			specific_u_above = u_hold[pressures_above500mb][-1]
			specific_u_below = u_hold[pressures_below500mb][0]

			specific_v_above = v_hold[pressures_above500mb][-1]
			specific_v_below = v_hold[pressures_below500mb][0]

			# Now we know f, want to solve for omega, v, and u
			omega500 = f * abs(specific_omega_below - specific_omega_above) + specific_omega_above
			u500 = f * abs(specific_u_below - specific_u_above) + specific_u_above
			v500 = f * abs(specific_v_below - specific_v_above) + specific_v_above
			
			# Add to the latitude band hold array
			omega_array_by_latzone_hold.append(omega500)
			u_array_by_latzone_hold.append(u500)
			v_array_by_latzone_hold.append(v500)

	# Add to the master array
	Omega = np.vstack([Omega, omega_array_by_latzone_hold])
	u = np.vstack([u, u_array_by_latzone_hold])
	v = np.vstack([v, v_array_by_latzone_hold])

# Get rid of the initialzation row (full of zeros)
Omega = Omega[1:,:]
u = u[1:,:]
v = v[1:,:]

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Main code block ends here --- #
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#####################################################################################################################################################################################

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Map creation code:
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bmlon, bmlat = np.meshgrid(lons, lats)

m = bm(projection = 'cea', llcrnrlat=-15,urcrnrlat=55, llcrnrlon=-35,urcrnrlon=120,resolution='c')
# m = bm(projection = 'cea', llcrnrlat=-90,urcrnrlat=90, llcrnrlon=-60,urcrnrlon=300,resolution='c')

m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,15.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='0.3')

clevs = np.linspace(-.2, .2, 13)
cs = m.contourf(bmlon, bmlat, Omega, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='bottom', pad="5%")
cbar.set_label("Vertical pressure velocity (Pa/s)")

if plot_wind_barbs:
	colsby4 = np.arange(0, bmlon.shape[1], wind_barb_density)
	rowsby4 = np.arange(0, bmlat.shape[0], wind_barb_density)

	barb_lon = bmlon[:,colsby4]
	barb_lon = barb_lon[rowsby4,:]
	
	barb_lat = bmlat[:,colsby4]
	barb_lat = barb_lat[rowsby4,:]	

	barb_u = u[rowsby4,:]
	barb_u = barb_u[:,colsby4]

	barb_v = v[rowsby4,:]
	barb_v = barb_v[:,colsby4]

	m.quiver(barb_lon, barb_lat, barb_u, barb_v, cmap=plt.cm.autumn, latlon = True)

plt.title("Vertical pressure velocity at " + str(p/100) + "mb")
plt.show()
