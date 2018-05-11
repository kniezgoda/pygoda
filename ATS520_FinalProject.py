#!/Users/kyleniezgoda/anaconda/bin/python
import os, glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame
from UsefulNetcdfFunctions import find_indices

# Set the root directory where the annually averaged files are located
root = "/data/niezgodk/"
nyears = 11

colnames = ['Year','Global', 'India', 'MaritimeContinent', 'EquatorialPacific', 'NorthPole', 'SouthPole','USA']
T2_DataFrame = PRECT_DataFrame = QFLX_DataFrame = DataFrame(columns = colnames)

# Box is the latitude and longitudes you're looking for
# lats and lons are the lists of latitude and longitude values
# def find_indices(box, lats, lons):
#     r = []
#     for N, b in enumerate(box):
#         if N < 2:
#             r.append(abs(lat-b).tolist().index(min(abs(lat-b))))
#         if N >= 2:
#             r.append(abs(lon-b).tolist().index(min(abs(lon-b))))
#     return r

for i in range(nyears):
    # Set a variable for the year with a nice name
    year = i+1
    
    # Locate the climo file for this year
    f = glob.glob(root+"*yr"+str(year)+"climo*")[0]
    
    # Open the netcdf file
    data = Dataset(f, "r")
    
    # Read in the variable data
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    t2 = data.variables['T'][:,-1,:,:].squeeze() # -1st level = ground = ~2m
    prect = data.variables['PRECT'][:,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day 
    qflx = data.variables['QFLX'][:,:,:].squeeze() * 1000 # Convert to g/m2/s 

    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #
    # Global averages:
    '''
    Steps:
        1) zonally average   
        2) average zones, weight by cos(lat)
    '''
    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- # 
    
    # Average zonally
    global_t2_zonal = np.mean(t2, axis = 1)
    global_prect_zonal = np.mean(prect, axis = 1)
    global_qflx_zonal = np.mean(qflx, axis = 1)

    # Average all the latitudes, weighting by cos(lat)
    wgt = np.cos(np.radians(lat))
    global_t2_avg = np.average(global_t2_zonal, weights = wgt)
    global_prect_avg = np.average(global_prect_zonal, weights = wgt)
    global_qflx_avg = np.average(global_qflx_zonal, weights = wgt)
    
    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #
    # Regional averages:
    '''
    Pretty sure I don't have to weight average here
    Safe to assume that the regional boxes are small enough 
    that differences in latitude won't effec the regional average 
    temperature significantly
    
    Steps: 
        1) extract data for the regions I am interested in 
            a) box = [bottom_lat, top_lat, left_lon, right_lon]
        2) average the values
    '''
    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #
    
    # India
    region = [5, 35, 65, 90]
    box = find_indices(region, lat, lon)
    # Extract the data for the region and average it all
    india_t2 = t2[box[0]:box[1], box[2]:box[3]]
    india_t2_avg = np.average(india_t2)
    india_prect = prect[box[0]:box[1], box[2]:box[3]]
    india_prect_avg = np.average(india_prect)
    india_qflx = qflx[box[0]:box[1], box[2]:box[3]]
    india_qflx_avg = np.average(india_qflx)
    
    # Maritime Continent
    region = [-10, 10, 90, 150]
    box = find_indices(region, lat, lon)
    mc_t2 = t2[box[0]:box[1], box[2]:box[3]]
    mc_t2_avg = np.average(mc_t2)
    mc_prect = prect[box[0]:box[1], box[2]:box[3]]
    mc_prect_avg = np.average(mc_prect)
    mc_qflx = qflx[box[0]:box[1], box[2]:box[3]]
    mc_qflx_avg = np.average(mc_qflx)
    
    # Equatorial Pacific
    region = [-5, 5, 160, 270]
    box = find_indices(region, lat, lon)
    ep_t2 = t2[box[0]:box[1], box[2]:box[3]]
    ep_t2_avg = np.average(ep_t2)
    ep_prect = prect[box[0]:box[1], box[2]:box[3]]
    ep_prect_avg = np.average(ep_prect)
    ep_qflx = qflx[box[0]:box[1], box[2]:box[3]]
    ep_qflx_avg = np.average(ep_qflx)
    
    # North Pole
    region = [80, 90, 0, 360]
    box = find_indices(region, lat, lon)
    np_t2 = t2[box[0]:box[1], box[2]:box[3]]
    np_t2_avg = np.average(np_t2)
    np_prect = prect[box[0]:box[1], box[2]:box[3]]
    np_prect_avg = np.average(np_prect)
    np_qflx = qflx[box[0]:box[1], box[2]:box[3]]
    np_qflx_avg = np.average(np_qflx)
    
    # South Pole
    region = [-90, 70, 0, 360]
    box = find_indices(region, lat, lon)
    sp_t2 = t2[box[0]:box[1], box[2]:box[3]]
    sp_t2_avg = np.average(sp_t2)
    sp_prect = prect[box[0]:box[1], box[2]:box[3]]
    sp_prect_avg = np.average(sp_prect)
    sp_qflx = qflx[box[0]:box[1], box[2]:box[3]]
    sp_qflx_avg = np.average(sp_qflx)
    
    # USA, just because
    region = [25, 50, 230, 290]
    box = find_indices(region, lat, lon)
    usa_t2 = t2[box[0]:box[1], box[2]:box[3]]
    usa_t2_avg = np.average(usa_t2)
    usa_prect = prect[box[0]:box[1], box[2]:box[3]]
    usa_prect_avg = np.average(usa_prect)
    usa_qflx = qflx[box[0]:box[1], box[2]:box[3]]
    usa_qflx_avg = np.average(usa_qflx)
    
    
    # Wrap everything up into a temporary DataFrame
    t2_tempdf = DataFrame({'Year' : [int(year)], \
                           'Global' : [global_t2_avg], \
                           'India' : [india_t2_avg], \
                           'MaritimeContinent' : [mc_t2_avg], \
                           'EquatorialPacific' : [ep_t2_avg], \
                           'NorthPole' : [np_t2_avg], \
                           'SouthPole' : [sp_t2_avg], \
                           'USA' : [usa_t2_avg]})
    prect_tempdf = DataFrame({'Year' : [int(year)], \
                           'Global' : [global_prect_avg], \
                           'India' : [india_prect_avg], \
                           'MaritimeContinent' : [mc_prect_avg], \
                           'EquatorialPacific' : [ep_prect_avg], \
                           'NorthPole' : [np_prect_avg], \
                           'SouthPole' : [sp_prect_avg], \
                           'USA' : [usa_prect_avg]})
    qflx_tempdf = DataFrame({'Year' : [int(year)], \
                           'Global' : [global_qflx_avg], \
                           'India' : [india_qflx_avg], \
                           'MaritimeContinent' : [mc_qflx_avg], \
                           'EquatorialPacific' : [ep_qflx_avg], \
                           'NorthPole' : [np_qflx_avg], \
                           'SouthPole' : [sp_qflx_avg], \
                           'USA' : [usa_qflx_avg]})
    
    # Stack it onto the bottom of the master DataFrame
    T2_DataFrame = T2_DataFrame.append(t2_tempdf)
    PRECT_DataFrame = PRECT_DataFrame.append(prect_tempdf)
    QFLX_DataFrame = QFLX_DataFrame.append(qflx_tempdf)
    
# Plot temps
plt.plot(T2_DataFrame['Year'], T2_DataFrame['Global']/max(T2_DataFrame['Global']), label = "Global")
plt.plot(T2_DataFrame['Year'], T2_DataFrame['USA']/max(T2_DataFrame['USA']), label = "USA")
plt.plot(T2_DataFrame['Year'], T2_DataFrame['India']/max(T2_DataFrame['India']), label = "India")
plt.plot(T2_DataFrame['Year'], T2_DataFrame['MaritimeContinent']/max(T2_DataFrame['MaritimeContinent']), label = "Maritime Continent")
plt.plot(T2_DataFrame['Year'], T2_DataFrame['NorthPole']/max(T2_DataFrame['NorthPole']), label = "North Pole")
plt.plot(T2_DataFrame['Year'], T2_DataFrame['SouthPole']/max(T2_DataFrame['SouthPole']), label = "South Pole")
plt.plot(T2_DataFrame['Year'], T2_DataFrame['EquatorialPacific']/max(T2_DataFrame['EquatorialPacific']), label = "Equatorial Pacific")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Model Year")
plt.ylabel("Variable normalized by maximum value")
plt.xlim(0,nyears+1)
plt.title("Temperature")
plt.show()

# Precipitation
plt.plot(PRECT_DataFrame['Year'], PRECT_DataFrame['Global']/max(PRECT_DataFrame['Global']), label = "Global")
plt.plot(PRECT_DataFrame['Year'], PRECT_DataFrame['USA']/max(PRECT_DataFrame['USA']), label = "USA")
plt.plot(PRECT_DataFrame['Year'], PRECT_DataFrame['India']/max(PRECT_DataFrame['India']), label = "India")
plt.plot(PRECT_DataFrame['Year'], PRECT_DataFrame['MaritimeContinent']/max(PRECT_DataFrame['MaritimeContinent']), label = "Maritime Continent")
plt.plot(PRECT_DataFrame['Year'], PRECT_DataFrame['NorthPole']/max(PRECT_DataFrame['NorthPole']), label = "North Pole")
plt.plot(PRECT_DataFrame['Year'], PRECT_DataFrame['SouthPole']/max(PRECT_DataFrame['SouthPole']), label = "South Pole")
plt.plot(PRECT_DataFrame['Year'], PRECT_DataFrame['EquatorialPacific']/max(PRECT_DataFrame['EquatorialPacific']), label = "Equatorial Pacific")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Model Year")
plt.ylabel("Variable normalized by maximum value")
plt.xlim(0,nyears+1)
plt.title("Precipitation")
plt.show()

# Evaporation
plt.plot(QFLX_DataFrame['Year'], QFLX_DataFrame['Global']/max(QFLX_DataFrame['Global']), label = "Global")
plt.plot(QFLX_DataFrame['Year'], QFLX_DataFrame['USA']/max(QFLX_DataFrame['USA']), label = "USA")
plt.plot(QFLX_DataFrame['Year'], QFLX_DataFrame['India']/max(QFLX_DataFrame['India']), label = "India")
plt.plot(QFLX_DataFrame['Year'], QFLX_DataFrame['MaritimeContinent']/max(QFLX_DataFrame['MaritimeContinent']), label = "Maritime Continent")
plt.plot(QFLX_DataFrame['Year'], QFLX_DataFrame['NorthPole']/max(QFLX_DataFrame['NorthPole']), label = "North Pole")
plt.plot(QFLX_DataFrame['Year'], QFLX_DataFrame['SouthPole']/max(QFLX_DataFrame['SouthPole']), label = "South Pole")
plt.plot(QFLX_DataFrame['Year'], QFLX_DataFrame['EquatorialPacific']/max(QFLX_DataFrame['EquatorialPacific']), label = "Equatorial Pacific")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Model Year")
plt.ylabel("Variable normalized by maximum value")
plt.xlim(0,nyears+1)
plt.title("Evaporation")
plt.show()