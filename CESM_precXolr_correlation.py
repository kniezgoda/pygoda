#!/Users/kyleniezgoda/anaconda/bin/python
#create OLR~rain correlation matrix for plotting in IDL from CESM model runs

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from openpyxl import load_workbook
from scipy.io import netcdf as nc
import os, sys
import datetime

home=os.path.expanduser("~")

#import cesm netcdf file
read = nc.netcdf_file(home+'/ysdata/1yr_sim/1yr_T31g37_B1850CN.cam.h1.all_data.nc', mode='r')
#extract the OLR data
OLRvals = read.variables['FLNT'][:].copy()
PRECvals = read.variables['PRECT'][:].copy()
LONvals = read.variables['lon'][:].copy()
LATvals = read.variables['lat'][:].copy()
#close netcdf file
read.close()

# #average over the first dimension (time)
# OLRvals = OLRvals.mean(0)

#find the latitude closest to singapore (1.3)
m = min(abs(LATvals-1.3))
lat = [n for n, x in enumerate(abs(LATvals-1.3)) if x == m][0]
#find the long closest to singapore (103.5)
m = min(abs(LONvals-103.5))
lon = [n for n, x in enumerate(abs(LONvals-103.5)) if x == m][0]

#extract a 1-dim precip vector from PRECvals using the singapore lat and lon
sing_prec = PRECvals[:,lat,lon]

#for loop through the cells of OLRvals, compute correlation
#coefficient for each cell vector (n = 30 for each cell)

#initiate empty array
corr_array = np.empty([48, 96])

for i in range(OLRvals.shape[1]) :
	for j in range(OLRvals.shape[2]) :
		corr_array[i,j] = np.corrcoef(sing_prec, OLRvals[:,i,j])[0,1]

#write the correlation array to a csv file
#np.savetxt('/Volumes/NiezFiles/CESM_files/2month_h0mnth_h1day_CESMsimulation/test.csv', corr_array, delimiter=',')

#code for netcdf file creation and export:
#open this netcdf file in panoply (ncview/anything else) and look at correlation variable. 
#Panoply will automatically add a geo layer which is nice
CorrelationMatrix = nc.netcdf_file(home+'/ysdata/1yr_sim/correlation_1yr_alldata.nc', mode='w')
nclat = CorrelationMatrix.createDimension('lat', len(LATvals))
nclon = CorrelationMatrix.createDimension('lon', len(LONvals))
latitudes = CorrelationMatrix.createVariable('latitude', 'd', ('lat',))
longitudes = CorrelationMatrix.createVariable('longitude', 'd', ('lon',))
prec_vs_olr_correlation = CorrelationMatrix.createVariable('prec_vs_olr_correlation', 'd', ('lat', 'lon'))
latitudes[:] = LATvals
longitudes[:] = LONvals
prec_vs_olr_correlation[:] = corr_array
CorrelationMatrix.close()


#test plotting, usually would use IDL for this but not sure about IDL access on OSU system yet
# im = pd.read_csv('/Volumes/NiezFiles/CESM_files/2month_h0mnth_h1day_CESMsimulation/test.csv', header=None) 
# fig = plt.figure(1)
# plt.imshow(im, interpolation='bilinear')
# plt.show()
