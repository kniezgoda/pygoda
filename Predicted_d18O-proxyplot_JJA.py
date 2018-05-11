import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.basemap import Basemap as bm
import os, glob
from netCDF4 import Dataset
import numpy as np
from scipy import stats

southern_lat = -50
northern_lat = 50
left_lon = 0
right_lon = 355

root = '/home/server/student/homes/kniezgod/model_runs/'
cntrldir = root+"F.C5.2deg.wiso.defaultSSTICE_kn002/"
testdir = root+"F.C5.2deg.wiso.obs6kSST_kn003/"

Test_JJA_file = glob.glob(testdir+"*JJA*")
Control_JJA_file = glob.glob(cntrldir+"*JJA*")

testdata = Dataset(Test_JJA_file[0], mode = 'r')
cntrldata = Dataset(Control_JJA_file[0], mode = 'r')

PRECT_control = cntrldata.variables['PRECT'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
PRECT_test = testdata.variables['PRECT'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
PRECT_units = "mm/day"
PRECT_long_name = "Total precipitation rate"

PRECT_H2O_control = cntrldata.variables['PRECT_H2O'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
PRECT_H2O_test = testdata.variables['PRECT_H2O'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
PRECT_H2O_units = "mm/day"
PRECT_H2O_long_name = cntrldata.variables['PRECT_H2O'].long_name

PRECT_H218O_control = cntrldata.variables['PRECT_H218O'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
PRECT_H218O_test = testdata.variables['PRECT_H218O'][0,:,:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
PRECT_H218O_units = "mm/day"
PRECT_H218O_long_name = cntrldata.variables['PRECT_H218O'].long_name

PRECT_d18O_control = (PRECT_H218O_control/PRECT_H2O_control - 1) * 1000
PRECT_d18O_test = (PRECT_H218O_test/PRECT_H2O_test - 1) * 1000
PRECT_d18O_units = "delta 18O (permil)"
PRECT_d18O_long_name = "delta 18O for PRECT"

# Test 
flist = glob.glob(testdir+"*cam.h0.????-06*")
flist += glob.glob(testdir+"*cam.h0.????-07*")
flist += glob.glob(testdir+"*cam.h0.????-08*")

# Create time series
for n, f in enumerate(flist):
	data = Dataset(f, mode = 'r')
	PRECT = data.variables['PRECT_H2O'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_H218O = data.variables['PRECT_H218O'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	d18O = (PRECT_H218O / PRECT - 1) * 1000
	# Fix d18O for values when PRECT is really small
	toolarge_d18O = d18O > 20
	d18O[toolarge_d18O] = float('NaN')
	toosmall_d18O = d18O < -30
	d18O[toosmall_d18O] = float("NaN")
	if n == 0:
		Precip = PRECT
		Isotope = d18O
		lats = data.variables['lat'][:]
		lons = data.variables['lon'][:]
	else:
		Precip = np.dstack((Precip, PRECT))
		Isotope = np.dstack((Isotope, d18O))


# Create linear model for each grid cell
Slope = np.zeros((Precip.shape[0], Precip.shape[1]))
Intercept = np.zeros((Precip.shape[0], Precip.shape[1]))
for i in range(Precip.shape[0]):
	for j in range(Precip.shape[1]):
		mask = ~np.isnan(Isotope[i,j,:])
		if len(Isotope[i,j,mask]) is not 0:
			Slope[i,j] = stats.linregress(Precip[i,j,mask], Isotope[i,j,mask])[0]
			Intercept[i,j] = stats.linregress(Precip[i,j,mask], Isotope[i,j,mask])[1]
		else:
			Slope[i,j] = float('NaN')
			Intercept[i,j] = float('NaN')


# Predict d18O using linear model and JJA precip
Test_predicted_JJAd18O = np.multiply(PRECT_test, Slope) + Intercept


# Control
flist = glob.glob(cntrldir+"*cam.h0.????-06*")
flist += glob.glob(cntrldir+"*cam.h0.????-07*")
flist += glob.glob(cntrldir+"*cam.h0.????-08*")
for n, f in enumerate(flist):
	data = Dataset(f, mode = 'r')
	PRECT = data.variables['PRECT_H2O'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	PRECT_H218O = data.variables['PRECT_H218O'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
	d18O = (PRECT_H218O / PRECT - 1) * 1000
	# Fix d18O for values when PRECT is really small
	toolarge_d18O = d18O > 20
	d18O[toolarge_d18O] = float('NaN')
	toosmall_d18O = d18O < -30
	d18O[toosmall_d18O] = float('NaN')
	if n == 0:
		Precip = PRECT[:]
		Isotope = d18O[:]
		lats = data.variables['lat'][:]
		lons = data.variables['lon'][:]
	else:
		Precip = np.dstack((Precip, PRECT))
		Isotope = np.dstack((Isotope, d18O))


Slope = np.zeros((Precip.shape[0], Precip.shape[1]))
Intercept = np.zeros((Precip.shape[0], Precip.shape[1]))
for i in range(Precip.shape[0]):
	for j in range(Precip.shape[1]):
		mask = ~np.isnan(Isotope[i,j,:])
		if len(Isotope[i,j,mask]) is not 0:
			Slope[i,j] = stats.linregress(Precip[i,j,mask], Isotope[i,j,mask])[0]
			Intercept[i,j] = stats.linregress(Precip[i,j,mask], Isotope[i,j,mask])[1]
		else:
			Slope[i,j] = float('NaN')
			Intercept[i,j] = float('NaN')


# Predict d18O using linear model and JJA precip
Control_predicted_JJAd18O = np.multiply(PRECT_control, Slope) + Intercept


# Proxy plotting
d18O_proxy = pd.DataFrame({'Site' : ["Dg", "Qf", "GB", "RGN", "LL",  "KNI",  "Bv"], \
    'lat' : [25.17,  17.07, 4, -5.36,  -8.32,  -15.18, -27.13], \
        'lon' : [108.5, 53.41,114, -37.44,  120.26,  128.37, -49.09], \
            '6kya_value' : [-8.51, -0.2,-9.65, -6,  -6.2,  -6.58, -2.34], \
                '1850_value' : [-7.74, -0.81,-9.16, -2.5, -6.07,  -7.55, -3.35]})

d18O_proxy['6kya-1850_difference'] = d18O_proxy['6kya_value'] - d18O_proxy['1850_value']

proxy_x, proxy_y = np.meshgrid(d18O_proxy['lon'], d18O_proxy['lat'])
num_proxies = d18O_proxy.shape[0]
bmlon, bmlat = np.meshgrid(lons, lats)

latind = []
for l in d18O_proxy['lat']:
	hold = [abs(LAT - l) for LAT in lats]
	holdmin = min(hold)
	minlatind = [n for n, x in enumerate(hold) if x == holdmin]
	latind.append(minlatind[0])


lonind = []
for l in d18O_proxy['lon']:
	hold = [abs(LON - l) for LON in lons]
	holdmin = min(hold)
	minlonind = [n for n, x in enumerate(hold) if x == holdmin]
	lonind.append(minlonind[0])

'''
proxy_x, proxy_y = np.meshgrid(d18O_proxy['lon'], d18O_proxy['lat'])
bmlon, bmlat = np.meshgrid(lons, lats)

# Plot test data
plt.subplot(3,1,1)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
min_1850 = -25
max_1850 = 0
clevs = np.linspace(min_1850,max_1850,21)
cs = m.contourf(bmlon, bmlat, Test_predicted_d18O, 20, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
proxy = m.scatter(proxy_x[0,:],proxy_y[:,0], c = d18O_proxy['6kya_value'], vmin = min_1850, vmax = max_1850, s= 100, cmap=plt.cm.RdBu_r, latlon=True)

# Plot control data
plt.subplot(3,1,2)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
clevs = np.linspace(min_1850,max_1850,21)
cs = m.contourf(bmlon, bmlat, Control_predicted_d18O, 20, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
proxy = m.scatter(proxy_x[0,:],proxy_y[:,0], c = d18O_proxy['1850_value'], vmin = min_1850, vmax = max_1850, s= 100, cmap=plt.cm.RdBu_r, latlon=True)


# Plot difference data 
plt.subplot(3,1,3)
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
min_diff = -4
max_diff = 4
clevs = np.linspace(min_diff,max_diff,21)
cs = m.contourf(bmlon, bmlat, Test_predicted_d18O-Control_predicted_d18O, clevs, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")
proxy = m.scatter(proxy_x[0,:],proxy_y[:,0], c = d18O_proxy['6kya-1850_difference'], vmin = min_diff, vmax = max_diff, s= 100, cmap=plt.cm.RdBu_r, latlon=True)


# Make things pretty
plt.subplots_adjust(hspace = .15)

plt.show()

'''
PRECT_d18O_test_JJA = PRECT_d18O_test
PRECT_d18O_control_JJA=PRECT_d18O_control


x = range(7)
plt.scatter(x,Test_predicted_ANNd18O[latind,lonind].tolist(), color = 'b', marker = 'o', label = "Amount-effect predicted, ANN")
plt.scatter(x,PRECT_d18O_test[latind, lonind].tolist(), color = 'b',  marker = '^', label = "Directly simulated, ANN")
plt.scatter(x,Test_predicted_JJAd18O[latind,lonind].tolist(), color = 'r', marker = 'o', label = "Amount-effect predicted, JJA")
plt.scatter(x,PRECT_d18O_test_JJA[latind, lonind].tolist(), color = 'r',  marker = '^', label = "Directly simulated, JJA")
plt.scatter(x,d18O_proxy['6kya_value'], color = 'k', marker = '*', label = 'Proxy value')
plt.xticks(x, d18O_proxy['Site'])
plt.legend(loc = 'upper center')
plt.xlabel("Speleothem location")
plt.ylabel("d18O")
plt.title("Mid-Holocene")
plt.savefig("JJA_test_scatter_simluated_d18O-performance.ps")
plt.show()

plt.scatter(x,Control_predicted_ANNd18O[latind,lonind].tolist(), color = 'b', marker = 'o', label = "Amount-effect predicted, ANN")
plt.scatter(x,PRECT_d18O_control[latind, lonind].tolist(), color = 'b',  marker = '^', label = "Directly simulated, ANN")
plt.scatter(x,Control_predicted_JJAd18O[latind,lonind].tolist(), color = 'r', marker = 'o', label = "Amount-effect predicted")
plt.scatter(x,PRECT_d18O_control_JJA[latind, lonind].tolist(), color = 'r',  marker = '^', label = "Directly simulated")
plt.scatter(x,d18O_proxy['1850_value'], color = 'k', marker = '*', label = 'Proxy value')
plt.xticks(x, d18O_proxy['Site'])
#plt.legend(loc = 'lower right')
plt.xlabel("Speleothem location")
plt.ylabel("d18O")
plt.title("Pre-industrial")
plt.savefig("JJA_control_scatter_simluated_d18O-performance.ps")
plt.show()

Diff_predicted_ANNd18O = np.subtract(Test_predicted_ANNd18O[latind,lonind].tolist(), Control_predicted_ANNd18O[latind,lonind].tolist())
Diff_simulated_ANNd18O = np.subtract(PRECT_d18O_test[latind,lonind].tolist(), PRECT_d18O_control[latind,lonind].tolist())
Diff_predicted_JJAd18O = np.subtract(Test_predicted_JJAd18O[latind,lonind].tolist(), Control_predicted_JJAd18O[latind,lonind].tolist())
Diff_simulated_JJAd18O = np.subtract(PRECT_d18O_test_JJA[latind,lonind].tolist(), PRECT_d18O_control_JJA[latind,lonind].tolist())
plt.scatter(x,Diff_predicted_ANNd18O, color = 'b', marker = 'o', label = "Amount-effect predicted, ANN")
plt.scatter(x,Diff_simulated_ANNd18O, color = 'b',  marker = '^', label = "Directly simulated, ANN")
plt.scatter(x,Diff_predicted_JJAd18O, color = 'r', marker = 'o', label = "Amount-effect predicted, JJA")
plt.scatter(x,Diff_simulated_JJAd18O, color = 'r',  marker = '^', label = "Directly simulated, JJA")
plt.scatter(x,d18O_proxy['6kya-1850_difference'], color = 'k', marker = '*', label = 'Proxy value')
plt.xticks(x, d18O_proxy['Site'])
#plt.legend(loc = 'lower right')
plt.xlabel("Speleothem location")
plt.ylabel("d18O")
plt.title("MH - PI difference")
plt.savefig("JJA_diff_scatter_simluated_d18O-performance.ps")
plt.show()


'''
# Plot as scatter plot
plt.scatter(Test_predicted_JJAd18O[latind,lonind].tolist(), d18O_proxy['6kya_value'], color = 'r', marker = "o", label = "Amount-effect predicted")
plt.scatter(Test_simulated_JJAd18O[latind, lonind].tolist(), d18O_proxy['6kya_value'],color = 'b',  marker = '^', label = "Directly simulated")
plt.xlabel("Model d18O")
plt.ylabel("Proxy d18O")
plt.legend(loc = 'lower right')
plt.title("Mid-Holocene")
#plt.savefig("test_scatter_simluated_d18O-performance.ps")
plt.show()

plt.scatter(Control_predicted_JJAd18O[latind,lonind].tolist(), d18O_proxy['1850_value'], color = "r", marker = 'o', label = "Amount-effect predicted")
plt.scatter(PRECT_d18O_control[latind, lonind].tolist(), d18O_proxy['1850_value'], color = 'b', marker = '^', label = "Directly simulated")
plt.xlabel("Model d18O")
plt.ylabel("Proxy d18O")
plt.legend(loc = 'lower right')
plt.title("Pre-industrial")
#plt.savefig("control_scatter_simluated_d18O-performance.ps")
plt.show()

Diff_predicted_JJAd18O = np.subtract(Test_predicted_JJAd18O[latind,lonind].tolist(), Control_predicted_JJAd18O[latind,lonind].tolist())
Diff_simulated_JJAd18O = np.subtract(Test_simulated_JJAd18O[latind,lonind].tolist(), Control_simulated_JJAd18O[latind,lonind].tolist())
plt.scatter(Diff_predicted_JJAd18O, d18O_proxy['6kya-1850_difference'], color = 'r', marker = "o", label = "Amount-effect predicted")
plt.scatter(Diff_simulated_JJAd18O, d18O_proxy['6kya-1850_difference'],color = 'b',  marker = '^', label = "Directly simulated")
plt.xlabel("Model d18O")
plt.ylabel("Proxy d18O")
plt.legend(loc = 'lower left')
plt.title("MH - PI difference")
#plt.savefig("diff_scatter_simluated_d18O-performance.ps")
plt.show()
'''
