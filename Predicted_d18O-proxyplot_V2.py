from pygoda import ncgoda, findClimoFile
import matplotlib.pyplot as plt
import pandas as pd
import os, glob
import numpy as np
from scipy import stats
import datetime

southern_lat = -50
northern_lat = 50
left_lon = 0
right_lon = 355

season = "DJF"
m1 = 1
m2 = 2
m3 = 12

root = '/home/server/student/homes/kniezgod/model_runs/'
cntrldir = root+"F.C5.2deg.wiso.defaultSSTICE_kn002/"
testdir = root+"F.C5.2deg.wiso.obs6kSST_kn003/"

Test_File = glob.glob(testdir + "*" + season + "*")
Control_File = glob.glob(cntrldir + "*" + season + "*")

testdata = ncgoda(Test_File[0])
cntrldata = ncgoda(Control_File[0])

PRECT_control = cntrldata.variable("PRECT")
PRECT_test = testdata.variable("PRECT")
PRECT_units = cntrldata.units
PRECT_long_name = cntrldata.long_name

PRECTd18O_control = cntrldata.PRECT_d18O()
PRECTd18O_test = testdata.PRECT_d18O()
PRECTd18O_units = cntrldata.units
PRECTd18O_long_name = cntrldata.long_name

# Create dates to cycle over files: 

first_model_year = 1
last_model_year = 29

start = datetime.date(year = 1950 + first_model_year, month = m1, day = 1)
end = datetime.date(year = 1950 + last_model_year, month = 12, day = 1)
delta = datetime.timedelta(days = 1)

# Create the timeseries of Precipitation and Precip isotopes
# Test files
d = start
while d <= end:
	year = str(d.year - 1950).zfill(4)
	if d.day != 1:
		# datetime.timedelta doesn't allow months to be the argument so 
		# have to use days and just skip the days != 1
		print "skipping " + year + "-" + d.strftime('%m-%d')
		d += delta
		continue
	#
	if (d.month != m1) & (d.month != m2) & (d.month != m3):
		print "skipping " + year + "-" + d.strftime('%m-%d')
		d += delta
		continue
	#
	# datetime doesn't allow any years < 1950, so we have to do this special thing to make it work
	print("*cam.h0." + year + d.strftime('-%m') + "*")
	f = findClimoFile("*cam.h0." + year + d.strftime('-%m') + "*", directory = testdir)[0]
	data = ncgoda(f)
	PRECT = data.variable('PRECT_H2O')
	d18O = data.PRECT_d18O()
	# Fix d18O for values when PRECT is really small
	toolarge_d18O = d18O > 20
	d18O[toolarge_d18O] = float('NaN')
	toosmall_d18O = d18O < -30
	d18O[toosmall_d18O] = float("NaN")
	if d == start:
		Precip = PRECT
		Isotope = d18O
		lats = data.lat
		lons = data.lon
	else:
		Precip = np.dstack((Precip, PRECT))
		Isotope = np.dstack((Isotope, d18O))
	d += delta

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

# Calculate predicted isotopes using the linear model
Test_predicted_ANNd18O = np.multiply(PRECT_test, Slope) + Intercept


# Control files
d = start
while d <= end:
	if d.day != 1:
		# datetime.timedelta doesn't allow months to be the argument so 
		# have to use days and just skip the days != 1
		print "skipping " + year + "-" + d.strftime('%m-%d')
		d += delta
		continue
	if (d.month != 6) & (d.month != 7) & (d.month != 8):
		print "skipping " + year + "-" + d.strftime('%m-%d')
		d += delta
		continue
	# datetime doesn't allow any years < 1950, so we have to do this special thing to make it work
	year = str(d.year - 1950).zfill(4)
	print("*cam.h0." + year + d.strftime('-%m') + "*")
	f = findClimoFile("*cam.h0." + year + d.strftime('-%m') + "*", directory = cntrldir)[0]
	data = ncgoda(f)
	PRECT = data.variable('PRECT_H2O')
	d18O = data.PRECT_d18O()
	# Fix d18O for values when PRECT is really small
	toolarge_d18O = d18O > 20
	d18O[toolarge_d18O] = float('NaN')
	toosmall_d18O = d18O < -30
	d18O[toosmall_d18O] = float("NaN")
	if d == start:
		Precip = PRECT
		Isotope = d18O
		lats = data.lat
		lons = data.lon
	else:
		Precip = np.dstack((Precip, PRECT))
		Isotope = np.dstack((Isotope, d18O))
	d += delta

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

# Calculate predicted isotopes using the linear model
Control_predicted_ANNd18O = np.multiply(PRECT_control, Slope) + Intercept

# Proxy info
d18O_proxy = pd.DataFrame({'Site' : ["Dg", "Qf", "GB", "RGN", "LL", "KNI", "Bv"], \
    'lat' : [25.17, 17.07, 4, -5.36, -8.32, -15.18, -27.13], \
        'lon' : [108.5, 53.41, 114, -37.44, 120.26, 128.37, -49.09], \
            '6kya_value' : [-8.51, -0.2, -9.65, -6, -6.2, -6.58, -2.34], \
                '1850_value' : [-7.74, -0.81, -9.16, -2.5, -6.07, -7.55, -3.35]})

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


x = range(7)
plt.scatter(x,Test_predicted_ANNd18O[latind,lonind].tolist(), color = 'r', marker = 'o', label = "Amount-effect predicted")
plt.scatter(x,PRECTd18O_test[latind, lonind].tolist(), color = 'b',  marker = '^', label = "Directly simulated")
plt.scatter(x,d18O_proxy['6kya_value'], color = 'k', marker = '*', label = 'Proxy value')
plt.xticks(x, d18O_proxy['Site'])
plt.legend(loc = 'upper center')
plt.xlabel("Speleothem location")
plt.ylabel("d18O")
plt.title("Mid-Holocene, " + season)
# plt.savefig("test_scatter_simluated_d18O-performance.ps")
plt.show()

plt.scatter(x,Control_predicted_ANNd18O[latind,lonind].tolist(), color = 'r', marker = 'o', label = "Amount-effect predicted")
plt.scatter(x,PRECTd18O_control[latind, lonind].tolist(), color = 'b',  marker = '^', label = "Directly simulated")
plt.scatter(x,d18O_proxy['1850_value'], color = 'k', marker = '*', label = 'Proxy value')
plt.xticks(x, d18O_proxy['Site'])
plt.legend(loc = 'best')
plt.xlabel("Speleothem location")
plt.ylabel("d18O")
plt.title("Pre-industrial, " + season)
# plt.savefig("control_scatter_simluated_d18O-performance.ps")
plt.show()

Diff_predicted_ANNd18O = np.subtract(Test_predicted_ANNd18O[latind,lonind].tolist(), Control_predicted_ANNd18O[latind,lonind].tolist())
Diff_simulated_ANNd18O = np.subtract(PRECTd18O_test[latind,lonind].tolist(), PRECTd18O_control[latind,lonind].tolist())
plt.scatter(x,Diff_predicted_ANNd18O, color = 'r', marker = 'o', label = "Amount-effect predicted")
plt.scatter(x,Diff_simulated_ANNd18O, color = 'b',  marker = '^', label = "Directly simulated")
plt.scatter(x,d18O_proxy['6kya-1850_difference'], color = 'k', marker = '*', label = 'Proxy value')
plt.xticks(x, d18O_proxy['Site'])
plt.legend(loc = 'best')
plt.xlabel("Speleothem location")
plt.ylabel("d18O")
plt.title("MH - PI difference, " + season)
# plt.savefig("diff_scatter_simluated_d18O-performance.ps")
plt.show()
