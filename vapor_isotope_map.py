#Python code for basic isotope graphs
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap as bm
from mpl_toolkits.basemap import shiftgrid 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Set the input file
directory = "/Users/kyleniezgoda/Yellowstone_Runs/AMWG_Diagnostics/"
fi = "isoCAM5CLM4_2deg_wdiddledobs6kSSTs_DJF_climo.nc"
fn = os.path.splitext(fi)[0]
#Set level ("ground", "850hPa", "500hPa", 200hPa)
level = "ground"

#Define what to plot ("d18OV", "dDV", "temp")
plotwhat = "d18OV"

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Read the netcdf file
data = Dataset(directory+fi, mode='r')

#Extract lat, lon, and lev data
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
bmlon, bmlat = np.meshgrid(lons, lats)

levs = data.variables['lev'][:]
lev850 = [n for n,i in enumerate([abs(l - 850) for l in levs]) if i == min([abs(l - 850) for l in levs])]
lev500 = [n for n,i in enumerate([abs(l - 500) for l in levs]) if i == min([abs(l - 500) for l in levs])]
levground = -1

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#User set level
if level == "ground":
	lev = levground
elif level == "850hPa":
	lev = lev850
elif level == "500hPa":
	lev = lev500
elif level == "200hPa":
	lev = lev200
else:
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "LEVEL NOT SET"
	print "Ground level being used, title of plot may be misleading!"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	lev = levground

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Extract data. Variable names are typical CAM5 output with isotopes enabled
temp = data.variables['T'][0,lev,:,:].squeeze()
HDOV = data.variables['HDOV'][0,lev,:,:].squeeze()
H218OV = data.variables['H218OV'][0,lev,:,:].squeeze()
H2OV = data.variables['H2OV'][0,lev,:,:].squeeze()

#Compute d18O and dD in delta values
dDV = (HDOV/H2OV - 1) * 1000
d18OV = (H218OV/H2OV - 1) * 1000

#Convert temp data to celcius
temp = temp - 273

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Creation of maps

#Regional box edge definitions
# box = pd.read_excel("/Users/kyleniezgoda/Documents/Noone_Group/Protocols/RegionalBoxes.xlsx")

#Set variable to plot here
if plotwhat == "d18OV":
	plot = d18OV
elif plotwhat == "dDV":
	plot = dDV
elif plotwhat == "temp":
	plot = temp
else:
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "PLOTWHAT NOT SET"
	print "d18OV will be plotted, title of plot may be misleading!"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	print "-----------------------------------------------------------------------------------------------------------"
	plot = d18OV


#Map code:

'''
#Plots regional box maps
fig = plt.figure()

ax = fig.add_subplot(241)
ax.set_title("SEA")
SEAmap = bm(projection = 'cea', llcrnrlat=box['latbottom_degN']['SEA'],urcrnrlat=box['lattop_degN']['SEA'], llcrnrlon=box['lonleft_degE']['SEA'],urcrnrlon=box['lonright_degE']['SEA'],resolution='c')
SEAmap.drawcoastlines()
SEAmap.drawparallels(np.arange(-90.,120.,30.))
SEAmap.drawmeridians(np.arange(0.,360.,60.))
SEAmap.drawmapboundary(fill_color='0.3')
cs = SEAmap.contourf(bmlon, bmlat, plotwhat, shading = 'flat', latlon=True)
SEAmap.colorbar(cs,location='bottom',pad="5%")


ax = fig.add_subplot(242)
ax.set_title("MC")
MCmap = bm(projection = 'cea', llcrnrlat=box['latbottom_degN']['MC'],urcrnrlat=box['lattop_degN']['MC'], llcrnrlon=box['lonleft_degE']['MC'],urcrnrlon=box['lonright_degE']['MC'],resolution='c')
MCmap.drawcoastlines()
MCmap.drawparallels(np.arange(-90.,120.,30.))
MCmap.drawmeridians(np.arange(0.,360.,60.))
MCmap.drawmapboundary(fill_color='0.3')
cs = MCmap.contourf(bmlon, bmlat, plotwhat, shading = 'flat', latlon=True)
MCmap.colorbar(cs,location='bottom',pad="5%")

ax = fig.add_subplot(243)
ax.set_title("IND")
INDmap = bm(projection = 'cea', llcrnrlat=box['latbottom_degN']['IND'],urcrnrlat=box['lattop_degN']['IND'], llcrnrlon=box['lonleft_degE']['IND'],urcrnrlon=box['lonright_degE']['IND'],resolution='c')
INDmap.drawcoastlines()
INDmap.drawparallels(np.arange(-90.,120.,30.))
INDmap.drawmeridians(np.arange(0.,360.,60.))
INDmap.drawmapboundary(fill_color='0.3')
cs = INDmap.contourf(bmlon, bmlat, plotwhat, shading = 'flat', latlon=True)
INDmap.colorbar(cs,location='bottom',pad="5%")

ax = fig.add_subplot(244)
ax.set_title("AFR")
AFRmap = bm(projection = 'cea', llcrnrlat=box['latbottom_degN']['AFR'],urcrnrlat=box['lattop_degN']['AFR'], llcrnrlon=box['lonleft_degE']['AFR'],urcrnrlon=box['lonright_degE']['AFR'],resolution='c')
AFRmap.drawcoastlines()
AFRmap.drawparallels(np.arange(-90.,120.,30.))
AFRmap.drawmeridians(np.arange(0.,360.,60.))
AFRmap.drawmapboundary(fill_color='0.3')
cs = AFRmap.contourf(bmlon, bmlat, plotwhat, shading = 'flat', latlon=True)
AFRmap.colorbar(cs,location='bottom',pad="5%")

ax = fig.add_subplot(245)
ax.set_title("WEP")
WEPmap = bm(projection = 'cea', llcrnrlat=box['latbottom_degN']['WEP'],urcrnrlat=box['lattop_degN']['WEP'], llcrnrlon=box['lonleft_degE']['WEP'],urcrnrlon=box['lonright_degE']['WEP'],resolution='c')
WEPmap.drawcoastlines()
WEPmap.drawparallels(np.arange(-90.,120.,30.))
WEPmap.drawmeridians(np.arange(0.,360.,60.))
WEPmap.drawmapboundary(fill_color='0.3')
cs = WEPmap.contourf(bmlon, bmlat, plotwhat, shading = 'flat', latlon=True)
WEPmap.colorbar(cs,location='bottom',pad="5%")

ax = fig.add_subplot(246)
ax.set_title("CEP")
CEPmap = bm(projection = 'cea', llcrnrlat=box['latbottom_degN']['CEP'],urcrnrlat=box['lattop_degN']['CEP'], llcrnrlon=box['lonleft_degE']['CEP'],urcrnrlon=box['lonright_degE']['CEP'],resolution='c')
CEPmap.drawcoastlines()
CEPmap.drawparallels(np.arange(-90.,120.,30.))
CEPmap.drawmeridians(np.arange(0.,360.,60.))
CEPmap.drawmapboundary(fill_color='0.3')
cs = CEPmap.contourf(bmlon, bmlat, plotwhat, shading = 'flat', latlon=True)
CEPmap.colorbar(cs,location='bottom',pad="5%")

ax = fig.add_subplot(247)
ax.set_title("EEP")
EEPmap = bm(projection = 'cea', llcrnrlat=box['latbottom_degN']['EEP'],urcrnrlat=box['lattop_degN']['EEP'], llcrnrlon=box['lonleft_degE']['EEP'],urcrnrlon=box['lonright_degE']['EEP'],resolution='c')
EEPmap.drawcoastlines()
EEPmap.drawparallels(np.arange(-90.,120.,30.))
EEPmap.drawmeridians(np.arange(0.,360.,60.))
EEPmap.drawmapboundary(fill_color='0.3')
cs = EEPmap.contourf(bmlon, bmlat, plotwhat, shading = 'flat', latlon=True)
EEPmap.colorbar(cs,location='bottom',pad="5%")

ax = fig.add_subplot(248)
ax.set_title("SA")
SAmap = bm(projection = 'cea', llcrnrlat=box['latbottom_degN']['SA'],urcrnrlat=box['lattop_degN']['SA'], llcrnrlon=box['lonleft_degE']['SA'],urcrnrlon=box['lonright_degE']['SA'],resolution='c')
SAmap.drawcoastlines()
SAmap.drawparallels(np.arange(-90.,120.,30.))
SAmap.drawmeridians(np.arange(0.,360.,60.))
SAmap.drawmapboundary(fill_color='0.3')
cs = SAmap.contourf(bmlon, bmlat, plotwhat, shading = 'flat', latlon=True)
SAmap.colorbar(cs,location='bottom',pad="5%")

plt.show()
'''

# for i,row in enumerate(plot):
# 	for j,e in enumerate(row):
# 		if e > 200:
# 			plot[i,j] = 200
# 		if e < -200:
# 			plot[i,j] = -200

#Plots tropics
m = bm(projection = 'cea', llcrnrlat=-35,urcrnrlat=35, llcrnrlon=-60,urcrnrlon=300,resolution='c')
m.drawcoastlines()
# m.drawparallels(np.arange(-90.,120.,15.))
# m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, plot, shading = 'flat', latlon = True, cmap=plt.cm.RdBu_r)
m.colorbar(cs, location='bottom', pad="5%")

plt.title(plotwhat + " at level " + level + "\n" + fn)
plt.show()
# plt.savefig("/Users/kyleniezgoda/Desktop/" + fn + ".ps", bbox_inches='tight', dpi = 500)







#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Old code:

'''
#Find the 850mb level index
lev850 = [n for n,i in enumerate([abs(l - 850) for l in levs]) if i == min([abs(l - 850) for l in levs])]
#Find the ground level index (last one = closest to ground)
levground = -1

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Find box index edges for each region

SEAbox = [[n for n,i in enumerate([abs(l - box['latbottom_degN']['SEA']) for l in lats]) if i == min([abs(l - box['latbottom_degN']['SEA']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lattop_degN']['SEA']) for l in lats]) if i == min([abs(l - box['lattop_degN']['SEA']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lonleft_degE']['SEA']) for l in lons]) if i == min([abs(l - box['lonleft_degE']['SEA']) for l in lons])], \
		[n for n,i in enumerate([abs(l - box['lonright_degE']['SEA']) for l in lons]) if i == min([abs(l - box['lonright_degE']['SEA']) for l in lons])]]

MCbox = [[n for n,i in enumerate([abs(l - box['latbottom_degN']['MC']) for l in lats]) if i == min([abs(l - box['latbottom_degN']['MC']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lattop_degN']['MC']) for l in lats]) if i == min([abs(l - box['lattop_degN']['MC']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lonleft_degE']['MC']) for l in lons]) if i == min([abs(l - box['lonleft_degE']['MC']) for l in lons])], \
		[n for n,i in enumerate([abs(l - box['lonright_degE']['MC']) for l in lons]) if i == min([abs(l - box['lonright_degE']['MC']) for l in lons])]]

INDbox = [[n for n,i in enumerate([abs(l - box['latbottom_degN']['IND']) for l in lats]) if i == min([abs(l - box['latbottom_degN']['IND']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lattop_degN']['IND']) for l in lats]) if i == min([abs(l - box['lattop_degN']['IND']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lonleft_degE']['IND']) for l in lons]) if i == min([abs(l - box['lonleft_degE']['IND']) for l in lons])], \
		[n for n,i in enumerate([abs(l - box['lonright_degE']['IND']) for l in lons]) if i == min([abs(l - box['lonright_degE']['IND']) for l in lons])]]

AFRbox = [[n for n,i in enumerate([abs(l - box['latbottom_degN']['AFR']) for l in lats]) if i == min([abs(l - box['latbottom_degN']['AFR']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lattop_degN']['AFR']) for l in lats]) if i == min([abs(l - box['lattop_degN']['AFR']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lonleft_degE']['AFR']) for l in lons]) if i == min([abs(l - box['lonleft_degE']['AFR']) for l in lons])], \
		[n for n,i in enumerate([abs(l - box['lonright_degE']['AFR']) for l in lons]) if i == min([abs(l - box['lonright_degE']['AFR']) for l in lons])]]

EEPbox = [[n for n,i in enumerate([abs(l - box['latbottom_degN']['EEP']) for l in lats]) if i == min([abs(l - box['latbottom_degN']['EEP']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lattop_degN']['EEP']) for l in lats]) if i == min([abs(l - box['lattop_degN']['EEP']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lonleft_degE']['EEP']) for l in lons]) if i == min([abs(l - box['lonleft_degE']['EEP']) for l in lons])], \
		[n for n,i in enumerate([abs(l - box['lonright_degE']['EEP']) for l in lons]) if i == min([abs(l - box['lonright_degE']['EEP']) for l in lons])]]

CEPbox = [[n for n,i in enumerate([abs(l - box['latbottom_degN']['CEP']) for l in lats]) if i == min([abs(l - box['latbottom_degN']['CEP']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lattop_degN']['CEP']) for l in lats]) if i == min([abs(l - box['lattop_degN']['CEP']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lonleft_degE']['CEP']) for l in lons]) if i == min([abs(l - box['lonleft_degE']['CEP']) for l in lons])], \
		[n for n,i in enumerate([abs(l - box['lonright_degE']['CEP']) for l in lons]) if i == min([abs(l - box['lonright_degE']['CEP']) for l in lons])]]

WEPbox = [[n for n,i in enumerate([abs(l - box['latbottom_degN']['WEP']) for l in lats]) if i == min([abs(l - box['latbottom_degN']['WEP']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lattop_degN']['WEP']) for l in lats]) if i == min([abs(l - box['lattop_degN']['WEP']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lonleft_degE']['WEP']) for l in lons]) if i == min([abs(l - box['lonleft_degE']['WEP']) for l in lons])], \
		[n for n,i in enumerate([abs(l - box['lonright_degE']['WEP']) for l in lons]) if i == min([abs(l - box['lonright_degE']['WEP']) for l in lons])]]

SAbox = [[n for n,i in enumerate([abs(l - box['latbottom_degN']['SA']) for l in lats]) if i == min([abs(l - box['latbottom_degN']['SA']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lattop_degN']['SA']) for l in lats]) if i == min([abs(l - box['lattop_degN']['SA']) for l in lats])], \
		[n for n,i in enumerate([abs(l - box['lonleft_degE']['SA']) for l in lons]) if i == min([abs(l - box['lonleft_degE']['SA']) for l in lons])], \
		[n for n,i in enumerate([abs(l - box['lonright_degE']['SA']) for l in lons]) if i == min([abs(l - box['lonright_degE']['SA']) for l in lons])]]

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#set level used for analysis
lev = levground

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Extract H2OV data (needed for computation of delta values)

#DJF
SEA_H2OV_DJF = DJF.variables['H2OV'][0,lev,SEAbox[0][0]:SEAbox[1][0],SEAbox[2][0]:SEAbox[3][0]]
MC_H2OV_DJF = DJF.variables['H2OV'][0,lev,MCbox[0][0]:MCbox[1][0],MCbox[2][0]:MCbox[3][0]]
IND_H2OV_DJF = DJF.variables['H2OV'][0,lev,INDbox[0][0]:INDbox[1][0],INDbox[2][0]:INDbox[3][0]]
AFR_H2OV_DJF = DJF.variables['H2OV'][0,lev,AFRbox[0][0]:AFRbox[1][0],AFRbox[2][0]:AFRbox[3][0]]
EEP_H2OV_DJF = DJF.variables['H2OV'][0,lev,EEPbox[0][0]:EEPbox[1][0],EEPbox[2][0]:EEPbox[3][0]]
CEP_H2OV_DJF = DJF.variables['H2OV'][0,lev,CEPbox[0][0]:CEPbox[1][0],CEPbox[2][0]:CEPbox[3][0]]
WEP_H2OV_DJF = DJF.variables['H2OV'][0,lev,WEPbox[0][0]:WEPbox[1][0],WEPbox[2][0]:WEPbox[3][0]]
SA_H2OV_DJF = DJF.variables['H2OV'][0,lev,SAbox[0][0]:SAbox[1][0],SAbox[2][0]:SAbox[3][0]]

#DJF
SEA_H2OV_JJA = JJA.variables['H2OV'][0,lev,SEAbox[0][0]:SEAbox[1][0],SEAbox[2][0]:SEAbox[3][0]]
MC_H2OV_JJA = JJA.variables['H2OV'][0,lev,MCbox[0][0]:MCbox[1][0],MCbox[2][0]:MCbox[3][0]]
IND_H2OV_JJA = JJA.variables['H2OV'][0,lev,INDbox[0][0]:INDbox[1][0],INDbox[2][0]:INDbox[3][0]]
AFR_H2OV_JJA = JJA.variables['H2OV'][0,lev,AFRbox[0][0]:AFRbox[1][0],AFRbox[2][0]:AFRbox[3][0]]
EEP_H2OV_JJA = JJA.variables['H2OV'][0,lev,EEPbox[0][0]:EEPbox[1][0],EEPbox[2][0]:EEPbox[3][0]]
CEP_H2OV_JJA = JJA.variables['H2OV'][0,lev,CEPbox[0][0]:CEPbox[1][0],CEPbox[2][0]:CEPbox[3][0]]
WEP_H2OV_JJA = JJA.variables['H2OV'][0,lev,WEPbox[0][0]:WEPbox[1][0],WEPbox[2][0]:WEPbox[3][0]]
SA_H2OV_JJA = JJA.variables['H2OV'][0,lev,SAbox[0][0]:SAbox[1][0],SAbox[2][0]:SAbox[3][0]]


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Extract isotope vapor data
#(time = 1, lev = 30, lat = 96, lon = 144)

#Compute delta values:
#dD = (HDOV/H2OV - 1) * 1000
#d18O = (H218OV/H2OV - 1) * 1000

#DJF
SEA_H218OV_DJF = DJF.variables['H218OV'][0,lev,SEAbox[0][0]:SEAbox[1][0],SEAbox[2][0]:SEAbox[3][0]]
SEA_d18OV_DJF = (SEA_H218OV_DJF/SEA_H2OV_DJF - 1) * 1000
SEA_HDOV_DJF = DJF.variables['HDOV'][0,lev,SEAbox[0][0]:SEAbox[1][0],SEAbox[2][0]:SEAbox[3][0]]
SEA_dDV_DJF = (SEA_HDOV_DJF/SEA_H2OV_DJF - 1) * 1000

MC_H218OV_DJF = DJF.variables['H218OV'][0,lev,MCbox[0][0]:MCbox[1][0],MCbox[2][0]:MCbox[3][0]]
MC_d18OV_DJF = (MC_H218OV_DJF/MC_H2OV_DJF - 1) * 1000
MC_HDOV_DJF = DJF.variables['HDOV'][0,lev,MCbox[0][0]:MCbox[1][0],MCbox[2][0]:MCbox[3][0]]
MC_dDV_DJF = (MC_HDOV_DJF/MC_H2OV_DJF - 1) * 1000

IND_H218OV_DJF = DJF.variables['H218OV'][0,lev,INDbox[0][0]:INDbox[1][0],INDbox[2][0]:INDbox[3][0]]
IND_d18OV_DJF = (IND_H218OV_DJF/IND_H2OV_DJF - 1) * 1000
IND_HDOV_DJF = DJF.variables['HDOV'][0,lev,INDbox[0][0]:INDbox[1][0],INDbox[2][0]:INDbox[3][0]]
IND_dDV_DJF = (IND_HDOV_DJF/IND_H2OV_DJF - 1) * 1000

AFR_H218OV_DJF = DJF.variables['H218OV'][0,lev,AFRbox[0][0]:AFRbox[1][0],AFRbox[2][0]:AFRbox[3][0]]
AFR_d18OV_DJF = (AFR_H218OV_DJF/AFR_H2OV_DJF - 1) * 1000
AFR_HDOV_DJF = DJF.variables['HDOV'][0,lev,AFRbox[0][0]:AFRbox[1][0],AFRbox[2][0]:AFRbox[3][0]]
AFR_dDV_DJF = (AFR_HDOV_DJF/AFR_H2OV_DJF - 1) * 1000

EEP_H218OV_DJF = DJF.variables['H218OV'][0,lev,EEPbox[0][0]:EEPbox[1][0],EEPbox[2][0]:EEPbox[3][0]]
EEP_d18OV_DJF = (EEP_H218OV_DJF/EEP_H2OV_DJF - 1) * 1000
EEP_HDOV_DJF = DJF.variables['HDOV'][0,lev,EEPbox[0][0]:EEPbox[1][0],EEPbox[2][0]:EEPbox[3][0]]
EEP_dDV_DJF = (EEP_HDOV_DJF/EEP_H2OV_DJF - 1) * 1000

CEP_H218OV_DJF = DJF.variables['H218OV'][0,lev,CEPbox[0][0]:CEPbox[1][0],CEPbox[2][0]:CEPbox[3][0]]
CEP_d18OV_DJF = (CEP_H218OV_DJF/CEP_H2OV_DJF - 1) * 1000
CEP_HDOV_DJF = DJF.variables['HDOV'][0,lev,CEPbox[0][0]:CEPbox[1][0],CEPbox[2][0]:CEPbox[3][0]]
CEP_dDV_DJF = (CEP_HDOV_DJF/CEP_H2OV_DJF - 1) * 1000

WEP_H218OV_DJF = DJF.variables['H218OV'][0,lev,WEPbox[0][0]:WEPbox[1][0],WEPbox[2][0]:WEPbox[3][0]]
WEP_d18OV_DJF = (WEP_H218OV_DJF/WEP_H2OV_DJF - 1) * 1000
WEP_HDOV_DJF = DJF.variables['HDOV'][0,lev,WEPbox[0][0]:WEPbox[1][0],WEPbox[2][0]:WEPbox[3][0]]
WEP_dDV_DJF = (WEP_HDOV_DJF/WEP_H2OV_DJF - 1) * 1000

SA_H218OV_DJF = DJF.variables['H218OV'][0,lev,SAbox[0][0]:SAbox[1][0],SAbox[2][0]:SAbox[3][0]]
SA_d18OV_DJF = (SA_H218OV_DJF/SA_H2OV_DJF - 1) * 1000
SA_HDOV_DJF = DJF.variables['HDOV'][0,lev,SAbox[0][0]:SAbox[1][0],SAbox[2][0]:SAbox[3][0]]
SA_dDV_DJF = (SA_HDOV_DJF/SA_H2OV_DJF - 1) * 1000

#JJA
SEA_H218OV_JJA = JJA.variables['H218OV'][0,lev,SEAbox[0][0]:SEAbox[1][0],SEAbox[2][0]:SEAbox[3][0]]
SEA_d18OV_JJA = (SEA_H218OV_JJA/SEA_H2OV_JJA - 1) * 1000
SEA_HDOV_JJA = JJA.variables['HDOV'][0,lev,SEAbox[0][0]:SEAbox[1][0],SEAbox[2][0]:SEAbox[3][0]]
SEA_dDV_JJA = (SEA_HDOV_JJA/SEA_H2OV_JJA - 1) * 1000

MC_H218OV_JJA = JJA.variables['H218OV'][0,lev,MCbox[0][0]:MCbox[1][0],MCbox[2][0]:MCbox[3][0]]
MC_d18OV_JJA = (MC_H218OV_JJA/MC_H2OV_JJA - 1) * 1000
MC_HDOV_JJA = JJA.variables['HDOV'][0,lev,MCbox[0][0]:MCbox[1][0],MCbox[2][0]:MCbox[3][0]]
MC_dDV_JJA = (MC_HDOV_JJA/MC_H2OV_JJA - 1) * 1000

IND_H218OV_JJA = JJA.variables['H218OV'][0,lev,INDbox[0][0]:INDbox[1][0],INDbox[2][0]:INDbox[3][0]]
IND_d18OV_JJA = (IND_H218OV_JJA/IND_H2OV_JJA - 1) * 1000
IND_HDOV_JJA = JJA.variables['HDOV'][0,lev,INDbox[0][0]:INDbox[1][0],INDbox[2][0]:INDbox[3][0]]
IND_dDV_JJA = (IND_HDOV_JJA/IND_H2OV_JJA - 1) * 1000

AFR_H218OV_JJA = JJA.variables['H218OV'][0,lev,AFRbox[0][0]:AFRbox[1][0],AFRbox[2][0]:AFRbox[3][0]]
AFR_d18OV_JJA = (AFR_H218OV_JJA/AFR_H2OV_JJA - 1) * 1000
AFR_HDOV_JJA = JJA.variables['HDOV'][0,lev,AFRbox[0][0]:AFRbox[1][0],AFRbox[2][0]:AFRbox[3][0]]
AFR_dDV_JJA = (AFR_HDOV_JJA/AFR_H2OV_JJA - 1) * 1000

EEP_H218OV_JJA = JJA.variables['H218OV'][0,lev,EEPbox[0][0]:EEPbox[1][0],EEPbox[2][0]:EEPbox[3][0]]
EEP_d18OV_JJA = (EEP_H218OV_JJA/EEP_H2OV_JJA - 1) * 1000
EEP_HDOV_JJA = JJA.variables['HDOV'][0,lev,EEPbox[0][0]:EEPbox[1][0],EEPbox[2][0]:EEPbox[3][0]]
EEP_dDV_JJA = (EEP_HDOV_JJA/EEP_H2OV_JJA - 1) * 1000

CEP_H218OV_JJA = JJA.variables['H218OV'][0,lev,CEPbox[0][0]:CEPbox[1][0],CEPbox[2][0]:CEPbox[3][0]]
CEP_d18OV_JJA = (CEP_H218OV_JJA/CEP_H2OV_JJA - 1) * 1000
CEP_HDOV_JJA = JJA.variables['HDOV'][0,lev,CEPbox[0][0]:CEPbox[1][0],CEPbox[2][0]:CEPbox[3][0]]
CEP_dDV_JJA = (CEP_HDOV_JJA/CEP_H2OV_JJA - 1) * 1000

WEP_H218OV_JJA = JJA.variables['H218OV'][0,lev,WEPbox[0][0]:WEPbox[1][0],WEPbox[2][0]:WEPbox[3][0]]
WEP_d18OV_JJA = (WEP_H218OV_JJA/WEP_H2OV_JJA - 1) * 1000
WEP_HDOV_JJA = JJA.variables['HDOV'][0,lev,WEPbox[0][0]:WEPbox[1][0],WEPbox[2][0]:WEPbox[3][0]]
WEP_dDV_JJA = (WEP_HDOV_JJA/WEP_H2OV_JJA - 1) * 1000

SA_H218OV_JJA = JJA.variables['H218OV'][0,lev,SAbox[0][0]:SAbox[1][0],SAbox[2][0]:SAbox[3][0]]
SA_d18OV_JJA = (SA_H218OV_JJA/SA_H2OV_JJA - 1) * 1000
SA_HDOV_JJA = JJA.variables['HDOV'][0,lev,SAbox[0][0]:SAbox[1][0],SAbox[2][0]:SAbox[3][0]]
SA_dDV_JJA = (SA_HDOV_JJA/SA_H2OV_JJA - 1) * 1000


min18O = np.amin(SEA_d18OV_DJF)
max18O = np.amax(SEA_d18OV_DJF)
minD = np.amin(SEA_dDV_DJF)
maxD = np.amax(SEA_dDV_DJF)
for i in [MC_d18OV_DJF,IND_d18OV_DJF,AFR_d18OV_DJF,EEP_d18OV_DJF,CEP_d18OV_DJF,WEP_d18OV_DJF,SA_d18OV_DJF, \
			SEA_d18OV_JJA,MC_d18OV_JJA,IND_d18OV_JJA,AFR_d18OV_JJA,EEP_d18OV_JJA,CEP_d18OV_JJA,WEP_d18OV_JJA,SA_d18OV_JJA]:
	if np.amin(i) < min18O:
		min18O = np.amin(i)
	if np.amax(i) > max18O:
		max18O = np.amax(i)
for i in [MC_dDV_DJF,IND_dDV_DJF,AFR_dDV_DJF,EEP_dDV_DJF,CEP_dDV_DJF,WEP_dDV_DJF,SA_dDV_DJF, \
			SEA_dDV_JJA,MC_dDV_JJA,IND_dDV_JJA,AFR_dDV_JJA,EEP_dDV_JJA,CEP_dDV_JJA,WEP_dDV_JJA,SA_dDV_JJA]:
	if np.amin(i) < minD:
		minD = np.amin(i)
	if np.amax(i) > maxD:
		maxD = np.amax(i)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Extract temp data

#DJF
SEA_T_DJF = DJF.variables['T'][0,lev,SEAbox[0][0]:SEAbox[1][0],SEAbox[2][0]:SEAbox[3][0]]
MC_T_DJF = DJF.variables['T'][0,lev,MCbox[0][0]:MCbox[1][0],MCbox[2][0]:MCbox[3][0]]
IND_T_DJF = DJF.variables['T'][0,lev,INDbox[0][0]:INDbox[1][0],INDbox[2][0]:INDbox[3][0]]
AFR_T_DJF = DJF.variables['T'][0,lev,AFRbox[0][0]:AFRbox[1][0],AFRbox[2][0]:AFRbox[3][0]]
EEP_T_DJF = DJF.variables['T'][0,lev,EEPbox[0][0]:EEPbox[1][0],EEPbox[2][0]:EEPbox[3][0]]
CEP_T_DJF = DJF.variables['T'][0,lev,CEPbox[0][0]:CEPbox[1][0],CEPbox[2][0]:CEPbox[3][0]]
WEP_T_DJF = DJF.variables['T'][0,lev,WEPbox[0][0]:WEPbox[1][0],WEPbox[2][0]:WEPbox[3][0]]
SA_T_DJF = DJF.variables['T'][0,lev,SAbox[0][0]:SAbox[1][0],SAbox[2][0]:SAbox[3][0]]

#JJA
SEA_T_JJA = JJA.variables['T'][0,lev,SEAbox[0][0]:SEAbox[1][0],SEAbox[2][0]:SEAbox[3][0]]
MC_T_JJA = JJA.variables['T'][0,lev,MCbox[0][0]:MCbox[1][0],MCbox[2][0]:MCbox[3][0]]
IND_T_JJA = JJA.variables['T'][0,lev,INDbox[0][0]:INDbox[1][0],INDbox[2][0]:INDbox[3][0]]
AFR_T_JJA = JJA.variables['T'][0,lev,AFRbox[0][0]:AFRbox[1][0],AFRbox[2][0]:AFRbox[3][0]]
EEP_T_JJA = JJA.variables['T'][0,lev,EEPbox[0][0]:EEPbox[1][0],EEPbox[2][0]:EEPbox[3][0]]
CEP_T_JJA = JJA.variables['T'][0,lev,CEPbox[0][0]:CEPbox[1][0],CEPbox[2][0]:CEPbox[3][0]]
WEP_T_JJA = JJA.variables['T'][0,lev,WEPbox[0][0]:WEPbox[1][0],WEPbox[2][0]:WEPbox[3][0]]
SA_T_JJA = JJA.variables['T'][0,lev,SAbox[0][0]:SAbox[1][0],SAbox[2][0]:SAbox[3][0]]

minT = np.amin(SEA_T_DJF)
maxT = np.amax(SEA_T_DJF)
for i in [MC_T_DJF,IND_T_DJF,AFR_T_DJF,EEP_T_DJF,CEP_T_DJF,WEP_T_DJF,SA_T_DJF,SEA_T_JJA,MC_T_DJF,IND_T_DJF,AFR_T_DJF,EEP_T_DJF,CEP_T_DJF,WEP_T_DJF,SA_T_DJF]:
	if np.amin(i) < minT:
		minT = np.amin(i)
	if np.amax(i) > maxT:
		maxT = np.amax(i)


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DJF.close()
JJA.close()
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Select what to plot ---> 1 = plot, 0 = don't plot
LMWL = 1
T_vs_d18OV = 0
T_vs_dDV = 0

plt.figure(1)
if LMWL:
	#Plot LMWLs
	plt.subplot(241)
	plt.plot(np.reshape(SEA_d18OV_DJF, (1,-1)), np.reshape(SEA_dDV_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(SEA_d18OV_JJA, (1,-1)), np.reshape(SEA_dDV_JJA, (1,-1)), 'ro')
	plt.axis([min18O, max18O, minD, maxD])
	plt.title("SEA")

	plt.subplot(242)
	plt.plot(np.reshape(MC_d18OV_DJF, (1,-1)), np.reshape(MC_dDV_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(MC_d18OV_JJA, (1,-1)), np.reshape(MC_dDV_JJA, (1,-1)), 'ro')
	plt.axis([min18O, max18O, minD, maxD])
	plt.title("MC")

	plt.subplot(243)
	plt.plot(np.reshape(IND_d18OV_DJF, (1,-1)), np.reshape(IND_dDV_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(IND_d18OV_JJA, (1,-1)), np.reshape(IND_dDV_JJA, (1,-1)), 'ro')
	plt.axis([min18O, max18O, minD, maxD])
	plt.title("IND")

	plt.subplot(244)
	plt.plot(np.reshape(AFR_d18OV_DJF, (1,-1)), np.reshape(AFR_dDV_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(AFR_d18OV_JJA, (1,-1)), np.reshape(AFR_dDV_JJA, (1,-1)), 'ro')
	plt.axis([min18O, max18O, minD, maxD])
	plt.title("AFR")

	plt.subplot(245)
	plt.plot(np.reshape(EEP_d18OV_DJF, (1,-1)), np.reshape(EEP_dDV_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(EEP_d18OV_JJA, (1,-1)), np.reshape(EEP_dDV_JJA, (1,-1)), 'ro')
	plt.axis([min18O, max18O, minD, maxD])
	plt.title("EEP")

	plt.subplot(246)
	plt.plot(np.reshape(CEP_d18OV_DJF, (1,-1)), np.reshape(CEP_dDV_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(CEP_d18OV_JJA, (1,-1)), np.reshape(CEP_dDV_JJA, (1,-1)), 'ro')
	plt.axis([min18O, max18O, minD, maxD])
	plt.title("CEP")

	plt.subplot(247)
	plt.plot(np.reshape(WEP_d18OV_DJF, (1,-1)), np.reshape(WEP_dDV_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(WEP_d18OV_JJA, (1,-1)), np.reshape(WEP_dDV_JJA, (1,-1)), 'ro')
	plt.axis([min18O, max18O, minD, maxD])
	plt.title("WEP")

	plt.subplot(248)
	plt.plot(np.reshape(SA_d18OV_DJF, (1,-1)), np.reshape(SA_dDV_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(SA_d18OV_JJA, (1,-1)), np.reshape(SA_dDV_JJA, (1,-1)), 'ro')
	plt.axis([min18O, max18O, minD, maxD])
	plt.title("SA")

if T_vs_d18OV:
	#Plot T vs d18OV
	plt.subplot(241)
	plt.plot(np.reshape(SEA_d18OV_DJF, (1,-1)), np.reshape(SEA_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(SEA_d18OV_JJA, (1,-1)), np.reshape(SEA_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("SEA")

	plt.subplot(242)
	plt.plot(np.reshape(MC_d18OV_DJF, (1,-1)), np.reshape(MC_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(MC_d18OV_JJA, (1,-1)), np.reshape(MC_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("MC")

	plt.subplot(243)
	plt.plot(np.reshape(IND_d18OV_DJF, (1,-1)), np.reshape(IND_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(IND_d18OV_JJA, (1,-1)), np.reshape(IND_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("IND")

	plt.subplot(244)
	plt.plot(np.reshape(AFR_d18OV_DJF, (1,-1)), np.reshape(AFR_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(AFR_d18OV_JJA, (1,-1)), np.reshape(AFR_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("AFR")

	plt.subplot(245)
	plt.plot(np.reshape(EEP_d18OV_DJF, (1,-1)), np.reshape(EEP_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(EEP_d18OV_JJA, (1,-1)), np.reshape(EEP_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("EEP")

	plt.subplot(246)
	plt.plot(np.reshape(CEP_d18OV_DJF, (1,-1)), np.reshape(CEP_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(CEP_d18OV_JJA, (1,-1)), np.reshape(CEP_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("CEP")

	plt.subplot(247)
	plt.plot(np.reshape(WEP_d18OV_DJF, (1,-1)), np.reshape(WEP_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(WEP_d18OV_JJA, (1,-1)), np.reshape(WEP_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("WEP")

	plt.subplot(248)
	plt.plot(np.reshape(SA_d18OV_DJF, (1,-1)), np.reshape(SA_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(SA_d18OV_JJA, (1,-1)), np.reshape(SA_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("SA")

if T_vs_dDV:
	#Plot T vs dD
	plt.subplot(241)
	plt.plot(np.reshape(SEA_dDV_DJF, (1,-1)), np.reshape(SEA_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(SEA_dDV_JJA, (1,-1)), np.reshape(SEA_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("SEA")

	plt.subplot(242)
	plt.plot(np.reshape(MC_dDV_DJF, (1,-1)), np.reshape(MC_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(MC_dDV_JJA, (1,-1)), np.reshape(MC_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("MC")

	plt.subplot(243)
	plt.plot(np.reshape(IND_dDV_DJF, (1,-1)), np.reshape(IND_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(IND_dDV_JJA, (1,-1)), np.reshape(IND_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("IND")

	plt.subplot(244)
	plt.plot(np.reshape(AFR_dDV_DJF, (1,-1)), np.reshape(AFR_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(AFR_dDV_JJA, (1,-1)), np.reshape(AFR_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("AFR")

	plt.subplot(245)
	plt.plot(np.reshape(EEP_dDV_DJF, (1,-1)), np.reshape(EEP_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(EEP_dDV_JJA, (1,-1)), np.reshape(EEP_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("EEP")

	plt.subplot(246)
	plt.plot(np.reshape(CEP_dDV_DJF, (1,-1)), np.reshape(CEP_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(CEP_dDV_JJA, (1,-1)), np.reshape(CEP_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("CEP")

	plt.subplot(247)
	plt.plot(np.reshape(WEP_dDV_DJF, (1,-1)), np.reshape(WEP_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(WEP_dDV_JJA, (1,-1)), np.reshape(WEP_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("WEP")

	plt.subplot(248)
	plt.plot(np.reshape(SA_dDV_DJF, (1,-1)), np.reshape(SA_T_DJF, (1,-1)), 'bo')
	plt.plot(np.reshape(SA_dDV_JJA, (1,-1)), np.reshape(SA_T_JJA, (1,-1)), 'ro')
	plt.axis([minD, maxD, minT, maxT])
	plt.title("SA")

plt.show()
plt.savefig("/Users/kyleniezgoda/Desktop/test.png")

'''