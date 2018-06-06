from pygoda import camgoda

# del(Rp) = (Re_bar*del_E + E_bar*del_Re + Rc_bar*del_C + C_bar*del_Rc - Rp_bar*del_P) / P_bar

cfile = 'fc5.2deg.wiso.piControl_kn028/AMWG.climoFiles/fc5.2deg.wiso.piControl_kn028_JAS_climo.nc'
tfile = 'fc5.2deg.wiso.mh6ka_kn032/AMWG.climoFiles/fc5.2deg.wiso.mh6ka_kn032_JAS_climo.nc'

cnc = camgoda(cfile)
tnc = camgoda(tfile)

<<<<<<< HEAD
box = [5, 30, 345, 45]
=======
box = (0, 35, 330, 35)
>>>>>>> d1e1c0fad5a708324398f76a9df51534676bf01a


# Done with "delta values" 
tPRECT_H218O = tnc.variable("PRECT_H218O", box)
tPRECT_H2O  = tnc.variable("PRECT_H2O", box)
tRp = (tPRECT_H218O / tPRECT_H2O - 1) * 1000

tQFLX_H218O = tnc.variable("QFLX_H218O", box)
tQFLX_H2O = tnc.variable("QFLX_H2O", box)
tRe = (tQFLX_H218O / tQFLX_H2O - 1) * 1000

cPRECT_H218O = cnc.variable("PRECT_H218O", box)
cPRECT_H2O  = cnc.variable("PRECT_H2O", box)
cRp = (cPRECT_H218O / cPRECT_H2O - 1) * 1000

cQFLX_H218O = cnc.variable("QFLX_H218O", box)
cQFLX_H2O = cnc.variable("QFLX_H2O", box)
cRe = (cQFLX_H218O / cQFLX_H2O - 1) * 1000


# 1) Re_bar
# Re_bar = (tRe + cRe) / 2
Re_bar = (tRe + cRe) / 2

# 2) del_E
del_E = tQFLX_H2O - cQFLX_H2O

# 3) E_bar
E_bar = (tQFLX_H2O + cQFLX_H2O) / 2

# 4) del_Re
del_Re = tRe - cRe

# 8) del_Rc
tRc = (tPRECT_H218O - tQFLX_H218O) / (tPRECT_H2O - tQFLX_H2O)
cRc = (cPRECT_H218O - cQFLX_H218O) / (cPRECT_H2O - cQFLX_H2O)
del_Rc = (tRc - cRc) * 1000

# 5) Rc_bar
Rc_bar = (((tRc - 1) * 1000) + ((cRc - 1) * 1000)) / 2 

# 6) del_C
tC = tPRECT_H2O - tQFLX_H2O
cC = cPRECT_H2O - cQFLX_H2O
del_C = tC - cC

# 7) C_bar
C_bar = (tC + cC) / 2

# 10) del_P
del_P = tPRECT_H2O - cPRECT_H2O

# 11) P_bar
P_bar = (tPRECT_H2O + cPRECT_H2O) / 2

# 9) Rp_bar
Rp_bar = (tRp + cRp) / 2 


##################
# Combine it all #
##################

one = Re_bar*del_E
two = E_bar*del_Re 
three = Rc_bar*del_C
four = C_bar*del_Rc 
five = -Rp_bar*del_P
SUM = one + two + three + four + five

print([x/SUM for x in [one, two, three, four, five]])

del_Rp = SUM / P_bar




<<<<<<< HEAD
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap as bm 

bmlon, bmlat = np.meshgrid(cnc.boxlon, cnc.boxlat)
southern_lat, northern_lat = cnc.boxlat[[0,-1]] 
left_lon, right_lon = cnc.boxlon[[0,-1]]
if 0 in cnc.boxlon[1:-2]: # if we cross the gml
	left_lon = cnc.boxlon[0]-360

fig = plt.figure()
clev = np.linspace(-25,25,11)

plt.subplot(211)
plt.title('del_Rp')
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, del_Rp, np.linspace(-4,4,17), shading = 'flat', latlon = True, cmap = plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")

plt.subplot(256)
plt.title('Re_bar*del_E')
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, one, clev, shading = 'flat', latlon = True, cmap = plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")

plt.subplot(257)
plt.title('E_bar*del_Re')
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, two, clev, shading = 'flat', latlon = True, cmap = plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")

plt.subplot(258)
plt.title('Rc_bar*del_C')
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, three, clev, shading = 'flat', latlon = True, cmap = plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")

plt.subplot(259)
plt.title('C_bar*del_Rc')
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, four, clev, shading = 'flat', latlon = True, cmap = plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")

plt.subplot(2,5,10)
plt.title('-Rp_bar*del_P')
m = bm(projection = 'cea', llcrnrlat=southern_lat,urcrnrlat=northern_lat, llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='c')
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
cs = m.contourf(bmlon, bmlat, five, clev, shading = 'flat', latlon = True, cmap = plt.cm.RdBu_r)
cbar = m.colorbar(cs, location='right', pad="5%")

plt.show()
=======
>>>>>>> d1e1c0fad5a708324398f76a9df51534676bf01a
#-------#



# Done with ratios 
tPRECT_H218O = tnc.variable("PRECT_H218O", box)
tPRECT_H2O  = tnc.variable("PRECT_H2O", box)
tRp = tPRECT_H218O / tPRECT_H2O

tQFLX_H218O = tnc.variable("QFLX_H218O", box)
tQFLX_H2O = tnc.variable("QFLX_H2O", box)
tRe = tQFLX_H218O / tQFLX_H2O

cPRECT_H218O = cnc.variable("PRECT_H218O", box)
cPRECT_H2O  = cnc.variable("PRECT_H2O", box)
cRp = cPRECT_H218O / cPRECT_H2O

cQFLX_H218O = cnc.variable("QFLX_H218O", box)
cQFLX_H2O = cnc.variable("QFLX_H2O", box)
cRe = cQFLX_H218O / cQFLX_H2O


# 1) Re_bar
# Re_bar = (tRe + cRe) / 2
Re_bar = (tRe + cRe) / 2

# 2) del_E
del_E = tQFLX_H2O - cQFLX_H2O

# 3) E_bar
E_bar = (tQFLX_H2O + cQFLX_H2O) / 2

# 4) del_Re
del_Re = tRe - cRe

# 8) del_Rc
tRc = (tPRECT_H218O - tQFLX_H218O) / (tPRECT_H2O - tQFLX_H2O)
cRc = (cPRECT_H218O - cQFLX_H218O) / (cPRECT_H2O - cQFLX_H2O)
del_Rc = tRc - cRc

# 5) Rc_bar
Rc_bar = (tRc + cRc) / 2 

# 6) del_C
tC = tPRECT_H2O - tQFLX_H2O
cC = cPRECT_H2O - cQFLX_H2O
del_C = tC - cC

# 7) C_bar
C_bar = (tC + cC) / 2

# 10) del_P
del_P = tPRECT_H2O - cPRECT_H2O

# 11) P_bar
P_bar = (tPRECT_H2O + cPRECT_H2O) / 2

# 9) Rp_bar
Rp_bar = (tRp + cRp) / 2 


##################
# Combine it all #
##################

one = Re_bar*del_E
two = E_bar*del_Re 
three = Rc_bar*del_C
four = C_bar*del_Rc 
five = -Rp_bar*del_P
SUM = one + two + three + four + five

print([x/SUM for x in [one, two, three, four, five]])

del_Rp = SUM / P_bar