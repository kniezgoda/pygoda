from pygoda import ncgoda
import os
import numpy as np
import matplotlib.pyplot as plt
os.chdir('/home/server/student/homes/kniezgod/model_runs/F.C5.2deg.wiso.defaultSSTICE_kn002')

c = ncgoda("F.C5.2deg.wiso.defaultSSTICE_kn002_JJA_climo.nc")
os.chdir('../F.C5.2deg.wiso.obs6kSST_kn003')
t = ncgoda("F.C5.2deg.wiso.obs6kSST_kn003_JJA_climo.nc")

##################
### 90E to 90W ###
##################
box = [-66.5,66.5,-90,90]

c.variable("OMEGA", box)
c.data = c.isobar(85000)
t.variable("OMEGA", box)
t.data = t.isobar(85000)
cZM = np.mean(c.data, axis = 1)
tZM = np.mean(t.data, axis = 1)

clat_interp = np.linspace(c.boxlat[0], c.boxlat[-1], num = 1000)
tlat_interp = np.linspace(t.boxlat[0], t.boxlat[-1], num = 1000)

cZM_interp = np.interp(clat_interp, c.boxlat, cZM)
tZM_interp = np.interp(tlat_interp, c.boxlat, tZM)

c_lat_of_max_cldhgh = clat_interp[list(cZM_interp*-1).index(max(list(cZM_interp*-1)))]
t_lat_of_max_cldhgh = tlat_interp[list(tZM_interp*-1).index(max(list(tZM_interp*-1)))]

print "\n\nRegion layout : [bottomLat, topLat, leftLon, rightLon]"
print "\n\nRegion : " + str(box)
print "Test max OMEGA lat : " + str(t_lat_of_max_cldhgh) # 8.45 deg
print "Control max OMEGA lat : " + str(c_lat_of_max_cldhgh) # 6.66 deg

print "Latitude difference : " + str(t_lat_of_max_cldhgh-c_lat_of_max_cldhgh)

plt.subplot(1,2,1)
plt.plot(clat_interp, cZM_interp, label = "PI", color = "r")
plt.plot(tlat_interp, tZM_interp, label = "MH", color = "b")
plt.ylabel("Vertical velocity at 500 mb (Pa/s)")
plt.xlabel("Latitude")
plt.legend()

##################
### 90E to 90W ###
##################
box = [-66.5, 66.5, 90, -90]

c.variable("OMEGA", box)
c.data = c.isobar(85000)
t.variable("OMEGA", box)
t.data = t.isobar(85000)
cZM = np.mean(c.data, axis = 1)
tZM = np.mean(t.data, axis = 1)

clat_interp = np.linspace(c.boxlat[0], c.boxlat[-1], num = 1000)
tlat_interp = np.linspace(t.boxlat[0], t.boxlat[-1], num = 1000)

cZM_interp = np.interp(clat_interp, c.boxlat, cZM)
tZM_interp = np.interp(tlat_interp, c.boxlat, tZM)

c_lat_of_max_cldhgh = clat_interp[list(cZM_interp*-1).index(max(list(cZM_interp*-1)))]
t_lat_of_max_cldhgh = tlat_interp[list(tZM_interp*-1).index(max(list(tZM_interp*-1)))]

print "\n\nRegion : " + str(box)
print "Test max OMEGA lat : " + str(t_lat_of_max_cldhgh) # 8.45 deg
print "Control max OMEGA lat : " + str(c_lat_of_max_cldhgh) # 6.66 deg

print "Latitude difference : " + str(t_lat_of_max_cldhgh-c_lat_of_max_cldhgh)

plt.subplot(1,2,2)
plt.plot(clat_interp, cZM_interp, label = "PI", color = "r")
plt.plot(tlat_interp, tZM_interp, label = "MH", color = "b")
plt.ylabel("Vertical velocity at 500 mb (Pa/s)")
plt.xlabel("Latitude")
plt.legend()
plt.show()


