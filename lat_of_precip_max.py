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
box = [-90,90,-90,90]

c.variable("PRECT", box)
t.variable("PRECT", box)
cZM = np.mean(c.data, axis = 1)
tZM = np.mean(t.data, axis = 1)

clat_interp = np.linspace(c.boxlat[0], c.boxlat[-1], num = 1000)
tlat_interp = np.linspace(t.boxlat[0], t.boxlat[-1], num = 1000)

cZM_interp = np.interp(clat_interp, c.boxlat, cZM)
tZM_interp = np.interp(tlat_interp, c.boxlat, tZM)

c_lat_of_max_precip = clat_interp[list(cZM_interp).index(max(list(cZM_interp)))]
t_lat_of_max_precip = tlat_interp[list(tZM_interp).index(max(list(tZM_interp)))]

print "\n\nRegion layout : [bottomLat, topLat, leftLon, rightLon]"
print "\n\nRegion : " + str(box)
print "Test max precip lat : " + str(t_lat_of_max_precip) # 8.45 deg
print "Control max precip lat : " + str(c_lat_of_max_precip) # 6.66 deg

print "Latitude difference : " + str(t_lat_of_max_precip-c_lat_of_max_precip)

plt.subplot(1,2,1)
plt.plot(clat_interp, cZM_interp, label = "PI", color = "r")
plt.plot(tlat_interp, tZM_interp, label = "MH", color = "b")
plt.ylabel("Precip")
plt.xlabel("Latitude")
plt.legend()

##################
### 90E to 90W ###
##################
box = [-90, 90, 90, -90]

c.variable("PRECT", box)
t.variable("PRECT", box)
cZM = np.mean(c.data, axis = 1)
tZM = np.mean(t.data, axis = 1)

clat_interp = np.linspace(c.boxlat[0], c.boxlat[-1], num = 1000)
tlat_interp = np.linspace(t.boxlat[0], t.boxlat[-1], num = 1000)

cZM_interp = np.interp(clat_interp, c.boxlat, cZM)
tZM_interp = np.interp(tlat_interp, c.boxlat, tZM)

c_lat_of_max_precip = clat_interp[list(cZM_interp).index(max(list(cZM_interp)))]
t_lat_of_max_precip = tlat_interp[list(tZM_interp).index(max(list(tZM_interp)))]

print "\n\nRegion : " + str(box)
print "Test max precip lat : " + str(t_lat_of_max_precip) # 8.45 deg
print "Control max precip lat : " + str(c_lat_of_max_precip) # 6.66 deg

print "Latitude difference : " + str(t_lat_of_max_precip-c_lat_of_max_precip)

plt.subplot(1,2,2)
plt.plot(clat_interp, cZM_interp, label = "PI", color = "r")
plt.plot(tlat_interp, tZM_interp, label = "MH", color = "b")
plt.ylabel("Precip")
plt.xlabel("Latitude")
plt.legend()
plt.show()


