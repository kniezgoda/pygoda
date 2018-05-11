from pygoda import ncgoda
from pygoda import findClimoFile
import numpy as np
import matplotlib.pyplot as plt

# Control file
cf = findClimoFile("F.C5.2deg.wiso.defaultSSTICE_kn002/*JJA*")[0]
# Test file
tf = findClimoFile("F.C5.2deg.wiso.obs6kSST_kn003/*JJA*")[0]
# Heating file
hf = findClimoFile("F.C5.2deg.wiso.defaultSSTICE.HeatingExperiments.branch.withHeating_kn012/*JJA*")[0]

box = [10,30,0,95]

# Open the files
c = ncgoda(cf)
t = ncgoda(tf)
h = ncgoda(hf)

# Calculate DTCOND vertical profiles
c_dtcond = c.variable("DTCOND", box)
t_dtcond = t.variable("DTCOND", box)
h_dtcond = h.variable("DTCOND", box)
c_dtcond_vert = np.mean(c_dtcond, (1,2))
t_dtcond_vert = np.mean(t_dtcond, (1,2))
h_dtcond_vert = np.mean(h_dtcond, (1,2))
#diff_dtcond_vert = t_dtcond_vert - c_dtcond_vert

# Calculate pressures
c.CalculatePressure(box = box)
c_pm = np.mean(c.P_m, (1,2))
t.CalculatePressure(box = box)
t_pm = np.mean(t.P_m, (1,2))
h.CalculatePressure(box = box)
h_pm = np.mean(h.P_m, (1,2))

# Plot
#fig = plt.figure()
plt.plot(c_dtcond_vert,c_pm/100, color = "black", label = "PI")
plt.plot(t_dtcond_vert, c_pm/100, color = "blue", label = "MH")
plt.plot(h_dtcond_vert, c_pm/100, color = "red", label = "Heating")
plt.legend(loc = 'best')
plt.xlabel('DTCOND (K/s)')
plt.ylabel('pressure (hPa)')
plt.ylim(1000,0)

plt.show()

