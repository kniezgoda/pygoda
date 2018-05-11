import numpy as np

f = open('/Users/kyleniezgoda/Desktop/env_temp_grad.csv', 'r')
data = np.genfromtxt(f, delimiter=',')
f.close()

env_pressure = np.delete(data[:,0], 0)
env_temp = np.delete(data[:,1], 0) - 273

f = open('/Users/kyleniezgoda/Desktop/parcel_profile.csv', 'r')
data = np.genfromtxt(f, delimiter=',')
f.close()

parcel_pressure = np.delete(data[:,0], 0)
parcel_temp = np.delete(data[:,1], 0) 

import matplotlib.pyplot as plt
# # Three-panel plot
# fig1, (ax1, ax2) = plt.subplots(1,2,sharey=True)
# # Temperature
# ax1.plot(env_temp,env_pressure,'o-')
# ax1.set_ylabel('Depth (m)')
# ax1.set_ylim(ax1.get_ylim()[::-1]) #this reverses the yaxis (i.e. deep at the bottom)
# ax1.set_xlabel('Environmental Temperature (C)')

# # Salinity
# ax2.plot(parcel_temp,parcel_pressure,'o-r')
# ax2.set_xlabel('Parcel Temp')
# ax2.yaxis.set_visible(False) # This erases the y ticks

# plt.show()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(env_temp,env_pressure,'o-')
 
# Draw x label
ax1.set_xlabel('Temperature (C)')
# Draw y label
ax1.set_ylabel('Depth (m)')
ax1.set_ylim(ax1.get_ylim()[::-1]) #this reverses the yaxis (i.e. deep at the bottom)
 
plt.show()

