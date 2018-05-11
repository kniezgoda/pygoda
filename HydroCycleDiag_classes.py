#####################################################
### Raw netcdf variables, no 'calculate' function ###
#####################################################

class PRECT_H2O:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['PRECT_H2O'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.test = testdata.variables['PRECT_H2O'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.diff = self.test - self.control
		self.units = "mm/day"
		self.long_name = cntrldata.variables['PRECT_H2O'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : np.linspace(0,16,17), "ANN_NA" : np.linspace(0,9,19), "ANN_MC" : np.linspace(0,9,19), "ANN_IM" : np.linspace(0,9,19), "ANN_EP" : np.linspace(0,9,19), \
					  "ANNdiff_" : np.linspace(-2,2,21), "ANNdiff_NA" : np.linspace(-1.5,1.5,21), "ANNdiff_MC" : np.linspace(-1.5,1.5,21), "ANNdiff_IM" : np.linspace(0,9,19), "ANNdiff_EP" : np.linspace(0,9,19), \
					  "DJF_" : np.linspace(0,16,17), "DJF_NA" : np.linspace(0,16,17), "DJF_MC" : np.linspace(0,16,17), "DJF_IM" : np.linspace(0,9,19), "DJF_EP" : np.linspace(0,9,19), \
					  "DJFdiff_" : np.linspace(-5,5,21), "DJFdiff_NA" : np.linspace(-5,5,21), "DJFdiff_MC" : np.linspace(-5,5,21), "DJFdiff_IM" : np.linspace(0,9,19), "DJFdiff_EP" : np.linspace(0,9,19), \
					  "JJA_" : np.linspace(0,30,21), "JJA_NA" : np.linspace(0,30,21), "JJA_MC" : np.linspace(0,30,21), "JJA_IM" : np.linspace(0,9,19), "JJA_EP" : np.linspace(0,9,19), \
					  "JJAdiff_" : np.linspace(-6,6,21), "JJAdiff_NA" : np.linspace(-6,6,21), "JJAdiff_MC" : np.linspace(-6,6,21), "JJAdiff_IM" : np.linspace(0,9,19), "JJAdiff_EP" : np.linspace(0,9,19), \
					  }.get(season,19)

class PRECT_H218O:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['PRECT_H218O'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.test = testdata.variables['PRECT_H218O'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.diff = self.test - self.control
		self.units = "mm/day"
		self.long_name = cntrldata.variables['PRECT_H218O'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)


class PRECT_HDO:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['PRECT_HDO'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.test = testdata.variables['PRECT_HDO'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.diff = self.test - self.control
		self.units = "mm/day"
		self.long_name = cntrldata.variables['PRECT_HDO'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

		
class PRECC:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['PRECC'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.test = testdata.variables['PRECC'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.diff = self.test - self.control
		self.units = "mm/day"
		self.long_name = testdata.variables['PRECC'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class PRECL:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['PRECL'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.test = testdata.variables['PRECL'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.diff = self.test - self.control
		self.units = "mm/day"
		self.long_name = cntrldata.variables['PRECL'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class PRECSC:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['PRECSC'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.test = testdata.variables['PRECSC'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.diff = self.test - self.control
		self.units = "mm/day"
		self.long_name = cntrldata.variables['PRECSC'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class PRECSL:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['PRECSL'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.test = testdata.variables['PRECSL'][:].squeeze() * 1000 * 60 * 60 * 24 # Convert to mm/day
		self.diff = self.test - self.control
		self.units = "mm/day"
		self.long_name = cntrldata.variables['PRECSL'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class PRECST:
	def __init__(self, cntrldata, testdata):
		precsc = PRECSC(cntrldata, testdata)
		precsl = PRECSL(cntrldata, testdata)
		self.control = precsc.control + precsl.control
		self.test = precsc.test + precsl.test
		self.diff = self.test - self.control
		self.units = "mm/day"
		self.long_name = "Total (convective + large-scale) snow rate"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)


class QFLX:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['QFLX'][:].squeeze() * 1000
		self.test = testdata.variables['QFLX'][:].squeeze() * 1000
		self.diff = self.test - self.control
		self.units = "g-H2O/m2/s"
		self.long_name = cntrldata.variables['QFLX'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class QFLX_H218O:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['QFLX_H218O'][:].squeeze()
		self.test = testdata.variables['QFLX_H218O'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['QFLX_H218O'].units
		self.long_name = cntrldata.variables['QFLX_H218O'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class QFLX_HDO:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['QFLX_HDO'][:].squeeze()
		self.test = testdata.variables['QFLX_HDO'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['QFLX_HDO'].units
		self.long_name = cntrldata.variables['QFLX_HDO'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class CLDHGH:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['CLDHGH'][:].squeeze()
		self.test = testdata.variables['CLDHGH'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['CLDHGH'].units
		self.long_name = testdata.variables['CLDHGH'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class PS:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['PS'][:].squeeze()
		self.test = testdata.variables['PS'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['PS'].units
		self.long_name = cntrldata.variables['PS'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class PSL:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['PSL'][:].squeeze()
		self.test = testdata.variables['PSL'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['PSL'].units
		self.long_name = cntrldata.variables['PSL'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
class H2OVsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['H2OV'][0,-1,:,:].squeeze()
		self.test = testdata.variables['H2OV'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['H2OV'].units
		self.long_name = cntrldata.variables['H2OV'].long_name + " at surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class H2OV850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "H2OV", 85000)
		self.test = extract_isobar(testdata, "H2OV", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['H2OV'].units
		self.long_name = cntrldata.variables['H2OV'].long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class H2OV500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "H2OV", 50000)
		self.test = extract_isobar(testdata, "H2OV", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['H2OV'].units
		self.long_name = cntrldata.variables['H2OV'].long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class H2OV200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "H2OV", 20000)
		self.test = extract_isobar(testdata, "H2OV", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['H2OV'].units
		self.long_name = cntrldata.variables['H2OV'].long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class H218OVsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['H218OV'][0,-1,:,:].squeeze()
		self.test = testdata.variables['H218OV'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['H218OV'].units
		self.long_name = cntrldata.variables['H218OV'].long_name + " at surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class H218OV850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "H218OV", 85000)
		self.test = extract_isobar(testdata, "H218OV", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['H218OV'].units
		self.long_name = cntrldata.variables['H218OV'].long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class H218OV500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "H218OV", 50000)
		self.test = extract_isobar(testdata, "H218OV", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['H218OV'].units
		self.long_name = cntrldata.variables['H218OV'].long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class H218OV200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "H218OV", 20000)
		self.test = extract_isobar(testdata, "H218OV", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['H218OV'].units
		self.long_name = cntrldata.variables['H218OV'].long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class HDOVsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['HDOV'][0,-1,:,:].squeeze()
		self.test = testdata.variables['HDOV'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['HDOV'].units
		self.long_name = cntrldata.variables['HDOV'].long_name + " at surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class HDOV850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "HDOV", 85000)
		self.test = extract_isobar(testdata, "HDOV", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['HDOV'].units
		self.long_name = cntrldata.variables['HDOV'].long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class HDOV500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "HDOV", 50000)
		self.test = extract_isobar(testdata, "HDOV", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['HDOV'].units
		self.long_name = cntrldata.variables['HDOV'].long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class HDOV200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "HDOV", 20000)
		self.test = extract_isobar(testdata, "HDOV", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['HDOV'].units
		self.long_name = cntrldata.variables['HDOV'].long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class Qsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['Q'][0,-1,:,:].squeeze() * 100
		self.test = testdata.variables['Q'][0,-1,:,:].squeeze() * 100
		self.diff = self.test - self.control
		self.units = "g-H2O/kg"
		self.long_name = cntrldata.variables['Q'].long_name + " at the surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class Q850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "Q", 85000) * 100
		self.test = extract_isobar(testdata, "Q", 85000) * 100
		self.diff = self.test - self.control
		self.units = "g-H2O/kg"
		self.long_name = cntrldata.variables['Q'].long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class Q500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "Q", 50000) * 100
		self.test = extract_isobar(testdata, "Q", 50000) * 100
		self.diff = self.test - self.control
		self.units = "g-H2O/kg"
		self.long_name = cntrldata.variables['Q'].long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class Q200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "Q", 20000) * 100
		self.test = extract_isobar(testdata, "Q", 20000) * 100
		self.diff = self.test - self.control
		self.units = "g-H2O/kg"
		self.long_name = cntrldata.variables['Q'].long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class Vsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['V'][0.-1,:,:].squeeze()
		self.test = testdata.variables['V'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['V'].units
		self.long_name = cntrldata.variables['V'].long_name + " at surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class V850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "V", 85000)
		self.test = extract_isobar(testdata, "V", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['V'].units
		self.long_name = "Meridional wind speed at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class V500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "V", 50000)
		self.test = extract_isobar(testdata, "V", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['V'].units
		self.long_name = "Meridional wind speed at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class V200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "V", 20000)
		self.test = extract_isobar(testdata, "V", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['V'].units
		self.long_name = "Meridional wind speed at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		
class VTsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['VT'][0,-1,:,:].squeeze()
		self.test = testdata.variables['VT'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['VT'].units
		self.long_name = cntrldata.variables['VT'].long_name + " at surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class VT850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "VT", 85000)
		self.test = extract_isobar(testdata, "VT", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['VT'].units
		self.long_name = VT.long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class VT500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "VT", 50000)
		self.test = extract_isobar(testdata, "VT", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['VT'].units
		self.long_name = VT.long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class VT200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "VT", 20000)
		self.test = extract_isobar(testdata, "VT", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['VT'].units
		self.long_name = VT.long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class VQsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['VQ'][0,-1,:,:].squeeze()
		self.test = testdata.variables['VQ'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['VQ'].units
		self.long_name = cntrldata.variables['VQ'].long_name + " at surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class VQ850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "VQ", 85000)
		self.test = extract_isobar(testdata, "VQ", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['VQ'].units
		self.long_name = VQ.long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class VQ500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "VQ", 50000)
		self.test = extract_isobar(testdata, "VQ", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['VQ'].units
		self.long_name = VQ.long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class VQ200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "VQ", 20000)
		self.test = extract_isobar(testdata, "VQ", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['VQ'].units
		self.long_name = VQ.long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class Usurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['U'][0,-1,:,:].squeeze()
		self.test = testdata.variables['U'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['U'].units
		self.long_name = "Zonal wind speed at surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class U850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "U", 85000)
		self.test = extract_isobar(testdata, "U", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['U'].units
		self.long_name = "Zonal wind speed at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class U500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "U", 50000)
		self.test = extract_isobar(testdata, "U", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['U'].units
		self.long_name = "Zonal wind speed at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class U200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "U", 20000)
		self.test = extract_isobar(testdata, "U", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['U'].units
		self.long_name = "Zonal wind speed at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class UTsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['UT'][0,-1,:,:].squeeze()
		self.test = testdata.variables['UT'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['UT'].units
		self.long_name = "Zonal temperature flux at the surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class UT850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "UT", 85000)
		self.test = extract_isobar(testdata, "UT", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['UT'].units
		self.long_name = UT.long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class UT500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "UT", 50000)
		self.test = extract_isobar(testdata, "UT", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['UT'].units
		self.long_name = UT.long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class UT200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "UT", 20000)
		self.test = extract_isobar(testdata, "UT", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['UT'].units
		self.long_name = UT.long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class UQsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['UQ'][0,-1,:,:].squeeze()
		self.test = testdata.variables['UQ'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['UQ'].units
		self.long_name = cntrldata.variables['UQ'].long_name + " at surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class UQ850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "UQ", 85000)
		self.test = extract_isobar(testdata, "UQ", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['UQ'].units
		self.long_name = UQ.long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class UQ500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "UQ", 50000)
		self.test = extract_isobar(testdata, "UQ", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['UQ'].units
		self.long_name = UQ.long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		

class UQ200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "UQ", 20000)
		self.test = extract_isobar(testdata, "UQ", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['UQ'].units
		self.long_name = UQ.long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class Tsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['T'][0,-1,:,:].squeeze()
		self.test = testdata.variables['T'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['T'].units
		self.long_name = cntrldata.variables['T'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class T850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "T", 85000)
		self.test = extract_isobar(testdata, "T", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['T'].units
		self.long_name = "Temperature at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class T500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "T", 50000)
		self.test = extract_isobar(testdata, "T", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['T'].units
		self.long_name = "Temperature at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		
class T200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "T", 20000)
		self.test = extract_isobar(testdata, "T", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['T'].units
		self.long_name = "Temperature at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class OMEGAsurface:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['OMEGA'][0,-1,:,:].squeeze()
		self.test = testdata.variables['OMEGA'][0,-1,:,:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['OMEGA'].units
		self.long_name = cntrldata.variables['OMEGA'].long_name + " at the surface"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class OMEGA850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "OMEGA", 85000)
		self.test = extract_isobar(testdata, "OMEGA", 85000)
		self.diff = self.test - self.control
		self.units = testdata.variables['OMEGA'].units
		self.long_name = OMEGA.long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		
class OMEGA500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "OMEGA", 50000)
		self.test = extract_isobar(testdata, "OMEGA", 50000)
		self.diff = self.test - self.control
		self.units = testdata.variables['OMEGA'].units
		self.long_name = OMEGA.long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
	
class OMEGA200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "OMEGA", 20000)
		self.test = extract_isobar(testdata, "OMEGA", 20000)
		self.diff = self.test - self.control
		self.units = testdata.variables['OMEGA'].units
		self.long_name = OMEGA.long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class Z3850:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "Z3", 85000) / 1000 # Convert to km
		self.test = extract_isobar(testdata, "Z3", 85000) / 1000 # Convert to km
		self.diff = self.test - self.control
		self.units = "km"
		self.long_name = cntrldata.variables['Z3'].long_name + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		
class Z3500:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "Z3", 50000) / 1000 # Convert to km
		self.test = extract_isobar(testdata, "Z3", 50000) / 1000 # Convert to km
		self.diff = self.test - self.control
		self.units = "km"
		self.long_name = cntrldata.variables['Z3'].long_name + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
		
class Z3200:
	def __init__(self, cntrldata, testdata):
		from pygoda import extract_isobar
		self.control = extract_isobar(cntrldata, "Z3", 20000) / 1000 # Convert to km
		self.test = extract_isobar(testdata, "Z3", 20000) / 1000 # Convert to km
		self.diff = self.test - self.control
		self.units = "km"
		self.long_name = cntrldata.variables['Z3'].long_name + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

'''
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class P:
	def __init__(self, cntrldata, testdata):
		from pygoda import PressureCalc
		self.control = PressureCalc(A.control, B.control, PS.control)
		self.test = PressureCalc(A.test, B.test, PS.test)
		self.diff = self.test - self.control
		self.units = "Pa"
		self.long_name = "Pressure"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class E:
	def __init__(self, cntrldata, testdata):
		self.control = Q.control/1000 * P.control / 0.622 # Vapor pressure, q*P/epsilon
		self.test = Q.test/1000 * P.test / 0.622
		self.diff = self.test - self.control
		self.units = "Pa"
		self.long_name = "Vapor pressure"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class E850:
	def __init__(self, cntrldata, testdata):
		self.control = Q850.control/1000 * 85000 / 0.622 # Vapor pressure, q*P/epsilon
		self.test = Q850.test/1000 * 85000 / 0.622 # Vapor pressure, q*P/epsilon
		self.diff = self.test - self.control
		self.units = "Pa"
		self.long_name = "Vapor pressure" + " at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class E500:
	def __init__(self, cntrldata, testdata):
		self.control = Q500.control/1000 * 50000 / 0.622 # Vapor pressure, q*P/epsilon
		self.test = Q500.test/1000 * 50000 / 0.622 # Vapor pressure, q*P/epsilon
		self.diff = self.test - self.control
		self.units = "Pa"
		self.long_name = "Vapor pressure" + " at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class E200:
	def __init__(self, cntrldata, testdata):
		self.control = Q200.control/1000 * 20000 / 0.622 # Vapor pressure, q*P/epsilon
		self.test = Q200.test/1000 * 20000 / 0.622 # Vapor pressure, q*P/epsilon
		self.diff = self.test - self.control
		self.units = "Pa"
		self.long_name = "Vapor pressure" + " at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class ES:
	def __init__(self, cntrldata, testdata):
		t = T(cntrldata, testdata)
		self.control = 2.53e8 * np.exp(-5.42e3 / t.control) # saturation vapor pressure, Ae^(-B/T)
		self.test = 2.53e8 * np.exp(-5.42e3 / t.test)
		self.diff = self.test - self.control
		self.units = "Pa"
		self.long_name = "Saturation vapor pressure"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class ES850:
	def __init__(self, cntrldata, testdata):
		t850 = T850(cntrldata, testdata)
		self.control = 2.53e8 * np.exp(-5.42e3 / t850.control) # saturation vapor pressure, Ae^(-B/T)
		self.test = 2.53e8 * np.exp(-5.42e3 / t850.test)
		self.diff = self.test - self.control
		self.units = "Pa"
		self.long_name = "Saturation vapor pressure"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class ES500:
	def __init__(self, cntrldata, testdata):
		t500 = T500(cntrldata, testdata)
		self.control = 2.53e8 * np.exp(-5.42e3 / t500.control) # saturation vapor pressure, Ae^(-B/T)
		self.test = 2.53e8 * np.exp(-5.42e3 / t500.test)
		self.diff = self.test - self.control
		self.units = "Pa"
		self.long_name = "Saturation vapor pressure"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class ES200:
	def __init__(self, cntrldata, testdata):
		t200 = T200(cntrldata, testdata)
		self.control = 2.53e8 * np.exp(-5.42e3 / t200.control) # saturation vapor pressure, Ae^(-B/T)
		self.test = 2.53e8 * np.exp(-5.42e3 / t200.test)
		self.diff = self.test - self.control
		self.units = "Pa"
		self.long_name = "Saturation vapor pressure"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class RH:
	def __init__(self, cntrldata, testdata):
		e = E(cntrldata, testdata)
		es = ES(cntrldata, testdata)
		self.control = e.control / es.control
		self.test = e.test / es.test
		self.diff = self.test - self.control
		self.units = "percent"
		self.long_name = "Relative humidity"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class RH850:
	def __init__(self, cntrldata, testdata):
		e850 = E850(cntrldata, testdata)
		es850 = ES850(cntrldata, testdata)
		self.control = e850.control / es850.control
		self.test = e850.test / es850.test
		self.diff = self.test - self.control
		self.units = "percent"
		self.long_name = "Relative humidity at 850 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class RH500:
	def __init__(self, cntrldata, testdata):
		e500 = E500(cntrldata, testdata)
		es500 = ES500(cntrldata, testdata)
		self.control = e500.control / es500.control
		self.test = e500.test / es500.test
		self.diff = self.test - self.control
		self.units = "percent"
		self.long_name = "Relative humidity at 500 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class RH200:
	def __init__(self, cntrldata, testdata):
		e200 = E200(cntrldata, testdata)
		es200 = ES200(cntrldata, testdata)
		self.control = e200.control / es200.control
		self.test = e200.test / es200.test
		self.diff = self.test - self.control
		self.units = "percent"
		self.long_name = "Relative humidity at 200 mb"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, \
					19, \.get19"JJAdiff199t19es19ueeze() * 60 * 6019rt to K/yr
		self.test = testdata.variables['TT19e(19365
		self.diff = self.test - self.co19 "K/ye19e = cn19TEND_TO1919clevs = {"ANN_" : 19, "ANN_NA" : 19, "A" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" 19,: 19, "ANNdiff_MC" : 19, "ANNdiff_IM"1919), "ANNdiff_EP"1919),\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "19EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class DCQ:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['DCQ'][:].squeeze()
		self.test = testdata.variables['DCQ'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['DCQ'].units
		self.long_name = cntrldata.variables['DCQ'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class DCQ:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['DCQ'][:].squeeze()
		self.test = testdata.variables['DCQ'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['DCQ'].units
		self.long_name = cntrldata.variables['DCQ'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class CMFDT:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['CMFDT'][:].squeeze()
		self.test = testdata.variables['CMFDT'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['CMFDT'].units
		self.long_name = cntrldata.variables['CMFDT'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class CMFDQR:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['CMFDQR'][:].squeeze()
		self.test = testdata.variables['CMFDQR'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['CMFDQR'].units
		self.long_name = cntrldata.variables['CMFDQR'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class DTCOND:
	def __init__(self, cntrldata, testdata):
		self.control = cntrldata.variables['DTCOND'][:].squeeze()
		self.test = testdata.variables['DTCOND'][:].squeeze()
		self.diff = self.test - self.control
		self.units = testdata.variables['DTCOND'].units
		self.long_name = cntrldata.variables['DTCOND'].long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)

class dRHDO_dt:
	def __init__(self, cntrldata, testdata):
		from pygoda import dR_dt
		self.control = dR_dt(cntrldatafname, "HDO")
		self.test = dR_dt(testdatafname, "HDO")
		self.diff = self.test - self.control
		self.units = "R / s"
		self.long_name = "Surface HDOV tendency due to horizontal advection"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
'''


######################################################################
### Variables that must be calculated using other netcdf variables ###
######################################################################

class PRECT_d18O:
	def __init__(self):
		self.units = "delta 18O (permil)"
		self.long_name = "delta 18O for PRECT"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : np.linspace(-24,0,21), "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : np.linspace(-5,5,21), "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : np.linspace(-24,0,21), "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : np.linspace(-5,5,21), "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
	def calculate(self, h218o, h2o):
		self.control = (h218o.control / h2o.control - 1) * 1000
		self.test = (h218o.test / h2o.test - 1) * 1000
		self.diff = self.test - self.control

class PRECT_dD:
	def __init__(self):
		self.units = "delta D (permil)"
		self.long_name = "delta D for PRECT"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : np.linspace(-150,0,21), "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : np.linspace(-25,25,21), "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : np.linspace(-150,0,21), "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : np.linspace(-30,30,21), "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
	def calculate(self, hdo, h2o):
		self.control = (hdo.control / h2o.control - 1) * 1000
		self.test = (hdo.test / h2o.test - 1) * 1000
		self.diff = self.test - self.control

class PRECT_dxs:
	def __init__(self):
		self.units = "d-excess (permil)"
		self.long_name = "d-excess (delta_D - 8*delta_18O) for PRECT"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : np.linspace(0,30,21), "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : np.linspace(-5,5,21), "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : np.linspace(0,35,21), "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : np.linspace(-5,5,21), "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
	def calculate(self, d18o, dd):
		self.control = dd.control - 8*d18o.control
		self.test = dd.test - 8*d18o.test
		self.diff = self.test - self.control

class dDV:
	def __init__(self):
		self.units = "permil"
		self.long_name = "delta D for VAPOR"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
	def calculate(self, hdov, h2ov):
		self.control = (hdov.control/h2ov.control - 1) * 1000
		self.test = (hdov.test/h2ov.test - 1) * 1000
		self.diff = self.test - self.control
		

class d18OV:
	def __init__(self):
		self.units = "permil"
		self.long_name = "delta 18O for VAPOR"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
	def calculate(self, h218ov, h2ov):
		self.control = (h218ov.control/h2ov.control - 1) * 1000
		self.test = (h218ov.test/h2ov.test - 1) * 1000
		self.diff = self.test - self.control
		

class dxsV:
	def __init__(self):
		self.units = "d-excess (permil)"
		self.long_name = "d-excess (delta_D - 8*delta_18O) for VAPOR"
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
	def calculate(self, d18ov, ddv):
		self.control = ddv.control - 8*d18ov.control
		self.test = ddv.test - 8*d18ov.test
		self.diff = self.test - self.control

class VQ_dD:
	def __init__(self, cntrldata, testdata):
		self.units = vq_hdo.units
		self.long_name = vq_hdo.long_name
		self.clevs = 19
	def set_clevs(self, season):
		import numpy as np
		self.clevs = {"ANN_" : 19, "ANN_NA" : 19, "ANN_MC" : 19, "ANN_IM" : 19, "ANN_EP" : 19,\
					  "ANNdiff_" : 19, "ANNdiff_NA" : 19, "ANNdiff_MC" : 19, "ANNdiff_IM" : 19, "ANNdiff_EP" : 19,\
					  "DJF_" : 19, "DJF_NA" : 19, "DJF_MC" : 19, "DJF_IM" : 19, "DJF_EP" : 19,\
					  "DJFdiff_" : 19, "DJFdiff_NA" : 19, "DJFdiff_MC" : 19, "DJFdiff_IM" : 19, "DJFdiff_EP" : 19,\
					  "JJA_" : 19, "JJA_NA" : 19, "JJA_MC" : 19, "JJA_IM" : 19, "JJA_EP" : 19,\
					  "JJAdiff_" : 19, "JJAdiff_NA" : 19, "JJAdiff_MC" : 19, "JJAdiff_IM" : 19, "JJAdiff_EP" : 19,\
					  }.get(season,19)
	def calculate(self, vq_h2o, vq_hdo):
		self.control = (vq_hdo.control/vq_h2o.control - 1) * 1000
		self.test = (vq_hdo.test/vq_h2o.test - 1) * 1000
		self.diff = self.test - self.control
