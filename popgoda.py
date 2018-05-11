class popgoda:
	def __init__(self, dataset):
		from netCDF4 import Dataset as ds
		from pygoda import ncdump
		self.dataset = dataset
		if type(dataset) is str:
			self.dataset = ds(dataset, "r")
		self.dims = [d for d in self.dataset.dimensions]
		self.vars = [v for v in self.dataset.variables]
		self.dimlen = [len(self.dataset.dimensions[d]) for d in self.dims]
		self.ulat = self.dataset.variables['ULAT'][:]
		self.ulon = self.dataset.variables['ULONG'][:]
		
		self.var = ''
		self.units = ''
		self.long_name = ''

		#self.levs = self.dataset.variables['lev'][:]
		self.isTime = any(d == "time" for d in self.dims)
		if self.isTime:
			self.ntime = self.dimlen[self.dims.index("time")]
			if self.ntime == 1:
				self.isTime = False

	def variable(self, var, box = None):
		import numpy as np
		DATA = self.dataset.variables[var][:]
		if np.ma.is_masked(DATA):
			DATA = np.ma.filled(DATA, fill_value = np.nan)
			DATA = DATA.squeeze()

		self.ndims = len(DATA.shape)

		self.mask = False
		if box is not None:
			latmask = (self.ulat>box[0]) & (self.ulat<box[1])
			lonmask = (self.ulon>box[2]) & (self.ulon<box[3])
			mask = latmask & lonmask
			self.mask = np.invert(mask)

		if self.ndims == 3:
			DATA = DATA[0,:,:]

		DATA = np.ma.array(DATA, mask = self.mask)
		DATA = np.ma.filled(DATA, fill_value = np.nan)


		self.data = DATA


