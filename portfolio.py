class P:
	def __init__(self):
		import os
		self.home = os.path.expanduser("~")
		import pandas as pd
		self.portfolio = pd.read_csv(self.home+"/Desktop/portfolio.csv", dtype = {'Stock':str,'PurchasedPrice':float,'Quantity':int})
		self.equity = pd.read_csv(self.home+"/Desktop/stocks.csv", dtype = {"date" : str})
		self.store()
		
	def store(self):
		# Store this as a dict
		self.dict = {}
		for n, stock in enumerate(self.portfolio.Stock):
			date = [d for d in self.equity['date']]
			currentPrice = [p for p in self.equity[stock]]
			gain = [c - self.portfolio.PurchasedPrice[n] for c in currentPrice]
			self.dict[stock] = {'date' : date, \
			'currentPrice' : currentPrice, \
			'gain' : gain
			}

	def plotPrice(self, stocks = None):
		import matplotlib.pyplot as plt
		from datetime import datetime as dt
		dateobj = [dt.strptime(d, '%m/%d/%y') for d in self.equity['date']]
		date = [dt.strftime(d, '%b %d') for d in dateobj]
		if stocks is None:
			stocks = [s for s in self.portfolio.Stock]
		
		for stock in stocks:
			plt.plot(dateobj, self.dict[stock]['currentPrice'], label = stock)
		
		plt.legend()
		plt.show()

	def plotGain(self, stocks = None):
		import matplotlib.pyplot as plt
		from datetime import datetime as dt
		dateobj = [dt.strptime(d, '%m/%d/%y') for d in self.equity['date']]
		date = [dt.strftime(d, '%b %d') for d in dateobj]
		if stocks is None:
			stocks = [s for s in self.portfolio.Stock]
		
		for stock in stocks:
			plt.plot(dateobj, self.dict[stock]['gain'], label = stock)
		
		plt.legend()
		plt.show()