'''
#################
### Objective ###
#################

Find a range of values that encompasses a pre-set fraction of the data and is centered on zero
'''
def diffColor(d, alpha = .95):
	import numpy as np
	numcells = d.size
	MAX = np.max(d)
	MIN = np.min(d)
	delta = MAX - MIN
	test = min(abs(MIN), abs(MAX))
	a = len(np.where(np.logical_and(d>-test, d<test))[0])
	frac = float(a)/float(numcells)
	trace = 0
	# print "try number: " + str(trace)
	# print "frac is: " + str(frac)
	while frac < alpha:
		test = test + delta/20
		a = len(np.where(np.logical_and(d>-test, d<test))[0])
		frac = float(a)/float(numcells)
		trace += 1
		# print "try number: " + str(trace)
		# print "frac is: " + str(frac)
	return test
