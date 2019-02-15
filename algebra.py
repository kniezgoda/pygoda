# Conversions for cam variables. 
# Using dict language, create a new variable with the following:
# "var" : [multiplicationFactor, additionFactor, newUnits]

algebra = { \
			"PRECT_H2O" : [60*60*24*1000, 0, "kg/m2/day"], \
			"PRECT_H218O" : [60*60*24*1000, 0, "kg/m2/day"], \
			"PRECT_HDO" : [60*60*24*1000, 0, "kg/m2/day"], \
			"PRECT" : [60*60*24*1000, 0, "kg/m2/day"], \
			"PRECC" : [60*60*24*1000, 0, "kg/m2/day"], \
			"PRECL" : [60*60*24*1000, 0, "kg/m2/day"], \
			"Q" : [1000, 0, "g-H2O/kg"], \
			"QFLX" : [60*60*24, 0, "kg/m2/day"], \
			"QFLX_H2O" : [60*60*24, 0, "kg/m2/day"], \
			"QFLX_H218O" : [60*60*24, 0, "kg/m2/day"], \
			"QFLX_HDO" : [60*60*24, 0, "kg/m2/day"], \
			"T" : [1, -273, "C"]
		}