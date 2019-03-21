# pygoda
K. Niezgoda python scripts, mostly for handling CAM netcdf files.

The idea here is for pygoda.py to contain methods and classes for manipulating data from isotope-enabled CESM netcdfs. Some methods will be broadly useful for other applications (e.g., sigmaFilter), hence why this file is called "pygoda" and not "ncgoda". Two important classes in this file are "camgoda" and "popgoda" which provide a class for creating smart instance of CAM and POP history files. 

[Aside: I will definitely want to merge popgoda and camgoda in the future, and just have one class called ncgoda that is able to handle pretty much any netcdf files. The problem there is that pop and cam have different "things" in them that will take some time to figure out how to do cleanly at the same time.]

The idea here is that things like precipitation stable isotope ratios, which have to be calculated from the history files, are all done manually with functions from the class (in this case, it would be camgoda.PRECT_d18O() or camgoda.ExtractData("PRECT_d18O")). Other functions in pygoda but outside popgoda and camgoda are various methods that I've built over time to do useful things, like find indices corresponding to lat and long ranges. Some methods are no longer used or are redunancies from legacy versions of the methods. This file needs some cleaning...

Outside of pygoda, the remaining .py files are for making graphs and maps. MOST files should use parser arguments, so that the code should be run from the command line with argument flags (e.g., -v VARIABLE, -lats LAT1 LAT2, etc..) The naming convention is pretty bad right now, but the idea is based on these ideas:
1) Code that creates figures of test, control, and test-control anomalies should be called Difference***.py. For instance, code that creates a test map, control map, and difference map of a variable field should be called DifferenceMap.py. Similarly, code that creates pressure versus latitude figures for control, test, and difference should be called DifferencePressureVsLat.py.
2) Code that creates one figure from one file should be named similar to above, but without "Difference", i.e. Map.py, PressureVsLat.py, etc..

The code should always be able to handle multiple input files and multiple variables. Standard flags that are included in most (if not all) code parser args are:

1) -v VARIABLE {no default}                    : variable name that can be used in camgoda.ExtractData()

2) -lats BOTTOM_LAT TOP_LAT {-90 90}           : latitude bounds for the "box" argument in camgoda and popgoda

3) -lons LEFT_LON RIGHT_LON {0 360}            : similar to -lats. LONS ARE ALWAYS > 0

4) -nosave, -show, -dev {False False, False}   : Provides functionality for saving and/or displaying plots to the screen. These flags do not take arguments, they are either set or not. -nosave will not write out figures to a file, -show will display the figure on the screen, and -dev will set -nosave and -show (so nothing is saved and figures are printed to the screen.)

5) -cdir CDIR -tdir TDIR {./ ./}                        : the control and test directory, where files will be searched for.

6) -grep GREP                               : a string to search for in cdir and tdir. This flag is used as the first argument in findClimoFile, and is surrounded by the "*" wildcard.

7) -t TFILE -c CFILE {None None}                           : a specified test and control file. If set, these will override any -cdir/-tdir and grep functionality. 

