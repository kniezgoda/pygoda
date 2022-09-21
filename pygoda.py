# A collection of functions useful for manipulating and working with netcdf files in python
# Other statistical and analysis tools have been added over time
# kniezgo@gmail.com

def runningMean(x, N, mode = 'same'):
    '''
    Computes the running mean of x over N steps
    '''
    import numpy as np
    return np.convolve(x, np.ones((N,))/N, mode=mode)

# =========================================================================================== #

def bandPass(data, axis = 0, low_bound = 0, high_bound = 1):
    '''
    function used for band-passing data. 
    Arguments:
    1) data; np.array of the data you wish to band pass
    2) axis (default 0); the axis over which to apply the FFT and IFFT
    3) low_bound (default 0); the lowest frequency to keep
    4) high_bound (default 1); the highest frequency to keep
    
    Returns:
    1) the band-passed data array in the same shape as the input array
    
    : Notes :
    For the low_ and high_bound, the frequencies are for each time step.
    E.g., if data is a monthly timeseries over N years, then to low-pass filter for
    signals of less than 1 cycle/year, let high_bound to 1./12.
    To band-pass for signals between 1 and 3 cycles/year set low_bound to 1./12. and high_bound to 3./12.
    To band-pass for ENSO signals, you should set low_bound to approximately 1./(12.*7.) and 
    high_bound to approximately 1./(12.*4.) (1 cycle every 4 - 7 years)
    
    Leaving axis, low_, and high_bound as defaults will return the same data
    To high-pass only, set low_bound to the cutoff freq, and leave high_bound as default
    To low-pass, set high_bound to cutoff freq, leave low_bound as default
    To band-pass, set both.
    '''
    import numpy as np
    ndims_data= len(data.shape)
    fft_freqs = np.abs(np.fft.fftfreq(data.shape[axis]))
    mask = (fft_freqs >= low_bound) * (fft_freqs <= high_bound)
    # make mask the same shape as data by adding new axes
    for i in range(ndims_data - 1):
        mask = mask[...,np.newaxis]
    return np.fft.ifft(np.fft.fft(data, axis = 0)*mask, axis = 0)

def boxOut(data, box, lat_axis = -2, lon_axis = -1, grid = "2deg", returnGrid = False):
    ''' 
    The idea here is to replace the existing technique for extracting data from a box.
    At first, to save time on computing, I would extract the data only from the box I wanted, and then do the analysis.
    However, it has become necessary to often change boxes, and this requires a complete reload of all the data, which is time consuming.
    This function allows for the user to read in all the spatial data (i.e. box = None), and then choose the box after the data has been read in
    which can be faster if the analysis requires a constant testing of different boxes. This function is NOT intended to be used on 
    data that has been spatially sliced (i.e. box != None) because it assumes the spatial extent is -90:90 and 0:360 nominally.
    This function works on a data array of any shape, as long as two of the axes represent lats and lons.
    Args: 
    1) data; a numpy array with axes corresponding to lat and lon. The lat and lon range should be the entire cam grid
    2) box; a list, np.array, or tuple of the form (bottom_lat, top_lat, left_lon, right_lon) with bounds
    (-90:90, -90:90, 0:360, 0:360)
    OPTIONAL
    3/4) lat_axis/lon_axis; the axes for latitude and longitude, defaulted to -2 (lat) -1 (lon)
    5) grid ["2deg" or "1deg"]; this defines the latitudes and longitudes. Default behavior is to set to the default CAM nominal 2-deg FV grid.
    Returns:
    The boxed-out data array 
    '''
    import numpy as np
    from pygoda import find_indices
    if grid == "2deg":
        boxlat = np.linspace(-90,90,96)
        boxlon = np.linspace(0,357.5,144)
    elif grid == "1deg":
        boxlat = np.linspace(-90,90,192)
        boxlon = np.linspace(0,358.75,288)
    else:
        print("Grid argument must be either '2deg' or '1deg' --- returning None")
        return None
    idxs = find_indices(box, boxlat, boxlon)
    RETURN = np.take(np.take(data, idxs[0], axis = lat_axis), idxs[1], axis = lon_axis)
    if returnGrid:
        RETURN = [RETURN, np.take(boxlat, idxs[0]), np.take(boxlon, idxs[1])]
    
    # This is the data
    RETURN = np.take(np.take(data, idxs[0], axis = lat_axis), idxs[1], axis = lon_axis)
    
    # This is included if the user wants the boxlat and boxlon data as well
    if returnGrid:
        RETURN = [RETURN, np.take(boxlat, idxs[0]), np.take(boxlon, idxs[1])]
    return RETURN


def aggregate(d, group, axis = 0, fun = "mean", rm_nan = True):
    '''
    *
    *** This function aggregates d along a specified axis based on a corresponding grouping array.
    *
    *** len(d) along the specified axis must equal len(group)
    i.e., each element of group should correspond to the index-aligned element of d along the axis.
    *
    *** fun can either be "mean" or "sum"
    *
    *** d is first indexed along the specified axis by each unique element of group, 
    the function (mean or sum) is applied, 
    and then the resultants are arranged into an array of shape (len(group), ...) 
    where the remaining axes are preserved from the original shape of d.
    *
    *** If axis != 0, then the shape will be (..., len(group), ...) where len(group) is 
    placed at the axis position specified in the arguments
    '''
    import numpy as np
    if len(group) != np.ma.size(d, axis = axis):
        print("group and d are not the same length, exiting...")
        return
    group = np.array(group)
    unique_groups = np.unique(group)
    group_idx = [np.where(group == g) for g in unique_groups]
    d_grouped = [np.take(d, idxs, axis = axis).squeeze() for idxs in group_idx]
    if fun == "mean":
        if rm_nan:
            ret = np.array([np.nanmean(dgrouped, axis = axis) for dgrouped in d_grouped])
        else:
            ret = np.array([np.mean(dgrouped, axis = axis) for dgrouped in d_grouped])
    elif fun == "sum":
        if rm_nan:
            ret = np.array([np.nansum(dgrouped, axis = axis) for dgrouped in d_grouped])
        else:
            ret = np.array([np.sum(dgrouped, axis = axis) for dgrouped in d_grouped])
    elif fun == "std":
        if rm_nan:
            ret = np.array([np.nanstd(dgrouped, axis = axis) for dgrouped in d_grouped])
        else:
            ret = np.array([np.std(dgrouped, axis = axis) for dgrouped in d_grouped])
            
    else: 
        print("fun argument must be 'mean', 'sum', or 'std', exiting...")
        return
    return ret, unique_groups


####################
### find_indices ### 
####################

def find_indices(box, lats, lons):
    '''
    Author : Kyle Niezgoda
    
    This function is useful for finding the indices of lat and lon lists
    that correspond to the latitude and longitudes in a box.
    
    For instance, maybe you want to extract data for India only, or within
    the lat range of 5N to 35N and lon range of 65E to 90E. To do this, you need
    to know the indices in the lat and lon list that are closest to the India 
    corrdinates. 
    
    To find these indices, you would do the following:
    
    lats = f.variables['lat'][:]
    lons = f.variables['lons'][:]
    region = [5, 35, 65, 90])
    bottom_index, top_index, left_index, right_index = find_indices(region, lats, lons)
    
    Arguments
    ---------
        box : list, [bottom_edge, top_edge, left_edge, right_edge] all in degrees
    lats : list of latitudes from the netcdf file
    lons : list of longitudes from the netcdf file
    
    Returns
    -------
    list, [bottom_index, top_index, left_index, right_index]
    '''
    import numpy as np

    if box is None:
        return [np.arange(len(lats)), np.arange(len(lons))]

    box = list(box)

    if box[2] < 0:
        box[2] = 360 + box[2]

    if box[3] < 0:
        box[3] = 360 + box[3]   

    def cat_idx():
        if l > r:
            return np.concatenate((np.arange(l, len(lons)), np.arange(0,r)))
        elif l == r:
            return [l]
        else:
            return np.arange(l,r)

    b = abs(lats - box[0]).tolist().index(min(abs(lats-box[0])))
    t = abs(lats - box[1]).tolist().index(min(abs(lats-box[1])))
    l = abs(lons-(box[2])).tolist().index(min(abs(lons-(box[2]))))
    r = abs(lons-(box[3])).tolist().index(min(abs(lons-(box[3]))))

    if b == t:
        return [[b], cat_idx()]

    return [np.arange(b,t), cat_idx()]



# =========================================================================================== #


def sigmaFilter(a, sigma = 3, n_passes=1):
    # Filters data based on the distance from the mean in units of SD
    # Main argument, a, must be a 1-d list or np array
    # Optional args: 
    # sigma (default = 3); the number of SD to screen for
    # n_passes (default = 1); the number of times to apply the filter
    # (average and sd are re-calculated after each pass)
    # Returns: the filtered array with np.nan replaced for filtered values
    import numpy as np
    avg = np.nanmean(a)
    stddev = np.nanstd(a)
    for i in range(n_passes):
        a = [np.nan if (x < avg - sigma*stddev) | (x > avg + sigma*stddev) else x for x in a]
        avg = np.nanmean(a)
        stddev = np.nanstd(a)
    return np.array(a)

def medianFilter(a, n = 3, n_passes=1):
    import scipy.stats as stats
    import numpy as np
    master = []
    # What is the lowest index you can start finding the median for
    lower_bound = (n-1)/2
    # What about the last index you can find the median for?
    upper_bound = len(a) - (n+1)/2
    for i in range(len(a)):
        if i < lower_bound:
            master.append(a[i])
        elif i > upper_bound:
            master.append(a[i])
        else:
            master.append(np.median(a[i-lower_bound:i+1+lower_bound]))
    return master

def corr(a, b, lag=0):
    '''
    Kyle Niezgoda, spring 2018

    Compute the correlation between two 1-dimensional arrays (python lists)
    Index lags can be introduced with the lag argument (default of zero = no lag)
    Correlations are normalized by N_tau, the number of non-nan samples for that particular lag.
    Returns the correlation value, R_xy.
    Auto-correlation are fine, just set a and b to the same list.

    Need to fix something wrong when doing an auto-correlation with time-lag = 0. 
    The R value returned (Rxx, in this case) is not equal to 1, which it should be.
    Instead, it is equal to something close to 0.997, not too far off but not entirely correct either.
    from scipy.stats import describe as desc
    
    ---> Negative lag implies that a1 leads a2!
    '''
    import numpy as np
    if lag == 0:
        # Have to do this because b[:-lag] doesn't work for lag == 0.
        a1 = np.array(a)
        a2 = np.array(b)
    elif lag > 0:
        a1 = np.array(a[lag:])
        a2 = np.array(b[:-lag])
    else:
        a1 = np.array(a[:lag])
        a2 = np.array(b[-lag:])
    a1_mu = np.nanmean(a1)
    a1_stddev = np.nanstd(a1)
    a2_mu = np.nanmean(a2)
    a2_stddev = np.nanstd(a2)
    # Find the nans 
    N = sum(~np.isnan(a1*a2))
    if N == 0:
        return np.nan
    return np.nansum((a1 - a1_mu)  * (a2- a2_mu)) / N / (a1_stddev * a2_stddev)


def corr_2d (a,b,axis=0,minN=0):
    # This function is used for computing timeseries correlation between two 3d arrays
    # The 0'th dimension should be time, and the 1st and 2nd should be space (i.e. lat, lon)
    # The timeseries for each grid point in a (i.e. a[:,i,j]) is correlated against the same gridpoint in b (i.e. b[:,i,j])
    # This is faster than looping over every grid cell because the means are computed over the entire time axis
    # However, I don't have a way to compute lagged corelations yet.
    import numpy as np
    import xarray as xr
    a = np.array(a)
    b = np.array(b)
    a_mu = np.nanmean(a, axis = axis)
    a_sd = np.nanstd(a, axis = axis)
    b_mu = np.nanmean(b, axis = axis)
    b_sd = np.nanstd(b, axis = axis)
    try:
        N = np.sum(~np.isnan(a*b), axis = axis)
        N = np.array(xr.DataArray(N).where(xr.DataArray(N) > minN))
    except ValueError:
        print("From pygoda.corr_2d: input arrays are not same dimension!")
        print("Fix inputs or add length-1 dimensions using np.expand_dims(array, -1)!")
        return
    return (np.nansum((a-a_mu) * (b-b_mu), axis = axis) / N / (a_sd*b_sd))


def Nstar_auto(a):
    N = np.sum(~np.isnan(a))
    hold = []
    for k in range(-(N-1), N-1):
        if not abs(k) < 0.8 * N:
            continue
        else:
            r = corr(a, a, lag=k)
            hold.append(((float(N)-float(abs(k)))/float(N)) * r**2)
    return N / np.nansum(hold)


def eof(d, removeMeans = True, verbose = False):
    '''
    Largely copied and modified from Jim Lerzcak
    d in the input data array, 2-d, where rows are times and columns are locations
    
    d.shape should be [N,M] with N times and M locations

    returns [a, F, l]
    a: the amplitude time series
    F: the EOFs
    l: the eigenvalues

    ANNOYINGLY, THIS DOES NOT WORK. USE THE EOF PACKAGE IN PYTHON INSTEAD
    '''
    import numpy as np
    import os
    N, M = d.shape
    if removeMeans:
        d = d - np.mean(d, axis = 0)

    # Calculate the mean square matrix (or the covariance matrix, if removeMeans)
    D = np.zeros(shape = (M,M))
    for ir in range(0,M):
        if verbose:
            percent = round(float(ir)/D.size/2 * 100)
            if percent % 10 == 0:
                print(str(percent) + "% done...")
        D[ir,ir] = np.nanmean(d[:,ir] * d[:,ir])
        # Doing it this way cuts the number of operations in half
        for ic in range(ir+1,M):
            D[ir, ic] = np.nanmean(d[:,ir] * d[:,ic])
            D[ic, ir] = D[ir,ic]

    home = os.path.expanduser("~") + "/"
    np.savetxt(home+".tempD.csv", D, delimiter = ',')
    # Compute the singular values of D
    # l is the variance (diagonals already removed)
    # F is the eigen function array
    W, l, F = np.linalg.svd(D)

    # Calculate the amplitude functions
    try:
        a = np.matmul(d, F)
    except(AttributeError):
        a = np.dot(d,F)

    return(a, F, l)


# =========================================================================================== #


####################
### PressureCalc ###
####################

def PressureCalc(A, B, PS):
    '''
    Author : Kyle Niezgoda
    Date : June 9, 2017

    This function creates a 3-d array of pressure values.
    For CESM, the 3-d pressure variable created and printed 
    out in the history files is not the true pressure. You 
    have to apply the coefficients A and B to the values, 
    along with the surface pressure, to get the real pressure
    since CESM has terrain following pycnals (or something
    like that.) The equation is P = A + B * PS w/ PS = surface
    pressure.

    Arguments
    ---------
    A : list of length n where n is the number of vertical levels (hyam)
    B : same as A (hybm)
    PS : 2-d array of surface pressure values

    Returns
    -------
    array : 3-d array of pressure values with dimensions =
    [n, PS.shape[0], PS.shape[1]] (lev, lat, lon)
    '''
    import numpy as np
    pressure = []
    for n in range(len(B)):
        pressure.append(np.array(A[n] * 100000 + B[n] * PS))
    return np.array(pressure)


# =========================================================================================== #


def camdates(start_year = None, end_year = None, months = [1,2,3,4,5,6,7,8,9,10,11,12], days = False, hours = None):
    import datetime

    # Compute the end day corresponding to the end month
    end_month = months[-1]
    if (end_month == 1) | (end_month == 3) | (end_month == 5) | (end_month == 7) | (end_month == 8) | (end_month == 10) | (end_month == 12):
        end_day = 31
    elif (end_month == 4) | (end_month == 6) | (end_month == 9) | (end_month == 11):
        end_day = 30
    else:
        end_day = 28

    # Handle for no years 
    # In this case, only months or month-days will be output
    strip_years = False
    if start_year is None or end_year is None:
        start_year = 0
        end_year = 0
        strip_years = True
        
    # Create the starting and ending date
    dates = []
    start_date = datetime.date(year = 1950 + start_year, month = months[0], day = 1)
    end_date = datetime.date(year = 1950 + end_year, month = end_month, day = end_day)

    d = start_date
    delta = datetime.timedelta(days = 1)
    #if hours is not None:
    #   delta = datetime.timedelta(hours = hours)

    # Loop until the end date is met
    while d <= end_date:
        # Datetime doesn't do dates before 1950 or something like that...
        # This is how I fix that
        YEAR = str(d.year - 1950)
        # Add leading zeros to the front of the string to match to cam history file naming structure
        # Years have 4 digits (YYYY), months and days have 2 digits (MM, DD)
        for i in range(4-len(YEAR)):
            YEAR = '0' + YEAR

        MONTH = str(d.month)
        if int(MONTH) in months:
            for i in range(2-len(MONTH)):
                MONTH = '0' + MONTH
        else:
            # Skip if we're not using this month
            d += delta
            continue

        if days:
            DAY = str(d.day)
            # Skip leap days because they don't exist in our cam runs
            if (int(MONTH) == 2) & (int(DAY) == 29):
                d += delta
                continue
            for i in range(2-len(DAY)):
                DAY = '0' + DAY
            if hours is not None:
                for iii in range(int(24/hours)):
                    HOUR = str(hours*iii*60*60)
                    for i in range(5-len(HOUR)):
                        HOUR = '0' + HOUR
                    dates.append(YEAR + '-' + MONTH + '-' + DAY + '-' + HOUR)
                d += delta
            else:
                dates.append(YEAR + '-' + MONTH + '-' + DAY)
                d += delta
        else:
            if d.day != 1:
                d += delta
                continue
            else:
                dates.append(YEAR + '-' + MONTH)
                d += delta
    if strip_years:
        dates = [d[5:]for d in dates]
        
    return dates


# =========================================================================================== #


#####################
### isnot_outlier ###
#####################

def isNotOutlier(points, thresh=3.5):
    """
    Copied and slightly modified from: 
    https://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
    
    Returns a boolean array with False if points are outliers and True 
    otherwise.

    Parameters:
    -----------
    points : An numobservations by numdimensions array of observations
    thresh : The modified z-score to use as a threshold. Observations with
        a modified z-score (based on the median absolute deviation) greater
        than this value will be classified as outliers.

    Returns:
    --------
    mask : A numobservations-length boolean array.

    References:
    ----------
    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
    Handle Outliers", The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    import numpy as np
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.nanmedian(points, axis=0)
    diff = np.nansum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score < thresh


# =========================================================================================== #


#####################
### findClimoFile ###
#####################

def findClimoFile(regexp, directory = '.'):
    """
    Finds a climo cam h0 file from a directory using glob
    """
    import os, glob, sys
    fname = glob.glob(directory+"/"+regexp)
    try:
        path = os.path.abspath(fname[0])
        name = os.path.splitext(os.path.split(fname[0])[1])[0]
        return [path, name]
    except IndexError:
        print("No file found! Looked for files like " + directory+"/"+regexp)
        #return [0,0]
        raise


# =========================================================================================== #


def niceClev(data, minnstep = 8, maxnstep = 22, alpha = .99, verb = False):
    import numpy as np
    import math

    data = data[~np.isnan(data)]
    numcells = data.size

    MAX = np.max(data)
    MIN = np.min(data)
    noZero = False
    if (MAX < 0) | (MIN > 0):
        if verb:
            print("No way to center color bar on zero")
        return 19
    else:
        delta = MAX - MIN
        test = min(abs(MIN), abs(MAX))
        a = len(np.where(np.logical_and(data>-test, data<test))[0])
        frac = float(a)/float(numcells)
        trace = 0
        while frac < alpha:
            test = test + delta/20
            a = len(np.where(np.logical_and(data>-test, data<test))[0])
            frac = float(a)/float(numcells)
            trace += 1
            if trace > 1000:
                if verb:
                    print("Infinite loop in the niceClev function.\nExiting function with return value of 19")
                return 19

    delta = test*2
    step = (1, 2, 5)
    power_ten = (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5)
    minnstep = minnstep
    maxnstep = maxnstep

    for p in power_ten:
        for s in step:
            nstep = math.ceil(delta/(s * (10**p)))
            if (minnstep <= nstep) & (nstep <= maxnstep) & (nstep % 2 == 0):
                step = s * (10**p)
                return np.linspace(-float(nstep)/2*step, float(nstep)/2*step, num = nstep+1)

    # Default return value if the above nested for loops don't work
    return np.linspace(-test, test, 21)


# =========================================================================================== #


def RegularClev(data, minnstep = 8, maxnstep = 22, alpha = .99, verb = False):
    import numpy as np
    import math

    data = data[~np.isnan(data)]
    numcells = data.size

    

    MAX = np.max(data)
    MIN = np.min(data)
    '''
    delta = MAX - MIN
    test = min(abs(MIN), abs(MAX))
    a = len(np.where(np.logical_and(data>-test, data<test))[0])
    frac = float(a)/float(numcells)
    trace = 0
    while frac < alpha:
        test = test + delta/20
        a = len(np.where(np.logical_and(data>-test, data<test))[0])
        frac = float(a)/float(numcells)
        trace += 1
        if trace > 1000:
            if verb:
                print("Infinite loop in the niceClev function.\nExiting function with return value of 19")
                return 19

    delta = test*2
    step = (1, 2, 5)
    power_ten = (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5)
    minnstep = minnstep
    maxnstep = maxnstep

    for p in power_ten:
        for s in step:
            nstep = math.ceil(delta/(s * (10**p)))
            if (minnstep <= nstep) & (nstep <= maxnstep) & (nstep % 2 == 0):
                step = s * (10**p)
                return np.linspace(-float(nstep)/2*step, float(nstep)/2*step, num = nstep+1)
    '''
    # Default return value if the above nested for loops don't work
    return np.linspace(MIN, MAX, 21)


# =========================================================================================== #

def HadleyCellInfo(psi500, lats):
    '''
    Using Kang and Polvani metrics to find Hadley Cell info
    '''
    import numpy as np
    zonal_mean = np.nanmean(psi500, 1)
    absmax = max(abs(zonal_mean))
    absmax_index = [n for n, m in enumerate(zonal_mean) if abs(m) == absmax][0]
    absmax_psival = zonal_mean[absmax_index]
    '''
    def running_mean(x, N):
        cumsum = np.cumsum(np.insert(x, 0, 0)) 
        return (cumsum[N:] - cumsum[:-N]) / N 

    zonal_mean_smoothed = running_mean(zonal_mean, 5)
    lats = running_mean(lats, 5)

    slope_psi = []
    for n in range(len(zonal_mean_smoothed)):
        del_lat = np.interp(n+.5, range(len(lats)), lats) - np.interp(n-.5, range(len(lats)), lats)
        del_psi = np.interp(n+.5, range(len(zonal_mean_smoothed)), zonal_mean_smoothed) - np.interp(n-.5, range(len(zonal_mean_smoothed)), zonal_mean_smoothed)
        slope_psi.append(del_psi/del_lat)

    slope_psi = np.array(slope_psi)

    absmax = max(abs(zonal_mean_smoothed))
    absmax_index = [n for n, m in enumerate(zonal_mean_smoothed) if abs(m) == absmax][0]
    absmax_psival = zonal_mean_smoothed[absmax_index]
    absmax_slope = slope_psi[absmax_index]

    if absmax_psival > 0:
        if absmax_slope > 0:
            fraction_toZeroSlope = np.interp(0, [slope_psi[absmax_index+1], absmax_slope], [0,1])
            truemax_psi = zonal_mean_smoothed[absmax_index] + ((1 - fraction_toZeroSlope) * (lats[absmax_index+1] - lats[absmax_index]) * slope_psi[absmax_index])
            truemax_lat = np.interp(0, [slope_psi[absmax_index+1], absmax_slope], [lats[absmax_index+1], absmax_slope])
        #
        if absmax_slope < 0:
            fraction_toZeroSlope = np.interp(0, [absmax_slope, slope_psi[absmax_index-1]], [0,1])
            truemax_psi =  zonal_mean_smoothed[absmax_index] - (fraction_toZeroSlope * (lats[absmax_index] - lats[absmax_index-1]) * slope_psi[absmax_index])
            truemax_lat = np.interp(0, [absmax_slope, slope_psi[absmax_index-1]], [lats[absmax_index], lats[absmax_index-1]])
        # Find the first positive slope going north
        first_pos_slope_goingNorth = np.argmax(slope_psi[absmax_index+1:] > 0) + absmax_index
        psi_zeroSlope = np.interp(0, [slope_psi[first_pos_slope_goingNorth-1], slope_psi[first_pos_slope_goingNorth]], [zonal_mean_smoothed[first_pos_slope_goingNorth-1], zonal_mean_smoothed[first_pos_slope_goingNorth]])
        zeroSlope_lat = np.interp(0, [slope_psi[first_pos_slope_goingNorth-1], slope_psi[first_pos_slope_goingNorth]], [lats[first_pos_slope_goingNorth-1], lats[first_pos_slope_goingNorth]])
        # Find the lat that corresponds to that psi value
        l = np.interp((psi_zeroSlope + truemax_psi)/2, [psi_zeroSlope, truemax_psi], [zeroSlope_lat, truemax_lat])
        HadleyCellNorthernEdge = l
        #
        # Find the first negative slope going south
        first_neg_slope_goingSouth = absmax_index - np.argmax(slope_psi[list(reversed(range(0,absmax_index)))] < 0) - 1
        psi_zeroSlope = np.interp(0, [slope_psi[first_neg_slope_goingSouth], slope_psi[first_neg_slope_goingSouth+1]], [psi[first_neg_slope_goingSouth], psi[first_neg_slope_goingSouth+1]])
        zeroSlope_lat = np.interp(0, [slope_psi[first_neg_slope_goingSouth], slope_psi[first_neg_slope_goingSouth+1]], [lats[first_neg_slope_goingSouth], lats[first_neg_slope_goingSouth+1]])
        l = np.interp((psi_zeroSlope + truemax_psi)/2, [psi_zeroSlope, truemax_psi], [zeroSlope_lat, truemax_lat])
        ITCZcenter = l
        #
        prev_pos_slope_goingSouth = first_neg_slope_goingSouth - np.argmax(slope_psi[list(reversed(range(0,first_neg_slope_goingSouth)))] > 0) - 1
        psi_zeroSlope2 = np.interp(0, [slope_psi[prev_pos_slope_goingSouth+1], slope_psi[prev_pos_slope_goingSouth]], [psi[prev_pos_slope_goingSouth+1], psi[prev_pos_slope_goingSouth]])
        zeroSlope_lat2 = np.interp(0, [slope_psi[prev_pos_slope_goingSouth+1], slope_psi[prev_pos_slope_goingSouth]], [lats[prev_pos_slope_goingSouth+1], lats[prev_pos_slope_goingSouth]])
        l = np.interp((psi_zeroSlope + psi_zeroSlope2)/2, [psi_zeroSlope, psi_zeroSlope2], [zeroSlope_lat2, zeroSlope_lat])
        HadleyCellSouthernEdge = l

    if absmax_psival < 0:
        if absmax_slope > 0:
            truemax_psi = np.interp(0, [slope_psi[absmax_index-1], absmax_slope], [zonal_mean_smoothed[absmax_index-1], zonal_mean_smoothed[absmax_index]])
            truemax_lat = np.interp(0, [slope_psi[absmax_index-1], absmax_slope], [lats[absmax_index-1], lats[absmax_index]])
        #
        if absmax_slope < 0:
            truemax_psi =  np.interp(0, [absmax_slope, slope_psi[absmax_index-1]], [zonal_mean_smoothed[absmax_index], zonal_mean_smoothed[absmax_index-1]])
            truemax_lat = np.interp(0, [absmax_slope, slope_psi[absmax_index-1]], [lats[absmax_index], lats[absmax_index-1]])
        #
        # Find the first positive slope going south
        first_pos_slope_goingSouth = absmax_index - np.argmax(slope_psi[list(reversed(range(0,absmax_index)))] > 0) - 1
        psi_zeroSlope = np.interp(0, [slope_psi[first_pos_slope_goingSouth+1], slope_psi[first_pos_slope_goingSouth]], [zonal_mean_smoothed[first_pos_slope_goingSouth+1], zonal_mean_smoothed[first_pos_slope_goingSouth]])
        zeroSlope_lat = np.interp(0, [slope_psi[first_pos_slope_goingSouth+1], slope_psi[first_pos_slope_goingSouth]], [lats[first_pos_slope_goingSouth+1], lats[first_pos_slope_goingSouth]])
        # Find the lat that corresponds to that psi value
        l = np.interp((psi_zeroSlope + truemax_psi)/2, [truemax_psi, psi_zeroSlope], [truemax_lat, zeroSlope_lat])
        HadleyCellSouthernEdge = l
        #
        first_neg_slope_goingNorth = np.argmax(slope_psi[absmax_index+1:] < 0) + absmax_index
        psi_zeroSlope = np.interp(0, [slope_psi[first_neg_slope_goingNorth], slope_psi[first_neg_slope_goingNorth-1]], [zonal_mean_smoothed[first_neg_slope_goingNorth], zonal_mean_smoothed[first_neg_slope_goingNorth-1]])
        zeroSlope_lat = np.interp(0, [slope_psi[first_neg_slope_goingNorth], slope_psi[first_neg_slope_goingNorth-1]], [lats[first_neg_slope_goingNorth], lats[first_neg_slope_goingNorth-1]])
        l = np.interp((psi_zeroSlope + truemax_psi)/2, [truemax_psi, psi_zeroSlope], [truemax_lat, zeroSlope_lat])
        ITCZcenter = l
        #
        next_pos_slope_goingNorth = np.argmax(slope_psi[first_neg_slope_goingNorth+1:] > 0) + first_neg_slope_goingNorth
        psi_zeroSlope2 = np.interp(0, [slope_psi[next_pos_slope_goingNorth-1], slope_psi[next_pos_slope_goingNorth]], [zonal_mean[next_pos_slope_goingNorth-1], zonal_mean[next_pos_slope_goingNorth]])
        zeroSlope_lat2 = np.interp(0, [slope_psi[next_pos_slope_goingNorth-1], slope_psi[next_pos_slope_goingNorth]], [lats[next_pos_slope_goingNorth-1], lats[next_pos_slope_goingNorth]])
        l = np.interp((psi_zeroSlope + psi_zeroSlope2)/2, [psi_zeroSlope, psi_zeroSlope2], [zeroSlope_lat2, zeroSlope_lat])
        HadleyCellNorthernEdge = l
    '''

    if absmax_psival > 0:
        #find the first negative value after the psimax and identify it's index
        #found argmax in a google search, not sure what it really does or why it works
        first_neg_index = np.argmax(zonal_mean[absmax_index:] < 0) + absmax_index
        #calculate the value at and above that index
        first_neg_value = zonal_mean[first_neg_index]
        last_pos_value = zonal_mean[first_neg_index-1]  
        #find the latitudes at those two points
        first_neg_lat = lats[first_neg_index]
        last_pos_lat = lats[first_neg_index-1]  
        #zero should fall in between first_neg_value and last_pos_value
        #calculate the fraction of the psi gap between first_neg and last_pos that 0 falls in 
        f = abs(last_pos_value) / abs(last_pos_value - first_neg_value) 
        #calculate the latitude using f and assuming lat changes are linear in this lat gap
        zerolat = last_pos_lat + f * abs(first_neg_lat - last_pos_lat)  
        #this latitude is the northern edge of the hadley cell
        HadleyCellNorthernEdge = zerolat    
        ################################################################################################
        # now go backwards from the absmax_index and find the next zero-crossing to find the ITCZ center
        # following the same basic algorithm as above
        ################################################################################################
        first_neg_index = absmax_index - np.argmax(zonal_mean[list(reversed(range(0,absmax_index)))] < 0)   
        first_neg_value = zonal_mean[first_neg_index]
        last_pos_value = zonal_mean[first_neg_index+1]
        first_neg_lat = lats[first_neg_index]
        last_pos_lat = lats[first_neg_index+1]  
        f = abs(last_pos_value) / abs(last_pos_value - first_neg_value) 
        #going up the data frame (decreasing latitude), so we subtract the fractional latitude change from the last positive latitude
        zerolat = last_pos_lat - f * abs(first_neg_lat - last_pos_lat)
        #this latitude is the center of the ITCZ
        ITCZcenter = zerolat    
        #################################################################################################
        # now go backwards from the ITCZ center and find the next zero-crossing to find the SH Hadley 
        # edge. Same basic algorithm. 
        #################################################################################################
        first_pos_index = first_neg_index - np.argmax(zonal_mean[list(reversed(range(0,first_neg_index)))] > 0) 
        first_pos_value = zonal_mean[first_pos_index]
        last_neg_value = zonal_mean[first_pos_index+1]
        first_pos_lat = lats[first_pos_index]
        last_neg_lat = lats[first_pos_index+1]  
        f = abs(last_neg_value) / abs(last_neg_value - first_pos_value) 
        #going up the data frame so subtract the fractional latitude change from the last positive latitude
        zerolat = last_neg_lat - f * abs(first_neg_lat - last_pos_lat)
        HadleyCellSouthernEdge = zerolat 
    ##############################################
    # Do practically the same thing for absmax < 0
    # Just kind of flip the operation around
    ##############################################
    if absmax_psival < 0:
        first_pos_index = absmax_index - np.argmax(zonal_mean[list(reversed(range(0,absmax_index)))] > 0)
        first_pos_value = zonal_mean[first_pos_index]
        last_neg_value = zonal_mean[first_pos_index+1]
        first_pos_lat = lats[first_pos_index]
        last_neg_lat = lats[first_pos_index+1]  
        f = abs(last_neg_value) / abs(last_neg_value - first_pos_value) 
        #going up the data frame so subtract the fractional latitude change from the last positive latitude
        zerolat = last_neg_lat - f * abs(last_neg_lat - first_pos_lat)
        HadleyCellSouthernEdge = zerolat    
        ################################################################################################
        # now go forwards from the absmax_index and find the next zero-crossing to find the ITCZ center
        ################################################################################################
        first_pos_index = np.argmax(zonal_mean[absmax_index:] > 0) + absmax_index
        first_pos_value = zonal_mean[first_pos_index]
        last_neg_value = zonal_mean[first_pos_index-1]
        first_pos_lat = lats[first_pos_index]
        last_neg_lat = lats[first_pos_index-1]  
        f = abs(last_neg_value) / abs(last_neg_value - first_pos_value) 
        zerolat = last_neg_lat + f * abs(last_neg_lat - first_pos_lat)
        ITCZcenter = zerolat
        ################################################################################################
        # now go forwards from the ITCZ center and find the next zero-crossing to find the NH Hadley edge
        ################################################################################################
        first_neg_index = np.argmax(zonal_mean[first_pos_index:] < 0) + first_pos_index
        first_neg_value = zonal_mean[first_neg_index]
        last_pos_value = zonal_mean[first_neg_index-1]
        first_neg_lat = lats[first_neg_index]
        last_pos_lat = lats[first_neg_index-1]  
        f = abs(first_neg_value) / abs(first_neg_value - last_pos_value)    
        zerolat = last_pos_lat + f * abs(last_pos_lat - first_neg_lat)
        HadleyCellNorthernEdge = zerolat

    return zonal_mean, {"ITCZcenter" : ITCZcenter, "HadleyCellNorthernEdge" : HadleyCellNorthernEdge, "HadleyCellSouthernEdge": HadleyCellSouthernEdge}


# =========================================================================================== #
def ensoDates(sy, ey, directory, returnSSTs = False):
    '''
    sy and ey are start year and end year. 

    returns: a list of two lists. index 0 is the nino dates, index 1 is the nina dates.

    the method for determining nino/nina dates is from:
    Trenberth, K. E. (1997) The Definition of El Nino. Bulletin of the American Meteorological Society, 78, 2771-2777.

    basically, compute nino3.4 sst monthly anomalies, compute the 5-month running mean of these anomalies, then look for
    instances when the 5-month running mean is > .3 deg c (nino) or < -.3 deg c (nina) for 6 months in a row at least.
    '''
    from pygoda import popgoda, camdates, findClimoFile, runningMean
    import numpy as np

    dates = camdates(sy, ey)

    nino34 = []
    for d in dates:
        f, fn = findClimoFile("*pop.h."+d+".nc", directory = directory)
        print(fn)
        pop = popgoda(f)
        nino34.append(pop.boxMean("TEMP", box = (-5,5,190,240)))


    anom_nino34 = nino34 - np.mean(nino34)
    rm5_anom_nino34 = runningMean(anom_nino34, 5)

    # Create a dummy list of the status of the nino34 anomaly running mean based on whether it is above or below .3 deg c, or neither
    status = []
    for t in rm5_anom_nino34:
        if t > .3:
            status.append('above')
        elif t < -.3:
            status.append('below')
        else:
            status.append('neither')

    # Cycle over the status list and look for instance of 6-in-a-row of the same
    current_status = ''
    counter = 1
    upper_idx = []
    count_idx = []
    for n, s in enumerate(status):
        if s is current_status:
            counter += 1
        else:
            if counter >= 6:
                upper_idx.append(n-1)
                count_idx.append(counter-1)
            counter = 1
        current_status = s

    # Do it one last time in case the last months are included
    if counter >= 6:
        upper_idx.append(n-1)
        count_idx.append(counter-1)

    lower_idx = np.array(upper_idx) - np.array(count_idx)

    # Now, find the nino and nina dates that correspond to these anomalies
    nino_dates = []
    nina_dates = []
    for n in range(len(lower_idx)):
        l = lower_idx[n]
        u = upper_idx[n]
        s = status[l]
        if s == 'above':
            nino_dates.append(dates[l:u+1])
        if s == 'below':
            nina_dates.append(dates[l:u+1])

    # Flatten the lists
    nino_dates = [item for sublist in nino_dates for item in sublist]
    nina_dates = [item for sublist in nina_dates for item in sublist]

    if returnSSTs:
        print("If all went well, the returned list is of nino then nina dates, then the raw SST timeseries")
        return [nino_dates, nina_dates, nino34]
    else:
        print("If all went well, the returned list is of nino then nina dates")
        return [nino_dates, nina_dates]


def nino34(ts_data, time_axis = 0, lat_axis = -2, lon_axis = -1, grid = "2deg"):
    from pygoda import boxOut
    import numpy as np
    nino34_box = (-5,5,190,240)
    data_nino34Box = boxOut(data, nino34_box, lat_axis = lat_axis, lon_axis = lon_axis, grid = grid)
    data_nino34Box_avg = np.nanmean(data_nino34Box, axis = (lat_axis, lon_axis))



###############
### camgoda ###
###############
    
class camgoda:
    """
*** Class camgoda ***

Class useful for working with iCAM5 history files. Instances are initiated with .dataset, .lat, and .lon
from the netcdf dataset. 

====================
!!!Need to update!!!
====================
1) Right now, class is not safely capable of handling time dimensions with length > 1.
That means camgoda should be used to open .nc files with only 1 time slice.

Not the best way to do it, and should incorporate time dimensions into isobar and variable methods for completion.

Right now, best way to deal with this is to just loop over files with 1 time slice and create a new
camgoda instance each time. 

---

2) Incorporate way to slice model level rather than isobar.

Along with this, might want to make a way to show all the 3d data rather than a slice.

---
=============
How it works:
=============

Initialize:
-----------
Arugments : [0] dataset (netCDF4 Dataset)

Methods:
--------
variable : returns 2 or 3d numpy array data. 
    Sets : .var, .box, .data, .boxlat, .boxlon, .units, and .long_name
    Arugments : [0] string variable name; [1] region box (oprtional) in form of [bottom, top, left, right]

isobar : returns 2d np array along an isobar
    Sets : .pressure, .data
    Arguments : [0] string variable name; [1] integer pressure in Pa; [2] region box (optional)

PRECT_d18O and PRECT_dD : returns 2d numpy array data. 
    Sets : .h2o, .h218o/.hdo, .var, .data, and .long_name
    Arguments : [0] region box (optional)

d18OV and dDV : returns 2d numpy array data. 
    Sets : .h2ov_p, .h218ov_p/.hdov_p, .var, .data, and .long_name
    Arguments : [0] integer pressure in Pa; [1] region box (optional)

    """
    def __init__(self, dataset):
        from netCDF4 import Dataset as ds
        self.dataset = dataset
        if type(dataset) is str:
            self.dataset = ds(dataset, "r")
        self.dims = [d for d in self.dataset.dimensions]
        self.vars = [v for v in self.dataset.variables]
        self.dimlen = [len(self.dataset.dimensions[d]) for d in self.dims]
        self.lat = self.dataset.variables['lat'][:]
        self.lon = self.dataset.variables['lon'][:]
        self.boxlat = self.lat
        self.boxlon = self.lon
        
        self.model = "CAM"
        if any(v == "levgrnd" for v in self.vars):
            self.model = "CLM"
            self.depths = self.dataset.variables['levgrnd'][:]
        if any(v == "z_t" for v in self.vars):
            self.model = "POP"
            self.depths = self.dataset.variables['z_t'][:]

        # Some variables that don't need to be caluclated more than once
        self.PressureCalculated = False
        self.PsiCalculated = False
        
        self.var = ''
        self.units = ''
        self.long_name = ''


        #self.levs = self.dataset.variables['lev'][:]
        self.isTime = any(d == "time" for d in self.dims)
        if self.isTime:
            self.ntime = self.dimlen[self.dims.index("time")]
            if self.ntime == 1:
                self.isTime = False
    
    def close(self):
        self.dataset.close()
    
    def setBox(self, box):
        from pygoda import find_indices
        self.box = box
        idxbox = find_indices(box, self.lat, self.lon)
        self.boxlat = self.lat[idxbox[0]]
        self.boxlon = self.lon[idxbox[1]]
        
    def variable(self, var, box = None, verb = False, setData = True, math = True, idx_add = 0):
        '''
        Main function for extracting data from netcdf file
        Performs conversions on data if there is documentation for it in algebra.py
        '''
        from pygoda import find_indices
        import xarray as xr
        import numpy as np

        # Read the variable 
        DATA = self.dataset.variables[var]
        data = DATA[:].squeeze()
        # Set some variables for the whole class to see
        self.box = box
        if setData:
            self.vartype = "2d"
            self.var = var
            try:
                self.units = DATA.units
            except AttributeError:
                self.units = var

        ndims = len(data.shape)
        
        # Strip the numbers down to the region from the box arg
        idxbox = find_indices(box, self.lat, self.lon)

        # Add and subtract extra values from the indices according to idx_add
        # This is used for calculating fluxes and divergences, where boundarie become a problem
        # If we are not at the edge of the world, it is best to add an extra lat and lon on the left and right 
        # in order to calculate fluxes and divergence at boundaries
        
        # This is causing problems when the differences between lats and lons is small (e.g., if box = (15,15,100,100)
        # Need to fix this still, until then, this will remain commented out.
        '''
        if (0 not in idxbox[0]) and ((len(self.lat)-1) not in idxbox[0]):
            idxbox[0] = np.append(idxbox[0][0]-idx_add, idxbox[0], idxbox[0][-1]+idx_add)
        if (0 not in idxbox[1]) and ((len(self.lon)-1) not in idxbox[1]):
            idxbox[1] = np.append(idxbox[1][0]-idx_add, idxbox[1], idxbox[1][-1]+idx_add)
        '''
        dataDims = "Undefined"
        # No lev, no time
        dataDims = "N/A"
        if ndims == 2:
            dataDims = ('lat', 'lon')
            data = np.array(xr.DataArray(data)[idxbox[0], idxbox[1]])
        # Either lev or time but not both
        if ndims == 3:
            if not self.isTime:
                dataDims = ('z', 'lat', 'lon')
                if setData:
                    self.vartype = "3d"
            data = np.array(xr.DataArray(data)[:, idxbox[0], idxbox[1]])
        if ndims == 4:
            dataDims = ("time", "z", "lat", "lon")
            if setData:
                self.vartype = "3d"
            data = np.array(xr.DataArray(data)[:, :, idxbox[0], idxbox[1]])

        # Do conversion (default of no conversion if there isn't documentation for a conversion in algebra.py)
        if math:
            try:
                from algebra import algebra
                # Convert the data using algebra.py file
                # Set new units to correspond to the algebra changes
                mult, add, units = algebra.get(var, [1, 0, self.units]) # default to no change if var is not found
                data = (data * mult) + add
                if setData:
                    self.units = units
                if verb:
                    if (mult != 1) | (add != 0):
                        print("From pygoda.variable(): Converted data as follows:\n\tmultiplied by " + str(mult) + "\n\tadded " + str(add))
            except ImportError:
                print("Could not find algebra.py, not able to do conversion!")

        self.boxlat = self.lat[idxbox[0]]
        self.boxlon = self.lon[idxbox[1]]
        if setData:
            self.long_name = DATA.long_name
            self.dataDims = dataDims
            self.data = data
        return data

    def CalculatePressure(self, box = None, setData = False, idx_add = 0):
        from pygoda import PressureCalc
        import numpy as np
        A_m = self.variable("hyam", box = box, setData = setData, idx_add = idx_add)
        A_i = self.variable("hyai", box = box, setData = setData, idx_add = idx_add)
        B_m = self.variable("hybm", box = box, setData = setData, idx_add = idx_add)
        B_i = self.variable("hybi", box = box, setData = setData, idx_add = idx_add)
        PS = self.variable("PS", box = box, setData = setData, idx_add = idx_add)
        # Create 3d pressure field
        self.P_m = PressureCalc(A_m,B_m,PS)
        self.P_i = PressureCalc(A_i,B_i,PS)
        pdel = np.zeros(shape = self.P_m.shape)
        for i in range(pdel.shape[1]):
            for j in range(pdel.shape[2]):
                pdel[:,i,j] = np.diff(self.P_i[:,i,j])
        self.Pdel = pdel
        self.PS = PS
        self.PressureCalculated = True

    def isobar(self, pressure, setData = True, verb = False):
        if self.var == '':
            print("No variable read in yet.\nUse self.variable() to read a variable first.")
            return
        if self.vartype == "2d":
            print("Current .data attribute is 2-dimensional, no pressure data.\nRead in a 3d variable first before running this method.")
            return
        import numpy as np
        data = self.data
        if not self.PressureCalculated:
            self.CalculatePressure(self.box)

        VAR = np.zeros(shape = (len(self.boxlat), len(self.boxlon)))
        for i in range(len(self.boxlat)):
            for j in range(len(self.boxlon)):
                VAR[i,j] = np.interp(pressure, P_m[:,i,j], data[:,i,j])

        self.pressure = pressure
        if setData:
            self.data = VAR
            self.long_name += " @" + str(int(pressure/100)) + "mb"

        return VAR

    def depth(self, d, setData = True, verb = False):
        '''
        Only works when the camgoda instance is of a CLM file.
        d should be in cm
        '''
        if self.model == "CAM":
            print("Can not run camgoda.depth() on a CAM file! Read in a CLM or POP file to run this method.")
            return
        if self.var == '':
            print("No variable read in yet.\nUse self.variable() to read a variable first.")
            return
        if self.vartype == "2d":
            print("Current .data attribute is 2-dimensional, no depth data.\nRead in a 3d variable first before running this method.")
            return
        import numpy as np
        d_m = d
        if self.model == "CLM":
            d_m = float(d) / 100. # Convert d to meters, which is the unit of the model
        data = self.data
        depths = self.depths
        VAR = np.zeros(shape = (len(self.boxlat), len(self.boxlon)))
        for i in range(len(self.boxlat)):
            for j in range(len(self.boxlon)):
                VAR[i,j] = np.interp(d_m, depths, data[:,i,j])

        if setData:
            self.data = VAR
            self.long_name += " @" + str(d) + "cm"

        return VAR
    
    def isotherm(self, val, setData = True):
        import numpy as np
        import scipy.interpolate as interp
        ret = np.zeros(shape = (len(self.boxlat), len(self.boxlon)))
        for i in range(len(self.boxlat)):
            for j in range(len(self.boxlon)):
                current_ij = self.data[:,i,j]
                mask = np.where(current_ij > 1e36, np.nan, 1)
                current_ij_masked = current_ij * mask
                if sum(~np.isnan(current_ij_masked)) == 0:
                    ret[i,j] = np.nan
                else:
                    f = interp.interp1d(current_ij_masked,self.depths)
                    try:
                        ret[i,j] = f(val)
                    except ValueError:
                        ret[i,j] = np.nan
        if setData:
            self.data = ret
            self.units = "cm"
            self.long_name = "depth of 15 degC isotherm"
        return ret
    
    # def isotherm_fast(self, temp, setData = True):
    #   import numpy as np
    #   if self.model is not "POP":
    #       print("Not a POP file, can not run isotherm_fast! Exiting...")
    #       return None
    #   d = self.data
    #   depths = self.depths
    #   y = np.zeros(shape = (d.shape[1], d.shape[2], 2))
    #   z = np.zeros(shape = (d.shape[1], d.shape[2], 2))
    #   for i in range(d.shape[1]):
    #       for j in range(d.shape[2]):
    #           dij = d[:,i,j]
    #           dij_mask = np.where(dij > 1e36, np.nan, 1)
    #           dij_masked = dij*dij_mask
    #           if sum(~np.isnan(dij_masked)) == 0:
    #           y[i,j,:] = np.nan
    #           z[i,j,:] = np.nan
    #           continue
    #           dij_nonan = dij_masked[~np.isnan(dij_masked)]
    #           try:
    #           y1 = dij_nonan[dij_nonan-15<0][0]
    #           y2 = dij_nonan[dij_nonan-15>0][-1]
    #           z1 = depths[np.where(dij_nonan == y1)[0][0]]
    #           z2 = depths[np.where(dij_nonan == y2)[0][0]]
    #           y[i,j,0] = y1
    #           y[i,j,1] = y2
    #           z[i,j,0] = z1
    #           z[i,j,1] = z2
    #           except IndexError:
    #           y[i,j,:] = np.nan
    #           z[i,j,:] = np.nan
    #   y1 = y[...,0]
    #   y2 = y[...,1]
    #   z1 = z[...,0]
    #   z2 = z[...,1]
    #   ret = ((temp-y1) / ((y1-y2)/(z1-z2))) + z1
    #   if setData:
    #       self.data = ret
    #       self.units = "cm"
    #       self.long_name = "depth of " + str(temp) + " degC isotherm"
    #   return ret


    def columnSum(self, box = None, setData = True):
        if self.vartype != "3d":
            print("Var type is not 3d, cannot compute column sum! Exiting...")
            return
        import numpy as np
        g = 9.8
        if not self.PressureCalculated:
            self.CalculatePressure(box, setData = False)
        # Weight the data by the absolute amount of mass in each layer
        csum = self.Pdel * self.data / g
        csum_vertsum = np.sum(csum, axis = 0)
        if setData:
            self.var += "_columnSum"
            self.long_name += ", column sum"
            self.vartype = "2d"
            self.data = csum_vertsum
        return csum_vertsum

    def columnMean(self, box = None, setData = True, idx_add = 0):
        if self.vartype != "3d":
            print("Var type is not 3d, cannot compute column sum! Exiting...")
            return
        import numpy as np
        if not self.PressureCalculated:
            self.CalculatePressure(box, setData = False, idx_add = idx_add)
        # Weight the data by the relative amount of mass in each layer
        cavg = self.Pdel / self.PS * self.data
        cavg_vertsum = np.sum(cavg, axis = 0)
        if setData:
            self.var += "columnMean"
            self.long_name += ", column mean"
            self.vartype = "2d"
            self.data = cavg_vertsum
        return cavg_vertsum

    def prep_map(self, season, region):
        try:
            from clevs import getClev, getCmap
            import matplotlib.pyplot as plt
            # pressure = ""
            # if self.vartype == "3d":
            #   pressure = self.pressure
            #self.clevs = getClev(self.var+str(pressure), season+"_"+region, self.data)
            self.cmap = getCmap(self.var, 'cmap')
            self.diffcmap = getCmap(self.var, 'diffcmap')
        except ImportError:
            return "Could not find clevs.py, not able to set clevs!"
            
    def mask(self, array, gt_lt, val, targetArray = None):
        '''
        This function is intended to mask an array for bad data.
        It can also be used for taking out small values of, e.g., QFLX, when used in the 
        denominator of ratios. 

        Viable options for 'gt_lt' argument are 'gt' (greater than) or 'lt' (less than) as of now.

        This method does not set ANY data, including self.data, self.long_name, self.units, or anything like that!

        It requires xarray.

        The basic idea is to take the 'array', and find all values
        that are 'bool' (less than or greater than) when compared to 'val'.
        
        Example:
        array = [[1,5,3,2,6,8,2345,3456,2,5,8,4],
                 [6,3,7,4,3,6,6,5,3,5,11234,6],
                 [6,4,3452345,3,6,5,3,5,67,8,6,4]]
        We want to remove the large values: 2345, 3456, 11234, 3452345, and 67
        as they are clearly outliers. To do so, we would run the function with the following arguments:

        masked_array = self.mask(array, 'gt', 50)

        This method also provides functionality for removing values from target arrays based on
        a mask generated from array, but most of the time this is not necessary. 

        To mask an array called targetArray with a mask generated from array, run

        masked_targetArray = self.mask(array, 'gt', 50, targetArray)
        '''
        if (gt_lt != "gt") and (gt_lt != "lt"):
            print("First argument (gt_lt) must be either 'gt' or 'lt'...Exiting method.")
            return 
        import xarray as xr
        import numpy as np
        
        # This is the array to be tested against the boolean expression in where()
        maskArray = xr.DataArray(array)

        # This is the array the will be masked and returned by the method
        if targetArray is not None:
            if targetArray.shape != array.shape:
                print("Target array and input array are not the same shape.\nCan not mask a target array of different shape than the input array!\nExiting...")
                return
            returnArray = xr.DataArray(targetArray)
        else:
            returnArray = xr.DataArray(array)
        if gt_lt == 'gt':
            # This is supposed to be the less than sign even though we're in the gt if statement.
            # This is because the masking function in xarray keeps values that 
            # satisfy the boolean in the where argument, and makes everything else nan
            ret = np.array(returnArray.where(maskArray < val))
        else:
            ret = np.array(returnArray.where(maskArray > val))
        return ret
    
    def PRECT_d18O(self, box = None):
        h2o = self.variable("PRECT_H2O", box, setData = False)
        h218o = self.variable("PRECT_H218O", box, setData = False)
        self.var = "PRECT_d18O"
        self.long_name = "d18O of PRECT"
        self.units = "permil"
        self.vartype = "2d"
        self.dataDims = ("lat", "lon")
        self.data = (h218o / h2o - 1) * 1000
        # Filter diriculous data
        self.data = self.mask(self.data, 'gt', 50)
        self.data = self.mask(self.data, 'lt', -50)
        return self.data

    def PRECT_dD(self, box = None):
        h2o = self.variable("PRECT_H2O", box, setData = False)
        hdo = self.variable("PRECT_HDO", box, setData = False)
        self.var = "PRECT_dD" 
        self.long_name = "dD of PRECT"
        self.units = "permil"
        self.vartype = "2d"
        self.dataDims = ("lat", "lon")
        self.data = (hdo / h2o - 1) * 1000
        return self.data

    def PRECT_dxs(self, box = None):
        dd = self.PRECT_dD(box)
        d18o = self.PRECT_d18O(box)
        self.var = "PRECT_dxs"
        self.long_name = "d-excess of PRECT"
        self.units = "permil"
        self.vartype = "2d"
        self.dataDims = ("lat", "lon")
        self.data = dd - 8 * d18o
        return self.data

    def d18OV(self, box = None):
        h218ov_p = self.variable("H218OV", box, setData = False) 
        h2ov_p = self.variable("H2OV", box, setData = False)
        self.var = "d18OV"
        self.long_name = "d18O of Vapor"
        self.units = "permil"
        self.vartype = "3d"
        self.dataDims = ("z", "lat", "lon")
        self.data = (h218ov_p / h2ov_p - 1) * 1000
        return self.data

    def dDV(self, box = None):
        hdov_p = self.variable("HDOV", box, setData = False) # only talk once
        h2ov_p = self.variable("H2OV", box, setData = False)
        self.var = "dDV" 
        self.long_name = "dD of Vapor"
        self.units = "permil"
        self.vartype = "3d"
        self.dataDims = ("z", "lat", "lon")
        self.data = (hdov_p / h2ov_p - 1) * 1000
        return self.data

    def dxsV(self, box = None):
        ddv_p = self.dDV(box)
        d18ov_p = self.d18OV(box)
        self.var = "dxsV" 
        self.long_name = "d-excess of VAPOR"
        self.units = "permil"
        self.vartype = "3d"
        self.dataDims = ("z", "lat", "lon")
        self.data = ddv_p - 8 * d18ov_p
        return self.data

    def VQ_d18O(self, box = None):
        import numpy as np
        vq_h218o = self.variable("VQ_H218O", box, setData = False)
        vq_h2o = self.variable("VQ_H2O", box, setData = False)
        vq_h2o, vq_h218o = [self.mask(np.abs(vq_h2o), 'lt', 0.001, x) for x in [vq_h2o, vq_h218o]]
        self.var = "VQ_d18O"
        self.long_name = "VQ_d18O"
        self.units = "m/s * delta"
        self.vartype = "3d"
        self.dataDims = ("z", "lat", "lon")
        data = (vq_h218o / vq_h2o - 1) * 1000
        self.data = data
        return self.data

    # def VQ_dD(self, box = None):
    #   vq_hdo = self.variable("VQ_HDO", box, setData = False)
    #   vq_h2o = self.variable("VQ_H2O", box, setData = False)
    #   self.var = "VQ_dD"
    #   self.long_name = "VQ_dD"
    #   self.units = "m/s * delta"
    #   self.vartype = "3d"
    #   self.data = (vq_hdo / vq_h2o - 1) * 1000
    #   return self.data

    # def UQ_d18O(self, box = None):
    #   uq_h218o = self.variable("UQ_H218O", box, setData = False)
    #   uq_h2o = self.variable("UQ_H2O", box, setData = False)
    #   self.var = "UQ_d18O"
    #   self.long_name = "UQ_d18O"
    #   self.units = "m/s * delta"
    #   self.vartype = "3d"
    #   self.data = (uq_h218o / uq_h2o - 1) * 1000
    #   return self.data

    # def UQ_dD(self, box = None):
    #   uq_hdo = self.variable("UQ_HDO", box, setData = False)
    #   uq_h2o = self.variable("UQ_H2O", box, setData = False)
    #   self.var = "UQ_dD"
    #   self.long_name = "UQ_dD"
    #   self.units = "m/s * delta"
    #   self.vartype = "3d"
    #   self.data = (uq_hdo / uq_h2o - 1) * 1000
    #   return self.data

    def QFLX_d18O(self, box = None):
        import xarray as xr
        import numpy as np
        qflx_h218o = self.variable("QFLX_H218O", box, setData = False, math = False)
        qflx_h2o = self.variable("QFLX_H2O", box, setData = False, math = False)
        qmask = xr.DataArray(qflx_h2o)
        qflx_h2o = np.array(qmask.where(abs(qmask) > .00000008))
        qflx_h218o = np.array(xr.DataArray(qflx_h218o).where(abs(qmask) > .00000008))
        # mask out large values of qflx_h218o
        qmask = xr.DataArray(qflx_h218o)
        qflx_h218o = np.array(qmask.where(abs(qmask) < 100))
        qflx_h2o = np.array(xr.DataArray(qflx_h2o).where(abs(qmask) < 100))
        self.var = "QFLX_d18O"
        self.long_name = "d18O of QFLX"
        self.units = "delta 18O"
        self.vartype = "2d"
        self.dataDims = ("lat", "lon")
        self.data = (qflx_h218o / qflx_h2o - 1) * 1000
        return self.data

    def QFLX_dD(self, box = None):
        import xarray as xr
        import numpy as np
        qflx_hdo = self.variable("QFLX_HDO", box, setData = False, math = False)
        qflx_h2o = self.variable("QFLX_H2O", box, setData = False, math = False)
        qmask = xr.DataArray(qflx_h2o)
        qflx_h2o = np.array(qmask.where(abs(qmask) > .00000005))
        qflx_hdo = np.array(xr.DataArray(qflx_hdo).where(abs(qmask) > .00000005))
        self.var = "QFLX_dD"
        self.long_name = "dD of QFLX"
        self.units = "delta D"
        self.vartype = "2d"
        self.dataDims = ("lat", "lon")
        self.data = (qflx_hdo / qflx_h2o - 1) * 1000
        return self.data

    def convergenceDelta(self, box = None):
        import numpy as np
        P = self.variable("PRECT_H2O", box = box)
        E  = self.variable("QFLX_H2O", box = box)
        RP = self.variable("PRECT_H218O", box = box)
        RE = self.variable("QFLX_H218O", box = box)

        # Mask everything by the weird values of RE near ice
        RE, RP, P, E = [self.mask(RE, 'lt', 1, x) for x in [RE, RP, P, E]]
        P_E = P-E

        # Mask for near-zero values of P-E
        RE, RP, P_E = [self.mask(np.abs(P_E), 'lt', 0.1, x) for x in [RE, RP, P_E]]

        data = ((RP-RE)/(P_E) - 1) * 1000
        self.var = "convergenceDelta"
        self.long_name = "d18O of moisture convergence"
        self.units = "delta 18O"
        self.vartype = "2d"
        self.dataDims = ("lat", "lon")
        self.data = data
        return data

    def psi(self, box = None):
        import numpy as np
        g = 9.8 #m/s2
        radius = 6371000 # radius of the earth in m
        V = self.variable("V", box = box, setData = False)
        nlats = len(self.boxlat)
        nlons = len(self.boxlon)
        if not self.PressureCalculated:
            self.CalculatePressure(box = box)
        psi = np.zeros(shape = [31, V.shape[1], V.shape[2]])
        for z in range(nlats):
            coslat = np.cos(np.radians(self.boxlat[z]))
            for m in range(nlons):
                P_hold = self.P_i[:,z,m]
                V_hold = V[:,z,m]
                del_psi = 0
                psi_hold = [0]
                for p in range(30): 
                    del_P = P_hold[p+1] - P_hold[p]
                    del_psi = V_hold[p] * del_P * 2*np.pi*radius * coslat / g
                    psi_hold.append(psi_hold[-1] + del_psi)
                    psi[p+1,z,m] = psi_hold[-1]

        self.long_name = "Top-down integrated psi"
        self.units = "Kg/s"
        self.var = "psi"
        self.vartype = "3d"
        self.dataDims = ("z", "lat", "lon")
        self.data = psi
        return self.data

    def divergence(self, box = None):
        u = self.variable("U", box)
        u_colMean = self.columnMean(box)
        v = self.variable("V", box)
        v_colMean = self.columnMean(box)
        return # not finished yet



    def RH(self, box = None):
        import numpy as np
        if not self.PressureCalculated:
            self.CalculatePressure(box = box)
        E = self.variable("Q", box = box, setData = False) / 1000 * self.P_m / 0.622
        Es = 610.78 * np.exp(self.variable("T", box = box, setData = False) / (self.variable("T", box = box, setData = False) + 238.3) * 17.2694)
        RH = E/Es * 100
        RH = self.mask(RH, 'gt', 120)
        self.data = RH
        self.long_name = "Relative Humidity"
        self.units = "%"
        self.var = "rh"
        self.vartype = "3d"
        self.dataDims = ("z", "lat", "lon")
        return RH
    
    def HadleyCellInfo(self, box = None):
        from pygoda import HadleyCellInfo
        self.psi(box = box)
        self.Hadley = HadleyCellInfo(self.isobar(50000), self.boxlat)
        return self.Hadley

    def LandMask(self, array = None):
        import numpy as np
        if array is None:
            array = self.data
        landfrac = self.variable("LANDFRAC", box = self.box, setData = False)
        return array*np.where(landfrac > .9, 1, np.nan)
    
    def OceanMask(self, array = None):
        import numpy as np
        if array is None:
            array = self.data
        ocnfrac = self.variable("OCNFRAC", box = self.box, setData = False)
        return array*np.where(ocnfrac > .9, 1, np.nan)
        
    def ExtractData(self, var_names, box = None, returnData = False):
        '''
        The main function for extracting data, everything else is just behind the scenes stuff
        Args: 
        V - variables to be read in, separated by a comma if more than one
        box - bottom lat (-90,90), top lat (-90,90), left lon (0,360), right lon (0,360)
        returnData - default to False but usually used as True, need to fix legacy code in order for default to be changed to true
        Returns:
        list of 3 is returnData = False, self.data if returnData = True
        Always set returnData = True if V is more than one variable!
        '''
        import numpy as np
        variables = [v.strip() for v in var_names.split(",")]
        self.setBox(box)
        RETURN = np.zeros(shape = (len(self.boxlat), len(self.boxlon), len(variables))).squeeze()
        for V_idx, V in enumerate(variables):
            var_is_3d = False
            var_is_ColumnMean = False
            if V[:3] == '3d_':
                var_is_3d = True
            if not var_is_3d:
                var = V
                vname = V
                pressure = None
            else:
                split = V.split("_")
                if 3 > len(split) > 4:
                    print("Something went wrong parsing the variable name.\nMake sure there are two or three underscores!")
                    return 
                if split[-1] == "ColumnMean":
                    var_is_3d = False
                    var_is_ColumnMean = True
                else:
                    pressure = float(split[-1]) * 100
                var = split[1]
                vname = split[1]
                if len(split) == 4: # Handle variables like VQ_H2O
                    var += "_" + split[2]
                    vname += "_" + split[2]
                vname += split[-1]
            # Extract the variable data
            # Special variables
            if var == "PRECT_d18O":
                self.PRECT_d18O(box)
            elif var == "PRECT_dD":
                self.PRECT_dD(box)
            elif var == "PRECT_dxs":
                self.PRECT_dxs(box)
            elif var == "QFLX_d18O":
                self.QFLX_d18O(box)
            elif var == "QFLX_dD":
                self.QFLX_dD(box)
            elif var == "fluxDelta":
                self.fluxDelta(box)
            elif var == "Column_d18OV":
                self.variable('H2OV', box)
                denom = self.columnSum(box, setData = False)
                self.variable('H218OV', box)
                num  = self.columnSum(box, setData = False)
                self.data = (num/denom - 1) * 1000
                self.long_name = "Column d18OV"
                self.units = "delta 18O"
                self.dataDims = ("lat", "lon")
            elif var == "Column_dDV":
                self.variable('Q', box)
                denom = self.columnSum(box, setData = False)
                self.variable('HDOV', box)
                num  = self.columnSum(box, setData = False)
                self.data = (num/denom - 1) * 1000
                self.long_name = "Column dDV"
                self.dataDims = ("lat", "lon")
                self.units = "delta D"
            elif var == "Column_Q":
                self.variable("Q", box)
                self.columnSum(box)
                self.long_name = "Vertically integrated water"
                self.units = "g-H2O"
                self.dataDims = ("lat", "lon")
            elif var == "P_E":
                data = self.variable('PRECT', box) - self.variable('QFLX', box)
                self.data = data
                self.units = "kg/m2/day"
                self.long_name = "Moisture convergence (P-E)"
                self.dataDims = ("lat", "lon")
            elif var == "PoverE":
                import numpy as np
                prect = self.variable('PRECT', box) # units of kg/m2/day
                qflx = self.variable('QFLX', box) # kg/m2/day
                qflx = self.mask(qflx, 'lt', .005) # remove qflx less than 5 g/m2/day ---> pretty much zero
                self.data = prect/qflx
                self.units = "unitless (ratio)"
                self.long_name = "Ratio of precipitation to evaporation P/E"
                self.dataDims = ("lat", "lon")
            elif var == "convergenceDelta":
                self.convergenceDelta(box)
            elif var == "d18OV":
                self.d18OV(box)
            elif var == "dDV":
                self.dDV(box)
            elif var == "dxsV":
                self.dxsV(box)
            elif var == "psi":
                self.psi(box)
            elif var == "RH":
                self.RH(box)
            elif var == "VQ_d18O":
                self.VQ_d18O(box)
            elif var == "VQ_dD":
                self.VQ_dD(box)
            elif var == "UQ_d18O":
                self.UQ_d18O(box)
            elif var == "UQ_dD":
                self.UQ_dD(box)
            elif var == "QFLX_d18O":
                self.QFLX_d18O(box)
            elif var == "VQ_d18O":
                self.VQ_d18O(box)
            elif var == "PRECTsquared":
                self.variable("PRECT", box)
                self.data = self.data**2
                self.long_name = "PRECT squared"
                self.units = "(mm/day)^2"
            elif var == "PRECTsqrt":
                self.variable("PRECT", box)
                self.data = self.data**(.5)
                self.long_name = "Square root of PRECT"
                self.units = "sqrt(mm/day)"
            elif var == "SST":
                import numpy as np
                ts = self.variable("TS", box, setData = False)
                ocnfrac = self.variable("OCNFRAC", box, setData = False)
                self.data = ts*np.where(ocnfrac < .9, np.nan, ocnfrac) - 273.15
                self.long_name = "SST"
                self.units = "deg C"
                self.dataDims = ("lat", "lon")
            elif var == "SOIL_d18O":
                soil_h2o = self.variable("H2OSOI", box = box)
                soil_h218o = self.variable("H2OSOI_H218O", box = box)
                self.long_name = "Soil d18O"
                self.units = "d18O"
                self.data = (soil_h218o / soil_h2o - 1) * 1000
            elif var == "SOIL_dD":
                soil_h2o = self.variable("H2OSOI", box = box)
                soil_hdo = self.variable("H2OSOI_HDO", box = box)
                self.long_name = "Soil dD"
                self.units = "dD"
                self.data = (soil_hdo / soil_h2o - 1) * 1000
            elif var == "RUNOFF_d18O":
                runoff_h2o = self.variable("QOVER_H2OTR", box = box)
                runoff_h218o = self.variable("QOVER_H218O", box = box)
                self.long_name = "Runoff d18O"
                self.units = "d18O"
                self.data = (runoff_h218o / runoff_h2o - 1) * 1000
            elif var == "INFILTRATION_d18O":
                infil_h2o = self.mask(self.variable("QINFL_H2OTR", box = box), 'lt', 0)
                infil_h218o = self.mask(infil_h2o, 'lt', 0, self.variable("QINFL_H218O", box = box))
                self.long_name = "Infiltration d18O"
                self.units = "d18O"
                self.data = (infil_h218o / infil_h2o - 1) * 1000
                
            # Regular variables inside the netcdf file
            else:
                try:
                    self.variable(var, box)
                except KeyError:
                    print("From camgoda.ExtractData:\nNot able to extract variable " + var + "...\nSkipping this variable.")
                    print("Is this a 3-spatial-dimension variable? If so, append 3d_ to the beginning of the variable name.")
                    return
            if var_is_3d:
                if self.model == "CAM":
                    self.isobar(pressure)
                elif self.model == "CLM":
                    # This functionality was bolted on and doesnt break old code.
                    # So, although it's not the prettiest way to approach things, there are no bugs.
                    pressure /= 100
                    self.depth(pressure)
            if var_is_ColumnMean:
                self.ColumnMean(box, setData = True)
            if returnData:
                if len(variables) == 1:
                    RETURN = self.data
                else:
                    RETURN[:,:,V_idx] = self.data
        # End of V loop - all data read by this point
        if returnData:
            return RETURN
        else:
            return [var_is_3d, var, pressure]


class popgoda:
    def __init__(self, dataset):
        from netCDF4 import Dataset as ds
        self.dataset = dataset
        if type(dataset) is str:
            self.dataset = ds(dataset, "r")
        self.dims = [d for d in self.dataset.dimensions]
        self.vars = [v for v in self.dataset.variables]
        self.dimlen = [len(self.dataset.dimensions[d]) for d in self.dims]
        
        # Keeping the self.ulat name because I don't feel like going back through the code to change ulat to lat even if the variable is tlat
        try:
            self.ulat = self.dataset.variables['ULAT'][:]
        except KeyError:
            self.ulat = self.dataset.variables['TLAT'][:]
        try:
            self.ulon = self.dataset.variables['ULONG'][:]
        except KeyError:
            self.ulon = self.dataset.variables['TLONG'][:]
    
        self.var = ''
        self.units = ''
        self.long_name = ''
        self.zCalculated = False

        #self.levs = self.dataset.variables['lev'][:]
        self.isTime = any(d == "time" for d in self.dims)
        if self.isTime:
            self.ntime = self.dimlen[self.dims.index("time")]
            if self.ntime == 1:
                self.isTime = False
    def boxMean(self, var, box = None, lev = None, setData = True):
        import numpy as np
        # Extract the data
        DATA = self.dataset.variables[var]
        data = DATA[:].squeeze()
        b,u,l,r = box
        if l < r:
            space_bool = (self.ulat > b) * (self.ulat < u) * (self.ulon < r) * (self.ulon > l)
        else:
            space_bool = (self.ulat > b) * (self.ulat < u) * ~((self.ulon < l) * (self.ulon > r))
        
        idx0, idx1 = np.where(space_bool)
        if lev is not None: # assume 3-d data 
            return np.nanmean(data[lev, idx0, idx1])
        else: # assume 2-d data (e.g. the "selected_*.nc" files)
            return np.nanmean(data[idx0,idx1])
        
        
    def surface(self, var, box = None, setData = True):
        from pygoda import find_indices
        import xarray as xr
        import numpy as np
        from scipy.interpolate import griddata
        # Extract the data
        DATA = self.dataset.variables[var]
        data = DATA[:].squeeze()

        # Add the continent mask with nan values
        if np.ma.is_masked(data):
            data = np.ma.filled(data, fill_value = np.nan)

        ndims = len(data.shape)

        if ndims ==3:
            data = data[0,:,:]
        
        lon = np.arange(0,360,.5)
        lat = np.arange(-90,90,.5)
        ulon = self.ulon.ravel()
        ulat = self.ulat.ravel()
        idxbox = find_indices(box, lat, lon)
        ilat = lat[idxbox[0]]
        ilon = lon[idxbox[1]]

        data = griddata((ulon, ulat), data.ravel(), (ilon[None,:], ilat[:,None]))
        
        if setData:
            self.lat = lat
            self.lon = lon
            self.boxlat = self.lat[idxbox[0]]
            self.boxlon = self.lon[idxbox[1]]
            self.long_name = DATA.long_name
            self.units = DATA.units
            self.data = data

        return data

    def fullvar(self, var, box = None, setData = True):
        #Still a work in progress
        from pygoda import find_indices
        import xarray as xr
        import numpy as np
        from scipy.interpolate import griddata
        # Extract the data
        DATA = self.dataset.variables[var]
        data = DATA[:].squeeze()

        # Add the continent mask with nan values
        if np.ma.is_masked(data):
            data = np.ma.filled(data, fill_value = np.nan)

        ndims = len(data.shape)
        # Regrid the data to a .5 degree rectangular grid
        # This regridding method taken from:
        # https://stackoverflow.com/questions/34408293/2-d-interpolation-ignoring-nan-values 
        # Also box the data out here
        lon = np.arange(0,360,.5)
        lat = np.arange(-90,90,.5)
        ulon = self.ulon.ravel()
        ulat = self.ulat.ravel()
        idxbox = find_indices(box, lat, lon)
        ilat = lat[idxbox[0]]
        ilon = lon[idxbox[1]]

        # Loop over the depth and regrid one layer at a time
        data_regrid = np.zeros(shape = (data.shape[0], len(ilat), len(ilon)))

        for i in range(data_regrid.shape[0]):
            if (float(i) / 60 * 100) % 10 == 0:
                print(var + ' data ' + str(float(i) / 60 * 100) + '% loaded...')
            temp = data[i,:,:].ravel()
            # Re-grid
            data_regrid[i,:,:] = griddata((ulon, ulat), temp, (ilon[None,:], ilat[:,None]))
        
        data = data_regrid
        
        if setData:
            self.lat = lat
            self.lon = lon
            self.boxlat = self.lat[idxbox[0]]
            self.boxlon = self.lon[idxbox[1]]
            self.long_name = DATA.long_name
            self.data = data

        return data

    def isobath(self, depth):
        # Need to finish fullvar first
        import numpy as np
        if not self.zCalculated:
            self.z = self.variable("z_t", setData = False)
            self.zCalculated = "True"

        VAR = np.zeros(shape = (len(self.boxlat), len(self.boxlon)))
        for i in range(len(self.boxlat)):
            for j in range(len(self.boxlon)):
                VAR[i,j] = np.interp(depth, self.z[:], self.data[:,i,j])

        return VAR

    def ExtractData(self, var, box):
        if var == "d18O":
            self.data = (self.surface("R18O", box) - 1) * 1000
            self.long_name = "d18O"
            self.units = "delta-18O"
        else:
            try:
                self.surface(var, box)
            except KeyError:
                print("Not able to plot variable " + var + "...\nSkipping this variable.")
                print("Is this a 3-spatial-dimension variable? If so, append 3d_ to the beginning of the variable name.")
                return

