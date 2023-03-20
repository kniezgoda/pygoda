import numpy as np
import matplotlib.pyplot as plt

def plotRigid(ax, x, y, rigid, alpha = 1, cmap = plt.cm.binary):
    '''
    This function returns a plt subplot for drawing the rigid objects on a matplotlib axis
    ax is the axis to plot the rigids on
    x and y are the 1-dimensional x and y coordinates for the grid
    Argument rigid should be the "rigid" variable from obstacles output
    rigid needs to be the correct shape (i.e. 2d and correct axes for plotting)
    '''
    rigid_bool = np.zeros_like(rigid, dtype = 'bool')
    rigid_bool[:] = rigid[:]
    rigid_mask = np.ma.masked_where(~rigid_bool.data, np.ones(rigid_bool.shape))
    return ax.pcolormesh(x,y,rigid_mask,shading = 'auto', cmap=cmap, alpha = alpha, vmin=0, vmax=1)


def digi3d(xp,yp,zp,xe,ye,ze):
    # Digitize the data into edge bins
    # Adapted from https://stackoverflow.com/questions/10686847/fast-categorization-binning
    nx = len(xe)+1
    ny = len(ye)+1
    nz = len(ze)+1
    hold = np.zeros(shape = (nx,ny,nz))
    digix = np.digitize(xp, xe)*np.where(np.isnan(xp), np.nan, 1)
    digix = digix[~np.isnan(digix)]
    digiy = np.digitize(yp, ye)*np.where(np.isnan(yp), np.nan, 1)
    digiy = digiy[~np.isnan(digiy)]
    digiz = np.digitize(zp, ze)*np.where(np.isnan(zp), np.nan, 1)
    digiz = digiz[~np.isnan(digiz)]
    
    # Sum up the digitized indices
    for i,j,k in zip(digix,digiy,digiz):
        hold[int(i),int(j),int(k)] += 1
       
    return hold


def digitize(points, edges):
    # Generalized digi3d to n-dimensional data
    # Digitize the data into edge bins
    # Adapted from https://stackoverflow.com/questions/10686847/fast-categorization-binning
    ndim = len(points)
    if ndim > 1:
        if len(points) != len(edges):
            print("points and edges are not same length, exiting")
            return None

    holdshape = []
    for p, e in zip(points, edges):
        holdshape.append(len(e)+1)
        
    hold = np.zeros(shape = holdshape)
    digis = []
    for p, e in zip(points, edges):
        digix = np.digitize(p, e)*np.where(np.isnan(p), np.nan, 1)
        digix = digix[~np.isnan(digix)]
        digis.append(digix)
    
    for i in range(len(digis[0])):
        indexer = []
        for dim in range(3):
            indexer.append(int(digis[dim][i]))
        hold[tuple(indexer)] += 1

    return hold
