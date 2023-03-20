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
    digix = np.digitize(particles_x, edgesx)*np.where(np.isnan(particles_x), np.nan, 1)
    digix = digix[~np.isnan(digix)]
    digiy = np.digitize(particles_y, edgesy)*np.where(np.isnan(particles_y), np.nan, 1)
    digiy = digiy[~np.isnan(digiy)]
    digiz = np.digitize(particles_z, edgesz)*np.where(np.isnan(particles_z), np.nan, 1)
    digiz = digiz[~np.isnan(digiz)]
    
    # Sum up the digitized indices
    for i,j,k in zip(digix,digiy,digiz):
        hold[int(i),int(j),int(k)] += 1
       
    return hold
