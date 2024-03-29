def plot_rigid(ax, x, y, rigid, alpha = 1, cmap = 'binary'):
    '''
    This function returns a plt subplot for drawing the rigid objects on a matplotlib axis
    ax is the axis to plot the rigids on
    x and y are the 1-dimensional x and y coordinates for the grid
    Argument rigid should be the "rigid" variable from obstacles output
    rigid needs to be the correct shape (i.e. 2d and correct axes for plotting)
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    rigid_bool = np.zeros_like(rigid, dtype = 'bool')
    rigid_bool[:] = rigid[:]
    cmap_put = plt.get_cmap(cmap)
    rigid_mask = np.ma.masked_where(~rigid_bool.data, np.ones(rigid_bool.shape))
    return ax.pcolormesh(xf,yf,rigid_mask,shading = 'auto', cmap=cmap_put, alpha = alpha, vmin=0, vmax=1)
