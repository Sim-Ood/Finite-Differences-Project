"""This file contains documentation for the initialising functions.

"""
import pydoc
import numpy as np

def init_grid(n,Rey):

    """This function initialises the discretised solution grids S and W with padding.
    
    S has ghost points and requires padding along the Inlet, Outlet, and Surface.
    W has ghost points and requires padding along the Outlet.
    
    S and W will be updated in the same loop, so the same padding is added to both to avoid indexing issues.
    Hence, both S and W must be padded on the top and outer edges, i.e. with 1 row and 2 columns.
    Because S and W will be transposed and flipped later, this function adds 2 rows and 1 column of padding.
    
    Args:
        n: Number of divisions in the y-direction between 0-1 inclusively.
        Rey: The Reynolds number.

    Returns:
        S_zero: Initialised padded stream grid of zeros - a 2D numpy array.
        W_zero: Initialised padded vorticity grid of zeros - a 2D numpy array.
        R: The Grid Reynolds number.
        x: The array of points along the x-axis.
        y: The array of points along the y-axis.
        h: The unit of equal grid spacing.

    """
    # Initialise axes and constants

    y = np.linspace(0,1,n)
    h = y[1] - y[0]
    x = np.linspace(0,2+h,2*n) # (2+h) ensures grid spacing in x = grid spacing in y
    R = Rey*h # We have set characteristic length L = y = 1

    # Initialise padded grids

    S_zero = np.zeros((len(x)+2,len(y)+1)) # Stream function grid
    W_zero = np.zeros((len(x)+2,len(y)+1)) # Vorticity function grid

    return S_zero, W_zero, R, x, y, h




def init_ghost(S_zero,W_zero,val_S,val_W):

    """This function sets the starting values of the ghost points in the initialised grids.
    These are the ghost points that will updated by the interior grid points during iteration.
    This allows more control over the initial conditions which effect how quickly the solutions converge. 
    
    Args:
        S_zero: Initialised padded stream grid of zeros - a 2D numpy array.
        W_zero: Initialised padded vorticity grid of zeros - a 2D numpy array.
        val_S: The numerical starting value the ghost points in the stream S grid will be set to.
        val_W: The numerical starting value the ghost points in the vorticity W grid will be set to.

    Returns:
        S_ghosted: The input stream grid with its ghost points set to val - a 2D numpy array.
        W_ghosted: The input vorticity grid with its ghost points set to val - a 2D numpy array.

    """
    # Inlet Ghost Points

    S_zero[0,:] = val_S
    # W inlet is a padding of 0s which are not ghost points

    # Outlet Ghost Points
    S_zero[-1,:] = val_S
    W_zero[-1,:] = val_W

    # Surface Ghost Points
    S_zero[:,-1] = val_S
    # W surface is a padding of 0s which are not ghost points

    S_ghosted = S_zero
    W_ghosted = W_zero

    return S_ghosted,W_ghosted

def initialise(n,Rey,val_S,val_W):
    """This function calls init_grid and init_ghost to initialise the solution grids and key constants.
    The initialised grids are padded and their ghost points are set to the starting values.
    This function also returns the grid Reynolds number R, the axis arrays x and y, and the grid spacing h.

    Args:
        n: Number of divisions in the y-direction between 0-1 inclusively.
        Rey: The Reynolds number.
        val_S: The numerical starting value the ghost points in the stream S grid will be set to.
        val_W: The numerical starting value the ghost points in the vorticity W grid will be set to.

    Returns:
        S: Fully initialised stream function grid - a 2D numpy array.
        W: Fully intialised vorticity function grid - a 2D numpy array.
        R: The Grid Reynolds number.
        x: The array of points along the x-axis.
        y: The array of points along the y-axis.
        h: The unit of equal grid spacing.
    
    """

    # Create padded grid of zeros and key constants
    S_zero,W_zero,R,x,y,h = init_grid(n,Rey)

    # Set the starting values of the ghost points
    S,W = init_ghost(S_zero,W_zero,val_S,val_W)

    return S,W,R,x,y,h

pydoc.writedoc("Initialising_functions")


 

