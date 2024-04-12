"""This file contains documentation for the updating functions.

These functions are written with the assumption that the solution grids will be transposed and flipped.
Hence for these functions, the (i,j) indexing of the solution grids should be treated as (x,y).

"""
import pydoc
import numpy as np

def apply_update_rules(S,W,R,n,h):
    """The function loops over the solution grids, applying the update rule to each point in the loop range.
    This excludes the boundaries (the top row, bottom row, and outermost columns).

    Args:
        S: Stream grid - a 2D numpy array.
        W: Vorticity grid - a 2D numpy array.
        R: Grid Reynolds number.
        n: Number of divisions in the y-direction between 0-1 inclusively.
        h: The unit of equal grid spacing.

    Returns:
        S: The input stream grid updated by one Gauss-Seidel iteration.
        W: The input vorticity grid updated by one Gauss-Seidel iteration.

    """

    # Find the stopping index for the loop range, S and W are the same shape
    stop_index_i,stop_index_j = S.shape

    for i in range(1,stop_index_i - 1):
        for j in range(1,stop_index_j - 1):
            S[i,j] = (1/4) * (S[i+1,j] + S[i-1,j] + S[i,j+1] + S[i,j-1] + ((h**2)*W[i,j]))
            W[i,j] = (1/4) * (W[i+1,j]+W[i-1,j]+W[i,j+1]+W[i,j-1]) - (
                (R/16)*(((S[i,j+1]-S[i,j-1])*(W[i+1,j]-W[i-1,j]))-((S[i+1,j]-S[i-1,j])*(W[i,j+1]-W[i,j-1])))
            )                                               
            
    return S,W



def grid_bound(S,W,h):
    """This function applies the boundary conditions to the solution grids.
    This updates the ghost points using the interior grid points and imposes fixed conditions.

    Args:
        S: Stream grid - a 2D numpy array.
        W: Vorticity grid - a 2D numpy array.
        h: The unit of equal grid spacing.
        U: The background flow velocity.

    Returns:
        S: The input stream grid with its grid boundaries updated.
        W: The input vorticity grid with its grid boundaries updated.

    """
                
    # Inlet Conditions (AB)
    S[0,:] = S[2,:]
    W[1,:] = 0

    # Outlet Conditions (CH)
    S[-1,:] = S[-3,:]
    W[-1,:] = W[-3,:]

    # Surface Conditions (BC)
    S[:,-1] = S[:,-3] + 2*h
    W[:,-2] = 0

    # Centreline Conditions (AH)
    S[:,0] = 0
    W[:,0] = 0

    return S,W

def get_beam(start_beam_at,prop_width,prop_height,x,y):
    """This function defines the index values where the beam is located in the solution grids.
    
    Args:
        start_beam_at: The proportion along the x-axis where the beamfront DE is located - a float between 0-1.
        prop_width: Width of the beam as a proportion of the x-axis - a float between 0-1.
        prop_height: Height of the beam as a proportion of the y-axis - a float between 0-1.
        x: The array of points along the x-axis.
        y: The array of points along the y-axis.
    
    Returns:
        beamfront: The integer index along the x-axis where the beamfront DE is located.
        beamback: The integer index along the x-axis there the beamback FG is located.
        beamtop: The integer index along the y-axis where the beamtop EF is located.

    """
    # Get axis points
    num_points_x = len(x)
    num_points_y = len(y)

    # Find beam indices - these must be integers.
    beamfront = (np.rint((start_beam_at*num_points_x))).astype(int)
    beamback = (beamfront + np.rint((prop_width*num_points_x))).astype(int)
    beamtop = (np.rint(prop_height*num_points_y)).astype(int)

    return beamfront,beamback,beamtop


def beam_bound(S,W,beamfront,beamback,beamtop,h):
    """This function applies the boundary conditions at the beam surfaces.

    Args:
        S: Stream grid - a 2D numpy array.
        W: Vorticity grid - a 2D numpy array.
        beamfront: The integer index along the x-axis where the beamfront DE is located.
        beamback: The integer index along the x-axis there the beamback FG is located.
        beamtop: The integer index along the y-axis where the beamtop EF is located.
        h: The unit of equal grid spacing.

    Returns:
        S: The input stream grid with its beam boundaries updated.
        W: The input vorticity grid with its beam boundaries updated.

    """

    # Beamfront conditions (DE)
    S[beamfront,0:beamtop+1] = 0
    W[beamfront,0:beamtop+1] = -(S[beamfront+1,0:beamtop+1]-(2*S[beamfront,0:beamtop+1])+S[beamfront-1,0:beamtop+1])/h**2

    # Beamback conditions (FG)
    S[beamback,0:beamtop+1] = 0
    W[beamback,0:beamtop+1] = -(S[beamback+1,0:beamtop+1]-(2*S[beamback,0:beamtop+1])+S[beamback-1,0:beamtop+1])/h**2

    # Beamtop conditions (Top EF)
    S[beamfront:beamback+1,beamtop] = 0
    W[beamfront:beamback+1,beamtop] = -(S[beamfront:beamback+1,beamtop+1]-(2*S[beamfront:beamback+1,beamtop])
                                        +S[beamfront:beamback+1,beamtop-1])/h**2
    
    return S,W


def apply_boundary_conditions(S,W,h,beamfront,beamback,beamtop):
    """This function calls grid_bound and beam_bound to apply the full set of boundary conditions.
    This implements the conditions at the Inlet, Outlet, Surface, Centreline, Beamfront, Beamback, and Beamtop.

    Args:
        S: Stream grid - a 2D numpy array.
        W: Vorticity grid - a 2D numpy array.
        h: The unit of equal grid spacing.
        beamfront: The integer index along the x-axis where the beamfront DE is located.
        beamback: The integer index along the x-axis there the beamback FG is located.
        beamtop: The integer index along the y-axis where the beamtop EF is located.

    Retuns:
        S: The input stream grid with all boundary conditions applied.
        W: The input vorticity grid with all boundary conditions applied.

    """
    # Enforce conditions on grid boundaries
    S_gridbound, W_gridbound = grid_bound(S,W,h)

    # Enforce conditions on beam surfaces
    S,W = beam_bound(S_gridbound,W_gridbound,beamfront,beamback,beamtop,h)

    return S,W


pydoc.writedoc("Updating_functions")
