"""This is the main python module where we call all the functions to initialise,
iteratively update, and then plot the stream S and vorticity W grids.
"""

#Please note the function files must be in the same folder as this Results file for the imports to run.
#They will need to be taken out of the module folders e.g. Upd_module and put into one folder.

#Alternatively, you can try ' from Upd_module.Updating_functions import...' to retrieve from folders.
#But this may depend on how the files are organised.


from Initialising_functions import initialise
from Updating_functions import get_beam, apply_update_rules, apply_boundary_conditions
from Plotting_functions import shape_sol, plot_flow, plot_errors
import numpy as np
import pydoc

def main(n, Rey, val_S,val_W, max_sweeps,tol, start_beam_at, prop_width, prop_height):
    """This is the main function to execute the full process and call the other functions.
    
    Args:
        n: Number of divisions in the y-direction between 0-1 inclusively.
        Rey: The Reynolds number.
        val_S: The numerical starting value the ghost points in the stream S grid will be set to.
        val_W: The numerical starting value the ghost points in the vorticity W grid will be set to.
        max_sweeps: The maximum number of Gauss-Seidel iterations the function will perform.
        tol: The tolerance for convergence, if relative error is less than or equal to tol, the iteration loop breaks.
        start_beam_at: The proportion along the x-axis where the beamfront DE is located - a float between 0-1.
        prop_width: Width of the beam as a proportion of the x-axis - a float between 0-1.
        prop_height: Height of the beam as a proportion of the y-axis - a float between 0-1.
        
    Results:
        Contour plot of stream function for the full solution space.
        Contour plot of vorticity function for the full solution space.
        Plot of relative errors in the stream function against iterations.
        Plot of relative errors in the vorticity function against iterations.

    """
    
    # Initialise grid, ghost points, and constants
    S,W,R,x,y,h = initialise(n, Rey, val_S,val_W)
    

    # Initialise empty grids to store the solution history
    S_history = np.empty_like(S)
    W_history = np.empty_like(W)

    # Initialise empty error arrays
    S_err = []
    W_err = []

    # Define beam placement
    beamfront,beamback,beamtop = get_beam(start_beam_at, prop_width, prop_height, x, y)

    # Update grid with Gauss-Seidel iteration
    for k in range(max_sweeps):
        
        # Take snapshots
        S_history[:] = S[:] # To avoid pointer issues
        W_history[:] = W[:] # To avoid pointer issues

        # Update Grids
        S, W = apply_update_rules(S, W, R, n, h)
        S, W = apply_boundary_conditions(S, W, h, beamfront, beamback, beamtop)

        # Find residual errors for plotting
        S_residual_err = np.linalg.norm(S-S_history)/np.linalg.norm(S)
        W_residual_err = np.linalg.norm(W-W_history)/np.linalg.norm(W)
        S_err.append(S_residual_err)
        W_err.append(W_residual_err)

        if S_residual_err <= tol and W_residual_err <= tol:
            print('Number of Sweeps to Convergence =', k)
            break
  
   
    # Transpose, flip and remove padding
    S_sol, W_sol = shape_sol(S,W)
    
    # Plot contours for stream and vorticity functions
    plot_flow(x, y, S_sol, W_sol, beamfront, beamback, beamtop)

    # Plot residual errors for stream and vorticity functions
    plot_errors(S_err, W_err)
    


if __name__ == "__main__":

    # Define constants
    n = 60
    Rey = 1000
    val_S = 0.5
    val_W = -1
    max_sweeps = 6000
    tol = 18e-5
    start_beam_at = 0.15
    prop_width = 0.08
    prop_height = 0.14

   

    # Call main function with arguments
    main(n, Rey, val_S, val_W, max_sweeps, tol, start_beam_at, prop_width, prop_height)

   

#pydoc.writedoc("Results_File")

