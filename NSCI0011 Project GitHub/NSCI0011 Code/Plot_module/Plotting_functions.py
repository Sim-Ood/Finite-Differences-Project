"""This file contains documentation for the plotting functions.

"""
import pydoc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def shape_sol(S,W):
    """This function transposes and flips the solution grids and then removes the padding.

    Args:
        S: Fully updated stream grid - a 2D numpy array.
        W: Fully updated vorticity grid - a 2D numpy array.

    Returns:
        S_sol: The input stream grid with correct indexing and padding removed.
        W_sol: The input vorticity grid with correct indexing and padding removed.
    """
    
    # Correct the indexing so that (i,j) = (x,y)
    S_sol = np.flipud(np.transpose(S))
    W_sol = np.flipud(np.transpose(W))

    # Remove the padding - top row and outermost columns
    S_sol = S_sol[1:,1:-1]
    W_sol = W_sol[1:,1:-1]

    return S_sol, W_sol

def plot_flow(x,y,S_sol,W_sol,beamfront,beamback,beamtop):
    """This function creates contour plots of the stream and vorticity functions.
    The contours are reflected in the centreline and added to the plot to visualise the full flow profile.
    The function superimposes a patch with the beam's dimensions onto the beam region in the contour plots.

    Args:
        x: The array of points along the x-axis.
        y: The array of points along the y-axis.
        S_sol: The solution stream grid, correctly indexed and with padding removed.
        W_sol: The solution vorticity grid, correctly indexed and with padding removed.
        beamfront: The integer index along the x-axis where the beamfront DE is located.
        beamback: The integer index along the x-axis there the beamback FG is located.
        beamtop: The integer index along the y-axis where the beamtop EF is located.

    Returns:
        Contour plot of stream function for the full solution space.
        Contour plot of vorticity function for the full solution space.
    
    """

    # Correct for matplotlib plotting the grid upside down
    S_plot = np.flipud(S_sol)
    W_plot = np.flipud(W_sol)

    # Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))

   
    # Define custom contours
    S_levels = np.array([0.0,0.000001,0.000008,
                         0.00001,0.00002,
                         0.0001,0.0003,0.0009,0.0012,
                         0.01,0.02,0.04,0.08,0.09,0.1,
                         0.2,0.4])
    #S_levels = np.linspace(np.min(S_sol),np.max(S_sol),70) - in case I need this later
    W_levels = np.linspace(np.min(W_sol), np.max(W_sol), 50)


    # Plot stream function contour lines
    cs1 = axs[0].contour(x, y, S_plot,levels=S_levels,colors='black',linewidths=0.9)
    cs2 = axs[0].contour(x, -y, S_plot,levels=S_levels,colors='black',linewidths=0.9)
    fig.colorbar(cs1, ax=axs[0], label='Stream Function')
    axs[0].set_title('Stream Function Contours')
    axs[0].set_xlabel('X')
    axs[0].set_ylabel('Y')

    # Place patch over the beam region in the stream plot
    rect_S = Rectangle((x[beamfront-1], -y[beamtop]), # coordinates of bottom left corner of block - corrected for removed padding
                     x[beamback] - x[beamfront], # width of beam
                     2*y[beamtop], # height of beam - doubled for full solution space
                     alpha=1,facecolor='grey',edgecolor='black',zorder=2)
    axs[0].add_patch(rect_S)

    # Plot vorticity contour lines
    cs3 = axs[1].contour(x, y, W_plot, levels=W_levels,colors='blue',linewidths=0.6,linestyles='solid')
    cs4 = axs[1].contour(x, -y, W_plot, levels=W_levels,colors='blue',linewidths=0.6,linestyles='solid')
    fig.colorbar(cs3, ax=axs[1], label='Vorticity')
    axs[1].set_title('Vorticity Contours')
    axs[1].set_xlabel('X')
    #axs[1].set_ylabel('Y')
    
    # Place patch over the beam region in the vorticity plot
    rect_W = Rectangle((x[beamfront-1], -y[beamtop]), # coordinates of bottom left corner of block - corrected for removed padding
                     x[beamback] - x[beamfront], # width of beam
                     2*y[beamtop], # height of beam - doubled for full solution space
                     alpha=1,facecolor='grey',edgecolor='black',zorder=2)
    axs[1].grid(False)
    axs[1].add_patch(rect_W)

    plt.show()

    return

def plot_errors(S_err,W_err):
    """This function plots the residual errors for the stream and vorticity functions.
    The residual error is the absolute difference in the solution grids between consecutive iterations.
    
    Args:
        S_err: Numpy array of the residuals for each iteration of the stream function grid.
        W_err: Numpy array of the residuals for each iteration of the vorticity function grid.

    Returns:
        Plot of relative errors in the stream function against iterations.
        Plot of relative errors in the vorticity function against iterations.
    
    """
    # Find number of Gauss-Seidel iterations to plot against
    iters = np.arange(1, len(S_err)+1,1)

    # Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))

    
    axs[0].plot(iters, S_err, label='Stream Function')
    axs[0].set_title('Residual Errors in Stream Function')
    axs[0].set_xlabel('Iterations')
    axs[0].set_ylabel('Residual Error')
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(iters, W_err, label='Vorticity Function', color='orange')
    axs[1].set_title('Residual Errors in Vorticity Function')
    axs[1].set_xlabel('Iterations')
    #axs[1].set_ylabel('Residual Error')
    axs[1].legend()
    axs[1].grid(True)

    plt.show()

    return

    

pydoc.writedoc("Plotting_functions")

