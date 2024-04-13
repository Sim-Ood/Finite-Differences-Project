# Finite-Differences-Project
This repository contains the Python files, html documentation, and report for the final project I submitted for the Topics in Scientific Computing course at UCL in 2024.

The following is a brief overview of the [project report](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Project%20GitHub/Finite%20Difference%20Method%20Project%20Final%20Report.pdf). 

## Methodology

### Problem Definition
Our project aim is to solve and visualise the two-dimensional flow profile of an oncoming stream with constant speed U around a submerged beam. The fluid behaviour is governed by two coupled partial differential equations, the stream function $\psi$ and vorticity $\omega$.

$$-\omega = \frac{\partial^{2}{\psi}}{\partial{x^2}}+\frac{\partial^{2}{\psi}}{\partial{y^2}}$$
$$\frac{\partial\psi}{\partial{y}}\frac{\partial\omega}{\partial{x}}-\frac{\partial\psi}{\partial{x}}\frac{\partial\omega}{\partial{y}}= \nu\left(\frac{\partial^{2}{\omega}}{\partial{x^2}} + \frac{\partial^{2}\omega}{\partial{y}^2}\right)$$

We discretise and solve the governing equations numerically to obtain steady solutions. From these, we can map contours showing the fluid particles' local rotation (vorticity) and paths of flow (streamlines) in the two-dimensional solution space to construct the flow profile. Our key challenge is implementing boundary conditions, derived from considering the physics of the system. Once solutions are obtained, the variation of the flow profile with its key parameters, such as the Reynolds number, is investigated. 

### Problem Geometry
Because the flow is symmetric about the centreline of the flow profile, computation is only performed over half the solution space which is later mirrored in the centreline to recover the full solution. We model the problem as a channel through which the fluid flows, entering through the inlet and leaving through the outlet.

![Diagram of the problem geometry](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Images/Half%20Flow%20Profile.pdf)

<div align="center">
    <img height="300" src="NSCI0011 Images/Half Flow Profile.pdf">
</div>


We parameterise the computational model with the dimensionless Reynolds number $R$, to support its scalability for differently sized problems.

$$R = \frac{UL}{\nu}\quad\\
\text{where U = background flow velocity,}\quad\nu = \text{viscosity,}\quad\text{L = characteristic length}$$

|           |                                                                                  |                               |
|-----------|----------------------------------------------------------------------------------|-------------------------------|
| Parameter | Description                                                                      | Value                         |
| $h$       | The uniform grid spacing between points,equal in both the $x$ and $y$ direction. | variable and dependent on $n$ |
| $n$       | The number of discrete points on the $y$-axis, determines grid resolution.       | variable                      |
| $L$       | The $y$-axis spanning 0-1, which is the characteristic length.                   | 1 $n$ grid points.            |
| $M$       | The $x$-axis spanning 0-2.                                                       | 2 $n$ grid points.            |
| $t$       | Height of the beam.                                                              | 0.14 $n$ grid points          |
| $w$       | Width of the beam.                                                               | 0.08 $n$ grid points          |
| $AD$      | Length between origin and beamfront $DE$, determining the beam placement.        | 0.15 $n$ grid points          |
| $U$       | Background flow velocity of the stream.                                          | 1                             |
| $R$       | The Reynolds number.                                                             | variable                      |


### Numerical Method

Due to the finite precision of computers, we cannot treat derivative terms in $\psi$ and $\omega$ as the rate of change over infinitesimally small intervals. Instead, we discretise the spatial domain of the problem into a grid of equally spaced points, and approximate derivates as the average rate of change over a finite region. We choose the central difference method which calculates the derivative at a point $y_k$ as the slope between the points either side, $y_{k-1}$ and $y_{k+1}$. For our two-dimensional functions of form $U_{i,j}$, this method allows each point to be expressed in terms of the 4 nearest neighbouring points and the grid spacing $h$.

![Diagram of the two-dimensional central difference method](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Images/5-point%20Method.pdf)

$$U_{i,j}=\frac{1}{4}\left[U_{i+1,j} + U_{i-1,j}+U_{i,j+1}+U_{i,j-1}-h^{2}f(x_{i},y_{j})\right]\\
\text{where}\quad\nabla^{2}U = f(x_i,y_j) = \frac{\partial^{2}}{\partial{x^2}}+\frac{\partial^{2}}{\partial{y^2}}$$

We discretise the governing equations with the central difference method, obtaining update rules where $R_G$ is the Reynolds number scaled to the size and resolution of the solution grids.

$$ \psi_{i,j}=\frac{1}{4}\left(\psi_{i+1,j} + \psi_{i-1,j}+\psi_{i,j+1}+\psi_{i,j-1}+h^{2}\omega_{i,j}\right)$$

$$\omega_{i,j}=\frac{1}{4}\left(\omega_{i+1,j} + \omega_{i-1,j}+\omega_{i,j+1}+\omega_{i,j-1}\right)
\\-\frac{R_G}{16}\left[\left(\psi_{i,j+1}-\psi_{i,j-1}\right)\left(\omega_{i+1,j}-\omega_{i-1,j}\right)-\left(\psi_{i+1,j}-\psi_{i-1,j}\right)\left(\omega_{i,j+1}-\omega_{i,j-1}\right)\right]$$

Separate solution grids are constructed for $\psi$ and $\omega$, which we populate with an initial guess of values. We employ Gauss-Seidel iteration to continually loop over the grids, updating the value at each grid point until the solutions converge. We exclude points along the grid boundaries from the loop as their neighbouring points fall outside the solution domain and would raise an error. To address this, we introduce "ghost" points around the grid boundaries, which allow the boundary points to be updated. The ghost points are removed after the solutions converge to restore the original domain.

![Diagram of Gauss-Seidel iteration updating a grid of $M$ rows and $N$ columns with a padding of ghost points (shown in orange)](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Images/Gauss%20Seidel.pdf)

For these equations, the indexing $(i,j)$ represents (columns,rows) with the origin $(0,0)$ at the bottom left corner of the grid. In Python however, $(i,j)$ represents (rows,columns) with the origin at the top left corner. For simplicity, we compute solutions using the equations with their original indexing and correct afterwards by transposing and flipping the grids.

![Diagram of axes transformation performed on the solution grids](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Images/Axes.pdf)


### Boundary Conditions 

Between loops, we implement boundary conditions to ensure the problem geometry is not overwritten by the update rules. We determine the boundary conditions by considering the physics of the system, and discretise as before. Please see the [project report](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Project%20GitHub/Finite%20Difference%20Method%20Project%20Final%20Report.pdf) for the full derivation.

|               |                                    |                                                                          |
|:-------------:|:----------------------------------:|:------------------------------------------------------------------------:|
|               |                                    |                                                                          |
| Boundary      | $\psi$                             | $\omega$                                                                 |
| Inlet AB      | $\psi_{i+1,j} = \psi_{i-1,j}$      | $\omega = 0$                                                             |
| Outlet CH     | $\psi_{i+1,j} = \psi_{i-1,j}$      | $\omega_{i+1,j} = \omega_{i-1,j}$                                        |
| Surface BC    | $\psi_{i,j+1} = \psi_{i,j-1} + 2h$ | $\omega = 0$                                                             |
| Centreline AH | $\psi  = 0$                        | $\omega = 0$                                                             |
| Beamfront DE  | $\psi = 0$                         | $\omega = \text{-}\frac{\psi_{i+1,j} - 2\psi_{i,j} + \psi_{i-1,j}}{h^2}$ |
| Beamback FG   | $\psi = 0$                         | $\omega = \text{-}\frac{\psi_{i+1,j} - 2\psi_{i,j} + \psi_{i-1,j}}{h^2}$ |
| Beamtop EF    | $\psi = 0$                         | $\omega = \text{-}\frac{\psi_{i,j+1}-2\psi_{i,j}+ \psi_{i,j-1}}{h^2}$    |
|               |                                    |                                                                          |

### Convergence

We confirm convergence once the relative error in the solutions has dropped to $\leqslant0.018%$. The relative error is found by evaluating the difference in the solution grids between iterations.

## Results


We find the flow profile for multiple Reynolds numbers. Below is a sample of the results in the [project report](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Project%20GitHub/Finite%20Difference%20Method%20Project%20Final%20Report.pdf) for $R = 1000$ and $R = 15000$.

![Flow Profile for R = 1000](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Images/R_1000.png)

![Flow Profile for R = 15000](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Images/R_15000.png)






















