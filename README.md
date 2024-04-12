# Finite-Differences-Project
This repository contains the Python files, html documentation, and report for my final project, submitted for the Topics in Scientific Computing course at UCL in 2024.

The following is a summary of this project, covered in detail in the [project report](https://github.com/Sim-Ood/Finite-Differences-Project/blob/main/NSCI0011%20Project%20GitHub/Finite%20Difference%20Method%20Project%20Final%20Report.pdf).

## Problem Definition
The project aim is to solve and visualise the two-dimensional flow profile of an oncoming stream with constant speed U around a submerged beam. The fluid behaviour is governed by two coupled partial differential equations, the stream function $\psi$ and vorticity $\omega$ are described by velocity components $\underline{u} = (u,v)$ as shown below.

  $$\frac{\partial{\psi}}{\partial{y}} = u,\quad \frac{\partial{\psi}}{\partial{x}} = -v,\quad \omega = \frac{\partial{v}}{\partial{x}}-\frac{\partial{u}}{\partial{y}}$$

Finite difference methods are employed to discretise and numerically solve these equations to find steady solutions that satisfy $\frac{\partial{\omega}}{\partial{t}} = 0$. From this, contours showing the fluid particles' local rotation (vorticity) and paths of flow (streamlines) are mapped in the two-dimensional solution space to construct the flow profile. The key challenge of this project is implementing boundary conditions, derived from considering the physics of the system. Once solutions are obtained, the variation of the flow profile with its key parameters, such as the Reynolds number, is also investigated. 

## Problem Geometry
Because the flow is symmetric about the centreline of the flow profile, computation is only performed over half the solution space which is later mirrored in the centreline to recover the full solution. The solution space is modelled as a channel through which the fluid flows, entering through the inlet and leaving through the outlet. The computational model is parameterised by the Reynolds number $R$, to support its scalability for differently sized problems.
$$R = \frac{UL}{\nu},\quad \text{where U = background flow velocity, \nu = viscosity, L = characteristic length}$$

## Numerical Method
To solve differential equations computationally, the gradient at a point $(x_k,y_k)$ on some function $y(x_k)$ is found from the slope over a finite interval. The central difference method finds the gradient at $(x_k,y_k)$ using the neighbouring points either side, i.e. over the interval $y_{k-1}$ to $y_{k+1}$. Similarly, the second order differential is obtained from considering the change of gradient over the finite interval. For a two-dimensional problem, the point $U_{i,j}$ in a solution grid for $U(x_i,y_j)$ is found using the 4 nearest neighbouring points.

$$\text{1st Order:}\quad\frac{dy}{dx}\approx\frac{\Delta{y}}{\Delta{x}} =\frac{y_{k+1}-y_{k+1}}{2h},\quad
\quad\text{2nd order:}\quad\frac{d^2y}{dx^2}\approx
    \frac{\frac{y_{k+1}-y_k}{h}-\frac{y_k-y_{k-1}}{h}}{h} = \frac{y_{k+1}-2y_k + y_{k-1}}{h^2}$$

  These formulae are applied to the governing equations to obtain the discretised update rules.


  

  
    










