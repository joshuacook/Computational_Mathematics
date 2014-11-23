from __future__ import division
import matplotlib
matplotlib.use('MacOSX') ## Change if you are using a different backend
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

def SolveLaplace(nx, ny, dx, epsilon = 1e-5, imax = 1000):

    ## Initialize the mesh with some values
    T = np.zeros((nx+1, ny+1))

    ## Set boundary conditions for the problem
    T[0,:] = 1 ## Top Boundary
    T[nx,:] = 0 ## Bottom Boundary
    T[:,0] = 1   ## Right Boundary
    T[:,ny] = 0 ## Left Boundary

    ## Store previous grid values to check against error tolerance
    TN = T + np.zeros((nx+1, ny+1))
    err = TN - T

    ## Constants
    k = 1          ## Iteration counter

    ## Iterative procedure
    while k <= imax:

        for i in np.arange(1., nx):

            for j in np.arange(1.,  ny):

                TN[i,j] = (T[i-1,j] + T[i+1,j] + T[i,j-1] + T[i,j+1])/4.
                err[i,j] = np.abs(TN[i,j] - T[i,j])

        T = TN + np.zeros((nx+1, ny+1))
        k += 1
        errmax = np.max(np.max(err))

        if errmax < epsilon:

            print("Convergence after ", k, " iterations.")
            return T

    print("No convergence after ", k, " iterations.")
    return False

def PlotSolution(nx,ny,dx,T):

    ## Set up x and y vectors for meshgrid
    x = np.linspace(0, nx * dx, nx+1)
    y = np.linspace(0, ny * dx, ny+1)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X, Y, T.transpose(), rstride=1, cstride=1, cmap=cm.cool, linewidth=0, antialiased=False)
    plt.xlabel("X")
    plt.ylabel("Y")
    #plt.zlabel("T(X,Y)")

    fig2 = plt.figure()
    cs = plt.contourf(X, Y, T.transpose(), 32, rstride=1, cstride=1, cmap=cm.cool)
    plt.colorbar()
    plt.xlabel("X")
    plt.ylabel("Y")

    plt.show()

## Size of plate and mesh
nx = 32
ny = 32
dx = 1/32

epsilon = 1e-2 ## Absolute Error tolerance
imax = 10000    ## Maximum number of iterations allowed

T = SolveLaplace(nx, ny, dx, epsilon, imax)
PlotSolution(nx, ny, dx, T)