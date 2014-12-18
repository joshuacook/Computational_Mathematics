# HenonMap.py: Plot out iterates of Henon's 2D map
#
# The Henon Map is given by
#     x_n+1 = f( x_n , y_n )
#     y_n+1 = g( x_n , y_n )
# with
#     f( x_n , y_n ) = y_n + 1 - a * x_n^2
#     g( x_n , y_n ) = b * x_n
#
# The state space is R^2
# and the control parameters range in
#       a in [0,2]
#       b in [0,1]
# Henon's original parameters for a chaotic attractor: (a,b) = (1.4,0.3)
# 

def HenonMap(a,b,x,y):
  return y + 1.0 - a *x*x, b * x

# Import plotting routines
from pylab import *

# Simulation parameters
#
# Control parameters:
a = 1.4
b = 0.3
# The number of iterations to throw away
nTransients = 100
# The number of iterations to generate
nIterates = 1000000

# Initial condition
xtemp = 0.1
ytemp = 0.3
for n in xrange(0,nTransients):
  xtemp, ytemp = HenonMap(a,b,xtemp,ytemp)
# Set up arrays of iterates (x_n,y_n) and set the initital condition
x = [xtemp]
y = [ytemp]
# The main loop that generates iterates and stores them
for n in xrange(0,nIterates):
  # at each iteration calculate (x_n+1,y_n+1)
  xtemp, ytemp = HenonMap(a,b,x[n],y[n])
  # and append to lists x and y
  x.append( xtemp )
  y.append( ytemp )

# Setup the plot
xlabel('x(n)') # set x-axis label
ylabel('y(n)') # set y-axis label
title('Henon with (a,b) = (' + str(a) + ',' + str(b) + ')') # set plot title
# Plot the time series
plot(x,y, 'r,')
# Make sure the iterates appear in a unit square
#axis('equal')
axis([-2.0, 2.0, -1.0, 1.0])

# Use command below to save figure
#savefig('HenonMapIterates', dpi=600)

# Display the plot in a window
show()
