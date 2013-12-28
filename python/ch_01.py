from pylab import *
from numpy import ma

X,Y = meshgrid( arange(0,2*pi,.2),arange(0,2*pi,.2) )
U = X
V = Y

#1
figure()
Q = quiver( U, V)
qk = quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W',
               fontproperties={'weight': 'bold'})
l,r,b,t = axis()
dx, dy = r-l, t-b
axis([l-0.1*dx, r+0.1*dx, b-0.1*dy, t+0.1*dy])

title('Minimal arguments, no kwargs')
show()