from pylab import * 

# Set limits and number of points in grid 

xmax = 1.0 
xmin = -xmax 
NX = 20 
ymax = 1.0 
ymin = -ymax 
NY = 20 

# Make grid and calculate vector components 
x = linspace(xmin, xmax, NX) 
y = linspace(ymin, ymax, NY) 
X, Y = meshgrid(x, y) 

U = Y
V = X

figure() 

QP = quiver(X,Y,U,V) 

# quiverkey(QP, 0.85, 1.05, 1.0, '1 mT', labelpos='N') 

# Set the left, right, bottom, top limits of axes 
dx = (xmax - xmin)/(NX - 1) # One less gap than points 
dy = (ymax - ymin)/(NY - 1) 
axis([xmin-dx, xmax+dx, ymin-dy, ymax+dy]) 
xlabel('x') 
ylabel('y') 
show() 
