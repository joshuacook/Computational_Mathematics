from pylab import * 
import os

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
R = sqrt(X**2 + Y**2)
U = X/R
V = Y/R
 

fig = figure() 

QP = quiver(X,Y,R**2,R**2) 

# quiverkey(QP, 0.85, 1.05, 1.0, '1 mT', labelpos='N') 

# Set the left, right, bottom, top limits of axes 
dx = (xmax - xmin)/(NX - 1) # One less gap than points 
dy = (ymax - ymin)/(NY - 1) 
axis([xmin-dx, xmax+dx, ymin-dy, ymax+dy]) 
xlabel('x') 
ylabel('y') 
# title('$\\frac{x\mathbf{i}+\mathbf{j}y}{\mathbf{r}}$', fontsize=24)
axhline(0, color='black')
axvline(0, color='black')

fig.savefig('dgc_ch01_ex_1_3b.png')
os.system('mv dgc_ch01_ex_1_3b.png ~/Dropbox/img/')
# show() 
