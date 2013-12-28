import numpy as np
import matplotlib.pyplot as plt

soa = np.array( [ [0,0,3,2], [0,0,1,1],[0,0,9,9]]) 
X,Y,U,V = zip(*soa)
plt.figure()
axes = plt.gca()
axes.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=1)
axes.set_xlim([-1,10])
axes.set_ylim([-1,10])
plt.draw()
plt.show()