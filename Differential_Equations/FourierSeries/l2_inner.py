import numpy as np
import scipy.linalg.blas as blas
import scipy.integrate as scint
import matplotlib.pyplot as plt
import seaborn as sls
import timeit

n    = 300.
h    = 1./(n-1)
dx=2*np.pi*h
X    = np.linspace(-np.pi,np.pi,n)
ONE  = np.ones(n)
ZERO = np.zeros(n) 

def L2_inner_trapz(f,g,dx):
    return np.trapz(f*g,dx=dx)/np.pi

def L2_inner_dot(f,g,dx):
    return f.dot(g)*dx/np.pi

def L2_inner(f,g,dx):
    return np.inner(f,g)*dx/np.pi

def L2_inner_ddot(f,g,dx):
    return blas.ddot(f,g)*dx/np.pi

def L2_inner_quad(f,g):
    return scint.quad(lambda x: f(x)*g(x),-np.pi,np.pi, epsabs=10e-16)[0]/np.pi  

f = lambda x: np.sin(x)
g = lambda x: np.exp(x)
F = f(X)
G = g(X)


setupstr1  = "import numpy as np;"
setupstr1 += "from __main__ import L2_inner,"
setupstr1 += "L2_inner_dot, L2_inner_ddot, L2_inner_trapz, f, g;"
setupstr1 += "n=%d; X=np.linspace(-np.pi,np.pi,n);"
setupstr1 += "F=f(X); G=g(X);"
setupstr1 += "dx=2*np.pi/(n-1)"

setupstr2  = "import numpy as np;"
setupstr2 += "from __main__ import L2_inner_quad, f, g;"
setupstr2 += "n=%d;"

def time(n):
    setup1 = setupstr1 % n
    setup2 = setupstr2 % n
    time1 = timeit.timeit( 'L2_inner(F,G,dx)', setup1, number=10)
    time2 = timeit.timeit( 'L2_inner_dot(F,G,dx)', setup1, number=10)
    time3 = timeit.timeit( 'L2_inner_ddot(F,G,dx)', setup1, number=10)
    time4 = timeit.timeit( 'L2_inner_trapz(F,G,dx)', setup1, number=10)
    time5 = timeit.timeit( 'L2_inner_quad(f,g)', setup2, number=10)
    return (time1, time2, time3, time4, time5)

## Perform timing for vector product.
times = np.zeros( (7,5) )
for i in range(7):
    times[i,:] = time( 10**(i+1) )

x = 10**np.arange(1,8,1)
f, ax = plt.subplots()
ax.set( xscale='log', yscale='log', title='Inner vs. BLAS vs. trapz', \
        ylabel='time [s]', xlabel='n')
ax.plot( x, times[:,0], label='numpy.inner' )
ax.plot( x, times[:,1], label='numpy.dot')
ax.plot( x, times[:,2], label='blas.ddot')
ax.plot( x, times[:,3], label='numpy.trapz')
ax.plot( x, times[:,4], label='int.quad')
plt.legend()
plt.show()
