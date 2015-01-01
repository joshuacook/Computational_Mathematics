
# Simple Power Method

1. Choose a starting vector $\mathbf{x}^{(0)}\in\mathcal{R}^n$ with
$||\mathbf{x}^{(0)}||=1$.
2. $k=0$
3. `while` some convergence criteria is not satisfied

    i. $k:=k+1$
    ii. $\mathbf{y}^{(k)}:=A\mathbf{x}^{(k-1)}$
    iii. $\mu_k:=||\mathbf{y}^{(k)}||$
    iv. $\mathbf{x}^{(k)}:=\mathbf{y}^{(k)}/\mu_k$


    import numpy.linalg, numpy.random, numpy as np, math
    from random import random as rand
    from numpy.linalg import eig


    B = numpy.array([[2,-12],[1,-5]])
    y = numpy.array([1,1])
    x = y
    for i in range(100):
        y = B.dot(x)
        mu = math.sqrt(y.dot(y))
        x = y/mu
    print numpy.transpose(eig(B)[1])
    print x

    [[ 0.9701425   0.24253563]
     [ 0.9486833   0.31622777]]
    [ 0.9486833   0.31622777]



    A = numpy.random.rand(3,3)
    x = numpy.random.rand(3)
    for i in range(20):
        y = A.dot(x)
        mu = math.sqrt(y.dot(y))
        x = y/mu
    np.transpose(eig(A)[1]),x




    (array([[ 0.51063895,  0.65312477,  0.55917429],
            [ 0.87551464, -0.17084359, -0.45198073],
            [-0.00620138, -0.53208609,  0.84666755]]),
     array([ 0.51063895,  0.65312477,  0.55917429]))


