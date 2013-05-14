#! /usr/bin/python3.3

import numpy as np
from numpy import *
from numpy import linalg as LA

# 1. Defining the Matrices

class SymNDArray(np.ndarray):
    def __setitem__(self, (i, j), value):
        np.ndarray.__setitem__(self, (i, j), value)
        np.ndarray.__setitem__(self, (j, i), value)

def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())

def symarray(input_array):
    "Returns a symmetrized version of input_array; further assignments to the array are automatically symmetrized."
    return symmetrize(np.array(input_array)).view(SymNDArray)

i=0
j=0

# Example:

A = mat(symarray(np.random.randint(8, size=(8, 8))))

print "\nA is"
print A
v=mat("5.;7.;9.;1.;0.;1.;0.;1.")

print "\nv is"
print v
U=pow(A,100)*v
print "\nU is"
print U

w=U/U[0]
print "\nw is"
print w
print "\nA*w is"
print A*w

e=w.I*A*w
print "\nThe calculated eigenvalue is "
print e

eVal = LA.eig(A)
print eVal
print eVal()[0]


error=(eVal[0]-e)/100
print error