'''
Created on Sep 29, 2016

@author: Keisuke
'''
from __future__ import division
from pylab import *
from numpy import *

class Solver_3D:
    def __init__(self,X):
        self.X=X

def solver_3D(a,b,c,d):
    b3a = b/(3*a)
    ca = c/a
    b3a2 = b3a*b3a
    Q = b3a2 - ca/3
    R = b3a2*b3a - b3a*ca/2 + d/(2*a)
    th = arccos(R/sqrt(Q*Q*Q)) + 2*pi
    x2 = -2*sqrt(Q)*cos(th/3) - b3a
    X1 = Solver_3D(real(x2))
    return(X1)



