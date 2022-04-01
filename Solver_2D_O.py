'''
Created on Sep 29, 2016

@author: Keisuke
'''
from __future__ import division
from pylab import *

class DSolver:
    def __init__(self, a, b, c):
        self.a = a
        self.b3a = b/(3*a)
        ca = c/a
        b3a2 = self.b3a*self.b3a
        self.R0 = b3a2*self.b3a - self.b3a*ca/2
        Q = b3a2 - ca/3
        self.rQ = sqrt(Q)
        self.rQ3 = sqrt(Q*Q*Q)
        
    def __call__(self, d):
        R = self.R0 + .5*d/self.a
        th = arccos(R/self.rQ3) + 2*pi
        x2 = -2*self.rQ*cos(th/3) - self.b3a
        return(x2)



