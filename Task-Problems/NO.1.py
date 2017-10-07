# -*- coding: utf-8 -*-
"""
Created on Mon Oct 02 16:17:28 2017

@author: Gerwyn Jones
"""

from __future__ import division
import numpy as np
import math


def nCr(n, r):
    
    f = math.factorial
    
    return f(n) / (f(r) * f(n-r))


def P(N, H, T): 
    
    """Probability of Heads H and Tails T from Number of goes N """
        
    Q = 2**N 
    
    t = nCr(N, N-T)
    
    p = nCr(N, H)
    
    if t == p:
        
        c = p
    
    return c/Q


a = P(3, 2, 1)

print a*100
