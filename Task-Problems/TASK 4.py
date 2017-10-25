# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt
from scipy.special import factorial


# Question 1


N = 14


P_new = (12/N)

P_old = (2/N)


P = 0.5  # makes no difference


def B(N, p, v):

    T = np.linspace(0, N, N+1)

    fact = factorial(N)/(factorial(T)*factorial(N-T))

    B = fact*(p**T)*((1-p)**(N-T))

    Bino = np.sum(B[T >= v])

    return Bino

print B(N, P, 12)*100
