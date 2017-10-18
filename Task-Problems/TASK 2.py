# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt


""" A health centre is introducing a new treatment for chronic back pain, which in clinical trials was able to cure 
symptoms in 70% of sufferers. The centre has a group 66 people with chronic back pain who have agreed to 
undergo the procedure. """

plt.close('all')

n = 66

p = 0.7

N = 100000

X = np.linspace(0, n, n + 1)

Y = sc.binom.pmf(X, n, p)

plt.plot(X, Y)


M_T = sc.binom.mean(n, p, loc=0)

S_T = sc.binom.std(n, p, loc=0)

Po50 = np.sum(Y[X>=50])

print Po50


plt.figure()


Data = np.random.rand(n, N)

TF = (Data < 0.7)*1.0

Result = np.sum(TF, axis=0)

M_F = np.mean(Result)

S_F = np.std(Result)

plt.hist(Result, bins=50)

print(M_T)
print(S_T)

print(M_F)
print(S_F)


plt.show()