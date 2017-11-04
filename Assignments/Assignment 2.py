# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt

plt.close('all')


#########################################################

# Question 1

def plikelihood(v, N):

    theta = np.linspace(0, 1, 100)

    p = (theta**v)*((1-theta)**(N-v))

    return p/np.sum(p)


def pprior(a, b):

    theta = np.linspace(0, 1, 100)

    p = (theta**(a-1))*((1-theta)**(b-1))

    return p

def pposterior(v, N, a, b):

    theta = np.linspace(0, 1, 100)

    post = (theta**((v + a)-1))*((1-theta)**((N - v + b)-1))

    norm_post = post/np.sum(post)

    return norm_post

v=0; N=3; a=1; b=1

theta = np.linspace(0, 1, 100)
likelihood = plikelihood(v, N)
prior = pprior(a, b)
posterior = pposterior(v, N, a, b)

print np.mean(likelihood)
print np.std(likelihood)

print np.sum(posterior)

plt.figure()

plt.plot(theta, likelihood)

plt.figure()

plt.plot(theta, prior)

plt.figure()

plt.plot(theta, posterior)

#########################################################

# Question 2



























plt.show()
