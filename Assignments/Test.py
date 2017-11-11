# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc

def norm(x, mu, sig):

    fact = (1./(sig*np.sqrt(2*np.pi)))

    p = fact*(np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))

    return p


X1 = 1.2  #  Initial data

GM1 = np.array([3, 0.2])  #  GM data
GM2 = np.array([1, 0.2])  #  GM data
O = np.array([2.2, 0.75])  #  Organic data


def likelihood_GM(X1, GM1, GM2):

    P_GM1 = norm(X1, GM1[0], GM1[1])
    P_GM2 = norm(X1, GM2[0], GM2[1])

    Tot = 0.5*P_GM1 + 0.5*P_GM2

    return Tot


def likelihood_O(X1, O):

    P_O = norm(X1, O[0], O[1])

    return P_O


prior_GM = 0.7

prior_O = 0.3

P_GM = likelihood_GM(X1, GM1, GM2)*prior_GM

P_O = likelihood_O(X1, O)*prior_O

P_OX = P_O/P_GM

print "Odds of tomato being organic", P_OX

print P_O
print P_GM