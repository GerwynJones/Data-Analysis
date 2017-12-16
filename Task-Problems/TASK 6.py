# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt


P = np.array([79, 82, 85, 88, 90])


T = np.array([8, 17, 30, 37, 52])


n = len(T)

N = 5000

Tmin = -500

Tmax = 0

sigma_T = 2.0


A = np.linspace(Tmin, Tmax, N)

B = np.linspace(3, 4, N)


X_2 = np.zeros((len(A), len(B)))


def X_sum(A, B):

    Num = (T - A - B*P)**2.0

    Den = sigma_T**2.0

    Sum = np.sum(Num/Den)

    return Sum


for i in range(len(A)):

    for j in range(len(B)):

        X_2[i, j] = X_sum(A[i], B[j])


for i in range(len(A)):

    found = False

    for j in range(len(B)):

        if X_2[i, j] == np.min(X_2):

            print A[i]

            print B[j]

            found = True

            break

        else:

            continue

    if found:

        break


Analytical_A = ((np.sum(P**2.0)*np.sum(T)) - (np.sum(P)*np.sum(P*T)))/((n*np.sum(P**2.0)) - (np.sum(P))**2.0)

Analytical_B = ((n*np.sum(T*P)) - (np.sum(P)*np.sum(T)))/((n*np.sum(P**2.0)) - (np.sum(P))**2.0)


print Analytical_A
print Analytical_B


plt.figure()

plt.plot(P, T, '.')

plt.ylabel('Temperature (C)')
plt.xlabel('Pressure (mm)')


plt.show()

