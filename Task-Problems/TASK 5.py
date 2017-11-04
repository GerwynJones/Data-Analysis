# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt
import random


np.random.seed(70)


# 7 island MCMC


Year = 365

t = 40

Walk = np.zeros(t*Year)

count = 1

island_no = np.array([1., 2., 3., 4., 5., 6., 7.])

N = np.zeros(len(island_no))

for j in island_no:

    N[np.int(j) - 1] = j*1000

decision = raw_input("Starting island?: ")

i = np.int(decision) - 1

while count < t*Year:

    Walk[0] = i + 1

    pspin = np.random.uniform(0, 1, 1)

    Coin = np.random.random_integers(0, 1)

    if Coin == 0:
        # head west

        if i == 0:

            Walk[count] = i + 1

        elif N[i] > N[i-1]:

            N_ratio = N[i-1]/N[i]

            if pspin < N_ratio:

                i -= 1

                Walk[count] = i + 1

            elif pspin > N_ratio:

                Walk[count] = i + 1

        else:
            Walk[count] = i + 1

    elif Coin == 1:
        # head east

        if i == 6:

            Walk[count] = i + 1

        elif N[i] < N[i+1]:

            i += 1

            Walk[count] = i + 1

        else:
            Walk[count] = i + 1

    print Walk[count]
    count += 1


Days = np.linspace(1, t*Year, t*Year)

V_f = np.zeros(len(island_no))

V_fa= np.zeros(len(island_no))

for z in island_no:

    a = np.int(z)

    V_f[a-1] = len(Walk[Walk == a])/len(Walk)

    V_fa[a-1] = N[a-1]/np.sum(N)

plt.figure()
plt.semilogy(Walk, Days)


plt.figure()
plt.plot(island_no, V_f)
plt.plot(island_no, V_fa)

plt.show()