# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc

x = np.linspace(0,10,11)

q = sc.binom.pmf(x, 10, .5)

p = 1 - q[0] - q[1] - q[2] - q[3] - q[4] - q[5]

print p


print np.min([101, 2])
