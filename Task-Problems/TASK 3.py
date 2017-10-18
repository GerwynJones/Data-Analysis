# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


""" Astronomers like to cling onto archaic measurements, and one classic example is the use of magnitudes to determine 
the brightness of celestial objects. Magnitudes arose from the eye’s logarithmic response to light, but now that 
computers are here to measure flux (and because we want to do physics!), these magnitudes need to be converted to 
actual fluxes. The magnitude scale is connected to the actual flux f from the object via, m = m0 - 2.5log10f """



# Question 1

""" 1) Assuming the flux is normally distributed, N(f0, σf), derive the expression for the mean, 
and standard deviation of the magnitude distribution that would arise from the flux distribution. """

# Mean of M

# M(f) = f0 + (5sigf^2)/(4f^2ln10)

# STD of M

# sigM = sqrt((2.5/fln10)^2 sigf^2)

# Question 2

""" Confirm your expression using a Gaussian random number generator (units are not  important). """

M0 = 10
f0 = 2
sigf = 0.2



f = np.random.uniform(f0, sigf, 1000000)

M = M0 - 2.5*np.log10(f)


plt.hist(M)

plt.show()