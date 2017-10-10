# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

"""
1) A new gun is expected to be able to fire its bullet to around 200m. However the manufacturers have identified 14 independent physical phenomena, each of which could affect the distance of the bullet by as much as +/- 8 meters. 
Using a uniform random number to mimic each of the physical phenomenon, show how one might expect the distribution of bullet distances to look after 1000 fires of the gun. (You need to make a plot!)
Calculate the mean and standard deviation of the distribution.
What happens if the number of random phenomena affecting the gun is actually only 4.
Discuss the shapes of the distribution with your classmates, and discuss why they look they way they do
Report what fraction of the distances are above 190m in each case.
"""

################################################

plt.close('all')

dist = 200
ierr = 8
ipp = 14
num = 1000

distribution = np.zeros(num)

for j in xrange(num):
    
    X = np.random.uniform (-ierr , ierr, ipp)
        
    L = np.sum(X)
    
    distribution[j] = dist + L

plt.figure()
plt.hist(distribution, bins = 20)


################################################

m = np.mean(distribution)
s = np.std(distribution)

################################################

new_distribution = np.zeros(num)

new_ipp = 4

for j in xrange(num):
        
    X = np.random.uniform (-ierr , ierr, new_ipp)
        
    L = np.sum(X)
    
    new_distribution[j] = dist + L

plt.figure()
plt.hist(new_distribution, bins = 20)

above = distribution[distribution > 190]

new_above = new_distribution[new_distribution > 190]

print len(above)
print len(new_above)

#################################################