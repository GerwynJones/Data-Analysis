# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Admin
"""

from __future__ import division


# Question 1

""" 30% and 10% of Republican and Independent voters are, respectively,
behind the change in the law, while 80% of the Democrat voters are in favour. You are visiting the
state, and ask a Police Officer what she thinks of the idea. She says she’s against the change to
the law. What is the probability that she votes Democrat? """

##########################################################

# Defining our functions

def P(A, B, C):
    """ Probability of disagreeing and voting democrat"""

    PB = (A + B + C)

    Probability = A/PB

    return Probability * 100


def Normalize(A, B, C, D, E, F):

    sum = (A + B + C + D + E + F)

    return A/sum, B/sum, C/sum, D/sum, E/sum, F/sum

########################################################

# Probability she votes democrat given that she disagrees

# Probability of agreeing democrats
dema = 0.8

# Probability of disagreeing democrats
demac = 1 - dema

# Probability of agreeing independents
inda = 0.1

# Probability of disagreeing independents
indac = 1 - inda

# Probability of agreeing republicans
repa = 0.3

# Probability of disagreeing republicans
repac = 1 - repa

################################################

DA, DAc, IA, IAc, RA, RAc = Normalize(dema, demac, inda, indac, repa, repac)

PovD = P(DAc, IAc, RAc)

print(PovD)

###############################################

# Question 2






