# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Admin
"""

from __future__ import division

""" 30% and 10% of Republican and Independent voters are, respectively,
behind the change in the law, while 80% of the Democrat voters are in favour. You are visiting the
state, and ask a Police Officer what she thinks of the idea. She says she’s against the change to
the law. What is the probability that she votes Democrat? """


def P(A, B, C):
    """ Probability of three union """

    Probability = A * B * C

    return Probability * 100


# Probability she agrees with the disagreeing democrats,
# Disagrees with the disagreeing republicans and disagrees with the disagreeing independents

# Probability of agreeing democrats
DA = 0.8

# Probability of disagreeing democrats
DAc = 1 - DA

# Probability of agreeing independents
IA = 0.1

# Probability of disagreeing independents
IAc = 1 - IA

# Probability of agreeing republicans
RA = 0.3

# Probability of disagreeing republicans
RAc = 1 - RA

PovD = P(DAc, IA, RA)

print(PovD)
