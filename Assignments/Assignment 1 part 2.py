# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc

# Question 1

""" 30% and 10% of Republican and Independent voters are, respectively,
behind the change in the law, while 80% of the Democrat voters are in favour. You are visiting the
state, and ask a Police Officer what she thinks of the idea. She says she’s against the change to
the law. What is the probability that she votes Democrat? """

##########################################################


# Finding the number of Democrat, Republican and Independent Voters from the State of Florida

# Source : https://www.nytimes.com/elections/results/florida


# Democrat

VoD = 4504975

# Republican

VoR = 4617886

# Independent

VoI = 81731

# 3 Party Total

To3 = VoD + VoR + VoI

# Others

VoO = 297025

# Total Overall

ToV = VoD + VoR + VoI + VoO


# From the Question we can find the number of voters who would be for and against the change in the law

# Voters Who are For and a Democrat

DVF = 0.8

DF = np.int(DVF*VoD)

# Voters Against and a Democrat

DVFc = (1-DVF)

DFc = np.int(DVFc*VoD)

# Voters Who are For and a Republican

RVF = 0.3

RF = np.int(RVF*VoR)

# Voters Against and a Republican

RVFc = (1-RVF)

RFc = np.int(RVFc*VoR)

# Voters Who are For and a Independent

IVF = 0.1

IF = np.int(IVF*VoI)

# Voters Against and a Independent

IVFc = (1-IVF)

IFc = np.int(IVFc*VoI)


# Average For

AF = (DF + RF + IF)/To3

# Average Against

AA = (DFc + RFc + IFc)/To3


# Others For

OF = np.int(AF*VoO)

# Others Against

OFc = np.int(AA*VoO)


# Total For

TF = DF + RF + IF + OF

# Total Against

TFc = DFc + RFc + IFc + OFc


# Defining our functions

def P(A, B):

    """ Probabibility of geting A out of B """

    return A/B


# Probability of voting Democrat

PD = P(VoD, ToV)

# Probability of voting Against the change

PFc = P(TFc, ToV)


# Probability of Voting against given you vote democrat is just the probability DVFc

PFcgD = DVFc


# Therefore using Bayes Theorem, the probability of voting democrat given your against the change is :-

Probability = (PFcgD*PD)/PFc

print Probability

###############################################

# Question 2

""" Roughly half of its latest batch of CPUs contains a flaw. How many CPUs from the batch would
they need to examine to know the probability that any given CPU is faulty to better than 5%?  """

p = 1/2
a = 0.05

N = 100000

x = np.linspace(0, N, N + 1)

x1 = sc.binom.pmf(x, N, p)

#b = np.min(x[x1 > a])/N

print TF/ToV
print TFc/ToV

