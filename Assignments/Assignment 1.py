# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sc
from scipy.special import factorial
from scipy.stats.stats import pearsonr

plt.close('all')

#########################################################

# Question 1

""" 30% and 10% of Republican and Independent voters are, respectively,
behind the change in the law, while 80% of the Democrat voters are in favour. You are visiting the
state, and ask a Police Officer what she thinks of the idea. She says she’s against the change to
the law. What is the probability that she votes Democrat? """

# Finding the number of Democrat, Republican and Independent Voters from the State of Florida

# Source : https://www.nytimes.com/elections/results/florida

print("Question 1 :- ")

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

print('P(D|Fc) =', Probability)

###################################################

# Question 2

""" Roughly half of its latest batch of CPUs contains a flaw. How many CPUs from the batch would
they need to examine to know the probability that any given CPU is faulty to better than 5%?  """

print("Question 2 :- ")

p = 1/2
a = 0.05   # error in our accuracy
t = 1.96   # t value (95% confidence)

n_cpu = (t/a)**2*(1-p)*(p)

print("Number of CPU for probability of at least one being faulty to better than 5% is", n_cpu)

##################################################

# Question 3

""" A group researching cancer have previously found that the genetic marker D3 is a useful
indication that a person will develop the more aggressive form of melanoma skin cancer, in that D3
is present in 65% of the aggressive cases. However the test is expensive. A rival group claim that
the marker M23 is more sensitive than D3, and works out considerably cheaper to test for. The
rival research team manage to get DNA samples from 7 patients with the aggressive form of the
disease, all of whom test positive for the genetic marker M23. Based on these results, is M23 a
better marker for the disease than D3?  """

print("Question 3 :- ")


def B(N, p, v):

    T = np.linspace(0, N, N+1)

    fact = factorial(N)/(factorial(T)*factorial(N-T))

    B = fact*(p**T)*((1-p)**(N-T))

    Bino = np.sum(B[T >= v])

    return Bino


print("P(7, 0.65, 7) =", B(7, 0.65, 7))

###################################################

# Question 4

""" Eight new recruits for a rugby team are timed in both the 100 meters and 1,500 to assess
their athletic abilities """

# Source : https://stackoverflow.com/questions/3949226/calculating-pearson-correlation-and-significance-in-python

# And : https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.corrcoef.html

# And : http://www.socscistatistics.com/pvalues/Default.aspx

print("Question 4 :- ")

x_Data = np.array([12, 11, 13, 14, 12, 15, 12, 16])
y_Data = np.array([280, 290, 220, 260, 270, 240, 250, 230])


def R(x_Data, y_Data):

    Mean_x_Data = np.mean(x_Data)
    Mean_y_Data = np.mean(y_Data)

    r = np.sum((x_Data - Mean_x_Data)*(y_Data - Mean_y_Data))/np.sqrt((np.sum((x_Data - Mean_x_Data)**2))*(np.sum((y_Data - Mean_y_Data)**2)))

    return r


r = R(x_Data, y_Data)

print("r = ", r)

print("P8(|r| >", np.abs(r), ") = 5.706 %, using analytical function for r and a p-value calculator")

# Or use an inbuilt function

pr = pearsonr(x_Data, y_Data)

print("P8(|r| >", np.abs(pr[0]), ") = ", pr[1]*100, "%, using an inbuilt function for both r and p-value")

# Plotting data points

p_coeff,residuals,_,_,_= np.polyfit(x_Data, y_Data, 1, full=True)

p = np.poly1d(p_coeff)

x_trial = np.linspace(np.min(x_Data), np.max(x_Data), 100)

plt.plot(x_trial, p(x_trial), color='blue')

plt.plot(x_Data, y_Data, color='black', linestyle='None', marker='.')
plt.xlabel("100m")
plt.ylabel("1500m")

#####################################################

# Question 5

""" Using only a uniform random number generator, compute your own table of significance
values for linear correlation coefficient r. Do not use the analytic expression for r """

print("Question 5 :- ")

N = np.linspace(3, 10, 8)

num = 10000


def P_r(N, num):

    r = np.zeros((len(N), num))

    for i in range(len(N)):

        for j in range(num):

            X = np.random.uniform(0, 1, np.int(N[i]))
            Y = np.random.uniform(0, 1, np.int(N[i]))

            r[i, j] = np.abs(R(X, Y))

    return r


r = P_r(N, num)


def table(N, r):
    print("Table of probability of correlation due to chance")

    print("           r     =     0,    0.1,    0.2,   0.3,   0.4,   0.5,   0.6,   0.7,   0.8,   0.9,   1")

    r0 = np.linspace(0, 1, 11)

    P = np.zeros((len(N), len(r0)))

    for i in range(len(N)):

        for j in range(len(r0)):

            P[i, j] = np.int64((len(r[i, :][r[i, :] > r0[j]])/num)*100)

        print("N = ", N[i], P[i, :])

    return P


P_Corr = table(N, r)



plt.show()