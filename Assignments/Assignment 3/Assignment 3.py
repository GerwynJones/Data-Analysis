# -*- coding: utf-8 -*-
"""
Created on Wed Nov 08 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt


######################################


##  Question 1

print "Question 1"

x1, y1, sigma_y1 = np.loadtxt("PX4128_CA3_Q1.dat", skiprows=1, unpack=True)


n1 = len(x1)


##  Analytical


def Analyticaleq_HomoA(x, y, n):

    ResultA = ((np.sum(x**2.0) * np.sum(y)) - (np.sum(x) * np.sum(x * y))) / ((n * np.sum(x**2.0)) - (np.sum(x))**2.0)

    return ResultA


def Analyticaleq_HomoB(x, y, n):

    ResultB = ((n*np.sum(y*x)) - (np.sum(x)*np.sum(y)))/((n*np.sum(x**2.0)) - (np.sum(x))**2.0)

    return ResultB


Analytical_A = Analyticaleq_HomoA(x1, y1, n1)

Analytical_B = Analyticaleq_HomoB(x1, y1, n1)

print "Analytical method for A : ", Analytical_A

print "Analytical method for B : ", Analytical_B


##  Brute-Force

N = 500

ymin = np.min(y1)

ymax = np.max(y1)

A = np.linspace(ymin, ymax, N)

B = np.linspace(0, 5, N)


X_2 = np.zeros((len(A), len(B)))


def X_sum(A, B):

    Num = (y1 - A - B*x1)**2.0

    Den = sigma_y1**2.0

    Sum = np.sum(Num/Den)

    return Sum


for i in xrange(len(A)):

    for j in xrange(len(B)):

        X_2[i, j] = X_sum(A[i], B[j])


for i in xrange(len(A)):

    found = False

    for j in xrange(len(B)):

        if X_2[i, j] == np.min(X_2):

            d = n1 - 2  ##  degrees fo freedom = 2

            print "Chi-Squared : ", (1/d)*np.min(X_2)

            print "Brute-Force method for A : ", A[i]

            print "Brute-Force method for B : ", B[j]

            found = True

            break

        else:

            continue

    if found:

        break


##  Uncertainties


def Uncertainty_HomoA(x, sigma_y, n):

    Un_A = sigma_y[0]*np.sqrt(np.sum(x**2.0) / ((n*np.sum(x**2.0)) - (np.sum(x))**2.0))

    return Un_A


def Uncertainty_HomoB(x, sigma_y, n):

    Un_B = sigma_y[0]*np.sqrt(n / ((n*np.sum(x**2.0)) - (np.sum(x))**2.0))

    return Un_B


sigma_A = Uncertainty_HomoA(x1, sigma_y1, n1)

sigma_B = Uncertainty_HomoB(x1, sigma_y1, n1)


print "Uncertainty in A : ", sigma_A

print "Uncertainty in B : ", sigma_B


##  Plotting graph with errorbars and best-fit line


x_line = np.linspace(np.min(x1), np.max(x1), 100)

y_line = Analytical_A + Analytical_B*x_line

ls = 'None'

fig, ax = plt.subplots()

# standard error bars
ax.errorbar(x1, y1, yerr=sigma_y1, marker=".", linestyle=ls, capsize=2)

ax.plot(x_line, y_line)

plt.xlabel("x-data")
plt.ylabel("y-data")


## Question 2

print "Question 2"

x2, y2, sigma_y2 = np.loadtxt("PX4128_CA3_Q2.dat", skiprows=1, unpack=True)


n2 = len(y2)

##  Analytical


def Analyticaleq_HeteroA(x, y, sigma):

    w = 1 / sigma**2.0

    ResultA = ((np.sum(w*x**2.0) * np.sum(w*y)) - (np.sum(w*x) * np.sum(w * x * y))) / ((np.sum(w) * np.sum(w*x**2.0)) - (np.sum(w*x))**2.0)

    return ResultA


def Analyticaleq_HeteroB(x, y, sigma):

    w = 1 / sigma**2.0

    ResultB = ((np.sum(w)*np.sum(w * y * x)) - (np.sum(w*x)*np.sum(w*y)))/((np.sum(w)*np.sum(w*x**2.0)) - (np.sum(w*x))**2.0)

    return ResultB


##  Actual best fit parameter

A_actual = Analyticaleq_HeteroA(x2, y2, sigma_y2)

B_actual = Analyticaleq_HeteroB(x2, y2, sigma_y2)

print "Analytical A : ", A_actual

print "Analytical B : ", B_actual

##  Uncertainty


def Uncertainty_HeteroA(x, sigma_y):

    w = 1 / sigma_y**2.0

    Un_A = np.sqrt(np.sum(w*x**2.0) / ((np.sum(w)*np.sum(w*x**2.0)) - (np.sum(w*x))**2.0))

    return Un_A


def Uncertainty_HeteroB(x, sigma_y):

    w = 1 / sigma_y ** 2.0

    Un_B = np.sqrt(np.sum(w) / ((np.sum(w)*np.sum(w*x**2.0)) - (np.sum(w*x))**2.0))

    return Un_B


sigma_A = Uncertainty_HeteroA(x2, sigma_y2)

sigma_B = Uncertainty_HeteroB(x2, sigma_y2)

print "Uncertainty in A : ", sigma_A

print "Uncertainty in B : ", sigma_B

##  Randomly sampling

N = 50000

A_sample = np.zeros(N)

B_sample = np.zeros(N)


for i in xrange(N):

    y_sample = np.random.normal(y2, sigma_y2)

    A_sample[i] = Analyticaleq_HeteroA(x2, y_sample, sigma_y2)

    B_sample[i] = Analyticaleq_HeteroB(x2, y_sample, sigma_y2)


print "Mean in A : ", np.mean(A_sample)

print "Std in A : ", np.std(A_sample)

print "Mean in B : ", np.mean(B_sample)

print "Std in B : ", np.std(B_sample)


##  Plotting graph with errorbars and best-fit line


x_line = np.linspace(np.min(x2), np.max(x2), 100)

y_line = A_actual + B_actual*x_line

y_line_sample = A_sample[-1] + B_sample[-1]*x_line

ls = 'None'

fig, ax = plt.subplots()

# standard error bars
ax.errorbar(x2, y2, yerr=sigma_y2, marker=".", linestyle=ls, capsize=2, label="True-Data")

ax.plot(x_line, y_line, label="Best-Fit")

# standard error bars
ax.plot(x2, y_sample, marker=".", linestyle=ls, label="Sample-Data")

ax.plot(x_line, y_line_sample, label="Sample-Fit")

plt.legend(loc="best")
plt.xlabel("x-data")
plt.ylabel("y-data")

plt.show()