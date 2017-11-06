# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt

plt.close('all')
np.random.seed(50)

#########################################################

# Question 1


def plikelihood(v, N):

    theta = np.linspace(0, 1, 100)

    p = (theta**v)*((1-theta)**(N-v))

    return p/np.sum(p)


def pprior(a, b):

    theta = np.linspace(0, 1, 100)

    p = (theta**(a-1))*((1-theta)**(b-1))

    return p


def pposterior(v, N, a, b):

    theta = np.linspace(0, 1, 100)

    post = (theta**((v + a)-1))*((1-theta)**((N - v + b)-1))

    norm_post = post/np.sum(post)

    return norm_post


v = 0; N = 3; a = 5; b = 7

theta = np.linspace(0, 1, 100)
likelihood = plikelihood(v, N)
prior = pprior(a, b)
posterior = pposterior(v, N, a, b)

mean = (v+a)/(N+b+a)
std = np.sqrt((mean*(1-mean))/(N+a+b))

print (1/len(posterior))*np.sum(posterior)

print "Mean =", mean
print "Standard dev =", std

#############################

# Plotting graphs

plt.figure()

plt.plot(theta, likelihood)

plt.figure()

plt.plot(theta, prior)

plt.figure()

plt.plot(theta, posterior)


#########################################################


# Question 2

age_data = np.array([2141.22, 1781.15, 1523.37, 1816.90, 1932.29, 1541.21, 720.782, 1026.22, 1687.55, 2460.59])

mean_age = np.mean(age_data)
std_age = np.std(age_data)

mean_prior = 1200
std_prior = 300


def norm(x, mu, sig):

    fact = (1./(sig*np.sqrt(2*np.pi)))

    n = fact*(np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))

    return n


def norm2(x, mu1, sig1, mu2, sig2):

    fact = (1./(sig1*np.sqrt(2*np.pi)))*(1./(sig2*np.sqrt(2*np.pi)))

    n = fact*(np.exp(-np.power(x - mu1, 2.) / (2 * np.power(sig1, 2.))))*(np.exp(-np.power(x - mu2, 2.) / (2 * np.power(sig2, 2.))))

    return n


def post_dist(age_data, mean_prior, std_prior):

    mean_age = np.mean(age_data)
    std_age = np.std(age_data)

    min_age = np.min(age_data)
    max_age = np.max(age_data)

    X = np.linspace(min_age, max_age, 5000)


    N_age = norm(X, mean_age, std_age)
    N_prior = norm(X, mean_prior, std_prior)

    N_post = N_age*N_prior/(np.sum(N_age*N_prior))

    return X, N_age, N_prior, N_post


def MCMC(time, step_size, data, prior, distribution):

    count = 1

    decision = raw_input("Starting point (betweem 720 and 2460) ?: ")

    i = np.int(decision)

    Walk = np.zeros(time)

    mean_data = np.mean(data)
    std_data = np.std(data)

    mean_prior = prior[0]
    std_prior = prior[1]

    while count < time:

        Walk[0]  = i

        urand = np.random.uniform(0, 1, 1)

        sign = np.random.random_integers(0, 1)

        WO = Walk[count - 1]  # Original theta

        PO = distribution(WO, mean_data, std_data, mean_prior, std_prior)  #  Original P(theta)

        if sign == 0:

            """ Move - step_size """

            WN = WO - step_size  # Proposed theta

            PN = distribution(WN, mean_data, std_data, mean_prior, std_prior)  #  Proposed P(theta)

            Pmove = np.min([PN / PO, 1])

            if WO <= np.min(data) :

                Walk[count] = WO

            elif Pmove < 1:

                if urand <= Pmove:

                    Walk[count] = WN

                elif urand > Pmove:

                    Walk[count] = WO

            elif Pmove >= 1:

                Walk[count] = WN

        elif sign == 1:

            """ Move + step_size """

            WN = WO + step_size  # Proposed theta

            PN = distribution(WN, mean_data, std_data, mean_prior, std_prior)  #  Proposed P(theta)

            Pmove = np.min([PN / PO, 1])

            if WO >= np.max(data) :

                Walk[count] = WO

            elif Pmove < 1:

                if urand <= Pmove:

                    Walk[count] = WN

                elif urand > Pmove:

                    Walk[count] = WO

            elif Pmove >= 1:

                Walk[count] = WN

        count += 1

    return Walk


def MCMC_int(X, data, prior, h, distribution1, distribution2):

    mean_data = np.mean(data)
    std_data = np.std(data)

    mean_prior = prior[0]
    std_prior = prior[1]

    mean_h = h[0]
    std_h = h[1]

    P = ((distribution1(X, mean_h, std_h))/(distribution2(X, mean_data, std_data, mean_prior, std_prior)))

    Sum = ((1/len(P))*(np.sum(P)))**-1.

    return Sum


X, N_age, N_prior, N_post = post_dist(age_data, mean_prior, std_prior)

plt.figure()

plt.plot(X, N_post)


print X[N_post == np.max(N_post)]


time = 200000
step_size = 10

prior = np.array([mean_prior, std_prior])

h = np.array([1345, 235])

Walk = MCMC(time, step_size, age_data, prior, norm2)

Sum = MCMC_int(Walk, age_data, prior, h, norm, norm2)

plt.figure()

plt.hist(Walk, bins=30, normed=True)

MCMC_mean = np.mean(Walk)
MCMC_std = np.std(Walk)

print MCMC_mean
print MCMC_std
print Sum


#########################################################


# Question 3


X1 = 1.2

GM1 = np.array([3, 0.2])
GM2 = np.array([1, 0.2])
O = np.array([2.2, 0.75])


def likelihood_GM(X1, GM1, GM2):

    P_GM1 = norm(X1, GM1[0], GM1[1])
    P_GM2 = norm(X1, GM2[0], GM2[1])

    Tot = 0.5*P_GM1 + 0.5*P_GM2

    return Tot

def likelihood_O(X1, O):

    P_O = norm(X1, O[0], O[1])

    return P_O

prior_GM = 0.7

prior_O = 0.3

P_GM = likelihood_GM(X1, GM1, GM2)*prior_GM

P_O = likelihood_O(X1, O)*prior_O

P_OX = P_O/P_GM

print P_OX


plt.show()
