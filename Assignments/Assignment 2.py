# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt
from scipy.integrate import quad

plt.close('all')
np.random.seed(50)

#########################################################


# Question 1


def p_binom(v, N, theta):

    p = (theta**v)*((1-theta)**(N-v))

    return p


def p_beta(a, b, theta):

    p = (theta**(a-1))*((1-theta)**(b-1))

    return p


def pposterior(v, N, a, b, theta):

    post = (theta**((v + a)-1))*((1-theta)**((N - v + b)-1))

    norm_post = post/np.sum(post)

    return norm_post


v = 0; N = 3; a = 2; b = 2

theta_1 = np.linspace(0, 1, 10000)

likelihood = p_binom(v, N, theta_1)

mean_likelihood = (v + 1)/(N + 2)
std_likelihood = np.sqrt((mean_likelihood*(1 - mean_likelihood))/(N + 3))

print "Mean of likelihood =", mean_likelihood
print "Standard dev of likelihood =", std_likelihood

prior = p_beta(a, b, theta_1)

mean_prior = a/(a + b)
std_prior = np.sqrt((mean_prior*(1 - mean_prior))/(a + b + 1))

print "Mean of prior =", mean_prior
print "Standard dev of prior =", std_prior

posterior = pposterior(v, N, a, b, theta_1)

Posterior_Beta = p_beta(a + v, N - v + b, theta_1)

mean = (v + a)/(N + b + a)
std = np.sqrt((mean*(1 - mean))/(N + a + b + 1))

print "Mean of posterior =", mean
print "Standard dev of posterior =", std

#############################

# Plotting graphs

plt.figure()

plt.plot(theta_1, likelihood)
plt.xlabel(r'$ m $')
plt.ylabel(r'$ P(v|N,m) $')
#plt.title("Likelihood")

plt.figure()

plt.plot(theta_1, prior)
plt.xlabel(r'$ m $')
plt.ylabel(r'$ P(m|a,b) $')
#plt.title("Prior")

plt.figure()

plt.plot(theta_1, posterior)
plt.xlabel(r'$ m $')
plt.ylabel(r'$ P(m|v,N) $')
#plt.title("Posterior")

plt.figure()

plt.plot(theta_1, Posterior_Beta, label='Beta_distribution')
plt.legend(loc='best')
plt.xlabel(r'$ m $')
plt.ylabel(r'$ P(m|v,N) $')
#plt.title("Posterior")


#########################################################


# Question 2


age_data = np.array([2141.22, 1781.15, 1523.37, 1816.90, 1932.29, 1541.21, 720.782, 1026.22, 1687.55, 2460.59])

n = len(age_data)
mean_age = np.sum(age_data)/n
std_age = (np.sqrt((1/(n-1))*np.sum((age_data - mean_age)**2)))/np.sqrt(n)

mean_prior = 1200
std_prior = 300


def norm(x, mu, sig):

    fact = (1./(sig*np.sqrt(2*np.pi)))

    p = fact*(np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))

    return p


def norm_2(x, mu1, sig1, mu2, sig2):

    fact = (1./(sig1*np.sqrt(2*np.pi)))*(1./(sig2*np.sqrt(2*np.pi)))

    p = fact*(np.exp(-np.power(x - mu1, 2.) / (2 * np.power(sig1, 2.))))*(np.exp(-np.power(x - mu2, 2.) / (2 * np.power(sig2, 2.))))

    return p


def post_dist(age_data, mean_age, std_age, mean_prior, std_prior, time):

    min_age = np.min(age_data)
    max_age = np.max(age_data)

    X = np.linspace(min_age, max_age, time)

    N_age = norm(X, mean_age, std_age)
    N_prior = norm(X, mean_prior, std_prior)

    Sum = (np.sum(N_age*N_prior))

    N_post = N_age*N_prior/Sum

    return X, N_age, N_prior, N_post, Sum


def MCMC(time, step_size, data, mean_data, std_data, prior, distribution):

    count = 1  #  Start the counter

    decision = (np.random.random_integers(72, 246)*10)   ##  Starting at a random point between 720 and 2460

    i = np.int(decision)

    Walk = np.zeros(time)
    Prob = np.zeros(time)

    mean_prior = prior[0]
    std_prior = prior[1]

    while count < time:

        Walk[0] = i

        urand = np.random.uniform(0, 1, 1)

        sign = np.random.random_integers(0, 1)

        WO = Walk[count - 1]  #  Original theta

        PO = distribution(WO, mean_data, std_data, mean_prior, std_prior)  #  Original P(theta)

        Prob[count - 1] = PO

        if sign == 0:

            """ Move - step_size """

            WN = WO - step_size  #  Proposed theta

            PN = distribution(WN, mean_data, std_data, mean_prior, std_prior)  #  Proposed P(theta)

            Pmove = np.min([PN / PO, 1])

            if WO <= np.min(data):

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

            WN = WO + step_size  #  Proposed theta

            PN = distribution(WN, mean_data, std_data, mean_prior, std_prior)  #  Proposed P(theta)

            Pmove = np.min([PN / PO, 1])

            if WO >= np.max(data):

                Walk[count] = WO

            elif Pmove < 1:

                if urand <= Pmove:

                    Walk[count] = WN

                elif urand > Pmove:

                    Walk[count] = WO

            elif Pmove >= 1:

                Walk[count] = WN

        count += 1

    return Walk, Prob


def MCMC_int(X, data, mean_data, std_data, prior, h, distribution1, distribution2):

    mean_prior = prior[0]
    std_prior = prior[1]

    mean_h = h[0]
    std_h = h[1]

    P = (distribution1(X, mean_h, std_h)/distribution2(X, mean_data, std_data, mean_prior, std_prior))

    P_D = ((1 / len(P)) * np.sum(P))**-1

    return P_D


time = 1000000
step_size = 5

X, N_age, N_prior, N_post, Post_Sum = post_dist(age_data, mean_age, std_age, mean_prior, std_prior, time)

plt.figure()

plt.plot(X, N_post)
plt.xlabel('theta')
plt.ylabel(r'$P(\theta|\hat X)$')

prior = np.array([mean_prior, std_prior])

h = np.array([1544, 150])

theta_2, Prob_post = MCMC(time, step_size, age_data, mean_age, std_age, prior, norm_2)

converge_time_l = np.int64(time - (3*time/4))

theta_last = theta_2[converge_time_l:]   ##  taking the last n values of our MCMC sample

P_D = MCMC_int(theta_last, age_data, mean_age, std_age, prior, h, norm, norm_2)

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.hist(theta_2, bins=30)
ax1.set_xlabel('theta')
ax2.hist(theta_last, bins=30)
ax2.set_xlabel('theta')

MCMC_mean = np.mean(theta_last)
MCMC_std = np.std(theta_last)

print "MCMC mean =", MCMC_mean
print "MCMC std =", MCMC_std
print "Evidence term = ", P_D


#########################################################


# Question 3


X1 = 1.2  #  Initial data

GM1 = np.array([3, 0.2])  #  GM data
GM2 = np.array([1, 0.2])  #  GM data
O = np.array([2.2, 0.75])  #  Organic data


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

print "Odds of tomato being organic", P_OX


plt.show()
