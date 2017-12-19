# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""
from __future__ import division

import numpy as np
import scipy as sp
import scipy.stats as sps
import matplotlib.pyplot as plt
import uncertainties.umath as umath
import uncertainties as u
from uncertainties import unumpy


from Assignment_4 import radial_vel, error_rad_vel_less, error_rad_vel_more, phase


plt.close('all')


Average_sig = (error_rad_vel_more + error_rad_vel_less)/2


w = 1/Average_sig**2

W = np.sum(w)

Wsin = np.sum(w*np.sin(2*np.pi*phase))

Wcos = np.sum(w*np.cos(2*np.pi*phase))

Wsin2 = np.sum(w*(np.sin(2*np.pi*phase))**2)

Wcos2 = np.sum(w*(np.cos(2*np.pi*phase))**2)

Wcossin = np.sum(w*np.cos(2*np.pi*phase)*np.sin(2*np.pi*phase))


Hessian_Matrix = np.array([W, Wsin, Wcos, Wsin, Wsin2, Wcossin, Wcos, Wcossin, Wcos2]).reshape(3, 3)

vel = unumpy.uarray(radial_vel[0, :], np.abs(Average_sig))


Wv = np.sum(w*vel)

Wvsin = np.sum(w*np.sin(2*np.pi*phase)*vel)

Wvcos = np.sum(w*np.cos(2*np.pi*phase)*vel)


Vel_Matrix = np.array([Wv, Wvsin, Wvcos]).reshape(3, 1)


Inverse_Hessian = np.linalg.inv(Hessian_Matrix)


Dot_Product = np.dot(Inverse_Hessian, Vel_Matrix)


def sinusoidal_curve(gamma, Kx, Ky, phase):

    return gamma + Kx*np.sin(2*np.pi*phase) + Ky*np.cos(2*np.pi*phase)


gamma = Dot_Product[0]/(1e3)

Kx = Dot_Product[1]/(1e3)

Ky = Dot_Product[2]/(1e3)


phase_trial = np.linspace(np.min(phase), np.max(phase), 100)

sinusoidal_vel = sinusoidal_curve(gamma, Kx, Ky, phase_trial)


K = (Kx**2 + Ky**2)**(1/2)



ls = 'None'

fig, ax = plt.subplots()

# standard error bars
ax.errorbar(phase, radial_vel[0, :]/(1e3), yerr=[error_rad_vel_less[0, :]/(1e3), error_rad_vel_more[0, :]/(1e3)],
            marker=".", linestyle=ls, capsize=2)

ax.legend(loc='best')

ax.set_ylabel(r"$Velocity \/ (Km/s)$")
ax.set_xlabel(r"$binary \/ phase$")

Curve_val = unumpy.nominal_values(sinusoidal_vel)

ax.plot(phase_trial, Curve_val)


def f_Mx_func(P, K):

    P_d = P*(24*60*60)

    G = 6.67e-11

    return (P_d*(K*(1e3))**3)/(2*np.pi*G)




def Mx_func(Mc, i, f_Mx, Decision=0):

    f_Mx_val = unumpy.nominal_values(f_Mx)

    c1 = (np.sin(i))**3

    c2 = -f_Mx_val

    c3 = -2*f_Mx_val

    c4 = -Mc**2*f_Mx_val

    Coeff_real = [c1, c2, c3, c4]

    Root = np.roots(Coeff_real)

    Root_very_real = Root[Root.imag == 0].real


    if Decision == 0:

        Result = Root_very_real

    else:

        f_Mx_std = unumpy.std_devs(f_Mx)

        sc1 = 0
        sc2 = f_Mx_std
        sc3 = 2*f_Mx_std
        sc4 = f_Mx_std

        C1 = u.ufloat((c1, sc1))
        C2 = u.ufloat((c2, sc2))
        C3 = u.ufloat((c3, sc3))
        C4 = u.ufloat((c4, sc4))

        Coeff = [C1, C2, C3, C4]

        @u.wrap
        def f(n=0, *P):
            ''' compute the nth root of the polynomial P and the uncertainty of the root'''
            p = np.array(P)
            N = len(p)

            M = np.diag(np.ones((N - 2,), p.dtype), -1)
            M[0, :] = -p[1:] / p[0]

            r = np.linalg.eigvals(M)
            r.sort()  # there is no telling what order the values come out in
            return r[n]


        for m in xrange(len(Coeff) - 1):

            Root = f(m, *Coeff)

            if Root.nominal_value == Root_very_real:

                print Root

                Result = Root

            else:

                continue

    return Result


P = 0.3440915


M_sun = 1.989e30


f_Mx = f_Mx_func(P, K)


incline = np.deg2rad(90)

Mc = 0.7*M_sun


Mx = Mx_func(Mc, incline, f_Mx, 1)

print f_Mx/M_sun

print Mx/M_sun

Ntrial = 1e4

# Creating standard deviations for each parameter


standard_Mc = 0.3 * M_sun
standard_K = unumpy.std_devs(K)
standard_P = 0.1

angle = np.linspace(np.deg2rad(50), np.deg2rad(90), 1e6)

# Creating arrays for paramaters and Mx

incline_trial = np.zeros(np.int(Ntrial))
Mc_trial = np.zeros(np.int(Ntrial))
K_trial = np.zeros(np.int(Ntrial))
P_trial = np.zeros(np.int(Ntrial))

f_Mx_trial = np.zeros(np.int(Ntrial))
Mx_trial = np.zeros(np.int(Ntrial))

inc = np.deg2rad(80)

for l in xrange(np.int(Ntrial)):
    incline_trial[l] = np.random.choice(angle)

    Mc_trial[l] = np.random.normal(Mc, standard_Mc)

    K_trial[l] = np.random.normal(unumpy.nominal_values(K), standard_K)

    P_trial[l] = np.random.normal(P, standard_P)

    f_Mx_trial[l] = f_Mx_func(P_trial[l], K_trial[l])

    Mx_trial[l] = Mx_func(Mc_trial[l], incline_trial[l], f_Mx_trial[l], 0)

print np.mean(Mx_trial) / M_sun

print np.std(Mx_trial) / M_sun

plt.figure()

plt.hist(Mx_trial / M_sun, bins=25)




hist = np.histogram(Mx_trial / M_sun, bins=50)
hist_dist = sps.rv_histogram(hist)


X = np.linspace(-1.0, 20.0, 50)



plt.figure()

plt.title("PDF from Template")
plt.hist(Mx_trial / M_sun, normed=True, bins=100)
plt.plot(X, hist_dist.pdf(X), label='PDF')
plt.plot(X, hist_dist.cdf(X), label='CDF')




plt.show()