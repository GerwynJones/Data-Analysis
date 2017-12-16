# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""
from __future__ import division

import numpy as np
import scipy as sp
import pandas as pd
import scipy.interpolate as spi
import matplotlib.pyplot as plt
import uncertainties.umath as umath
from uncertainties import ufloat
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


print K


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


def Mx_func(Mc, i, f_Mx):

    M_sun = 1.989e30

    Mx = np.linspace(0.5, 10, 1e6)*M_sun

    for j in xrange(len(Mx)):

        q = Mc/Mx[j]

        result = (Mx[j]*(np.sin(i))**3)/(1 + q)**2

        if result >= f_Mx:

            Mx_actual = Mx[j]

            break

        else:

            continue

    return Mx_actual


P = 0.3440915


M_sun = 1.989e30


f_Mx = f_Mx_func(P, K)


incline = np.deg2rad(90)

Mc = 0.7*M_sun


Mx = Mx_func(Mc, incline, f_Mx)

print f_Mx/M_sun

print Mx/M_sun






plt.show()