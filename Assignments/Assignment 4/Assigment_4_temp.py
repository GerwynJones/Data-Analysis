# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import numpy as np
import scipy as sp
import scipy.interpolate as spi
import scipy.stats as sps
import matplotlib.pyplot as plt


plt.close('all')


def f(x, x_data, y_data, s_data, t_data, per_data):

    """ Cubic Spline """

    tck = spi.splrep(x_data, y_data, k=3, s=s_data, t=t_data, per=per_data)

    return spi.splev(x, tck)


def smooth(y_data, err, window_width):

    """ Smoothing data """

    y_smooth = np.zeros(len(y_data))

    for i in range(len(y_data)):

        if i > len(y_data) - window_width / 2:

            S = y_data[i - int(window_width / 2):]

            y_smooth[i] = np.sum(S)/len(S)

        elif i < window_width / 2:

            S = y_data[: i + int(window_width / 2) + 1]

            y_smooth[i] = np.sum(S)/len(S)

        else:

            S = y_data[i - int(window_width / 2): i + int(window_width / 2) + 1]

            y_smooth[i] = np.sum(S)/len(S)

    return y_smooth


def doppler_shift(lambda_original, velocity):

    c = 3e8

    """ Velocity in terms of c """

    lambda_new = lambda_original*(1 + velocity/c)

    return lambda_new


def sinusoidal_curve(gamma, Kx, Ky, Vel, phase):

    return Vel - gamma - Kx*np.sin(2*np.pi*phase) - Ky*np.cos(2*np.pi*phase)



template_list = ['g5', 'g9', 'k0', 'k1', 'k2', 'k4', 'k5', 'k7', 'k8', 'm0']

gs200_id = np.arange(1, 14)

gs200_size = gs200_id.size



W = str(template_list[6])

print W

lambda_temp_orig, flux_temp_orig, err_temp_orig = np.loadtxt(
    '/home/gerwyn/Documents/Physics/Computing Year 4/Data-Analysis/Assignments/Assignment 4/MiniProjectAllData/KeckTemplates/keck_' + W + '.txt',
    unpack=True)

# Define h:
window_width_temp = 10

# Create an array of smoothed data:
flux_temp_smooth = smooth(flux_temp_orig, err_temp_orig, window_width_temp)

# Define the number of knots:
nknots = 4

# Create the array of knots:
knots = np.arange(lambda_temp_orig[1], lambda_temp_orig[-1],
                  (lambda_temp_orig[-1] - lambda_temp_orig[1]) / np.double(nknots))

# Use cubic spline:
spl_temp = f(lambda_temp_orig, lambda_temp_orig, flux_temp_smooth, s_data=None, t_data=knots, per_data=0)

# Continuum subtracted template:
continuum_subtracted_temp = spl_temp - flux_temp_smooth

plt.figure()

plt.plot(lambda_temp_orig, spl_temp)

plt.plot(lambda_temp_orig, flux_temp_orig, label=W)

plt.legend(loc='best')

plt.figure()

plt.plot(lambda_temp_orig, continuum_subtracted_temp)


Q = str(2)

lambda_keck_gs2000_orig, flux_keck_gs2000_orig, err_keck_gs2000_orig = np.loadtxt(
    '/home/gerwyn/Documents/Physics/Computing Year 4/Data-Analysis/Assignments/Assignment 4/MiniProjectAllData//GS2000/keck_gs2000_0' + Q + '.txt',
    unpack=True)

# Define h:
window_width_keck = 10

flux_keck_gs2000_smooth = smooth(flux_keck_gs2000_orig, err_keck_gs2000_orig, window_width_keck)

# Create the array of knots:
knots = np.arange(lambda_keck_gs2000_orig[1], lambda_keck_gs2000_orig[-1],
                  (lambda_keck_gs2000_orig[-1] - lambda_keck_gs2000_orig[1]) / np.double(nknots))

# Use cubic spline:
spl_keck = f(lambda_keck_gs2000_orig, lambda_keck_gs2000_orig, flux_keck_gs2000_smooth, s_data=None, t_data=knots,
             per_data=0)

# Continuum subtracted gs200:
continuum_subtracted_keck = spl_keck - flux_keck_gs2000_smooth

plt.figure()

plt.plot(lambda_keck_gs2000_orig, spl_keck)

plt.plot(lambda_keck_gs2000_orig, flux_keck_gs2000_smooth, label=Q)

plt.legend(loc='best')

plt.figure()

plt.plot(lambda_keck_gs2000_orig, continuum_subtracted_keck)



Vel = np.linspace(-1e6, -3e5, 1001)

Chi_squared = np.zeros(len(Vel))

for i in xrange(len(Vel)):

    # Doppler shift

    doppler_lambda_temp = doppler_shift(lambda_temp_orig, Vel[i])


    # Interpolating template using cubic spline

    continuum_temp_shifted = f(lambda_keck_gs2000_orig, doppler_lambda_temp, continuum_subtracted_temp, s_data=None, t_data=None, per_data=0)


    # Reducing the arrays so the may be interpolated correctly through cubic splines and exclude difference in patterns above 6400

    # Keck arrays

    Wavelength_min = 5650

    Wavelength_max = 6900

    reduced_lambda_keck = lambda_keck_gs2000_orig[(lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

    reduced_continuum_keck = continuum_subtracted_keck[(lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

    reduced_err_keck_gs2000 = err_keck_gs2000_orig[(lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

    # Template arrays

    reduced_continuum_temp_shifted = continuum_temp_shifted[(lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]


    # Scaling parameter

    A_parameter = np.sum((reduced_continuum_keck*reduced_continuum_temp_shifted)/(reduced_err_keck_gs2000**2)) / np.sum((reduced_continuum_temp_shifted**2)/(reduced_err_keck_gs2000**2))

    sigma_A = np.sum((reduced_continuum_temp_shifted**2)/(reduced_err_keck_gs2000**2))

    # Chi Squared

    Chi_squared[i] = np.sum((reduced_continuum_keck - reduced_continuum_temp_shifted*A_parameter)**2/reduced_err_keck_gs2000**2)



plt.figure()

plt.plot(Vel, Chi_squared)

plt.figure()

plt.plot(reduced_lambda_keck, reduced_continuum_keck, label='GS2000')


plt.plot(reduced_lambda_keck, reduced_continuum_temp_shifted, label='Template')

plt.legend(loc='best')


plt.show()

