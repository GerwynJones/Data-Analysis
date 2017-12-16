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


def cubic_spline(x, x_data, y_data, s_data, t_data, per_data):

    """ Cubic Spline """

    tck = spi.splrep(x_data, y_data, k=3, s=s_data, t=t_data, per=per_data)

    return spi.splev(x, tck)


def smooth(y_data, window_width):

    """ Smoothing data """

    # window = np.ones(window_width)/window_width
    # y_smooth = np.convolve(y_data, window, mode='same')

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


# template_list = ['g5', 'g9', 'k0', 'k1', 'k2', 'k4', 'k5', 'k7', 'k8', 'm0']

template_list = ['k5']

gs200_id = np.arange(1, 14)

gs200_size = gs200_id.size


v = 1e6

Vel = np.linspace(-v, v, 101)

radial_vel = np.zeros((len(template_list), gs200_size))

error_rad_vel_less = np.zeros((len(template_list), gs200_size))

error_rad_vel_more = np.zeros((len(template_list), gs200_size))


phase = np.array([-0.1405, -0.0583, 0.0325, 0.0998, 0.1740, 0.2310, 0.3079, 0.3699, 0.4388, 0.5008, 0.5698, 0.6371, 0.7276])


Wavelength_min = 5700

Wavelength_max = 6300


# fig, ax = plt.subplots()

for k in range(len(template_list)):

    W = str(template_list[k])

    print W

    lambda_temp_orig, flux_temp_orig, err_temp_orig = np.loadtxt(
        '/home/gerwyn/Documents/Physics/Computing Year 4/Data-Analysis/Assignments/Assignment 4/MiniProjectAllData/KeckTemplates/keck_' + W + '.txt',
        unpack=True)

    # Define h:
    window_width = 5

    # Create an array of smoothed data:
    flux_temp_smooth = smooth(flux_temp_orig, window_width)

    # Define the number of knots:
    nknots = 4

    # Create the array of knots:
    knots = np.arange(lambda_temp_orig[1], lambda_temp_orig[-1],
                      (lambda_temp_orig[-1] - lambda_temp_orig[1]) / np.double(nknots))

    # Use cubic spline:
    spl_temp = cubic_spline(lambda_temp_orig, lambda_temp_orig, flux_temp_smooth, s_data=None, t_data=knots, per_data=0)

    # Continuum subtracted template:
    continuum_subtracted_temp = spl_temp - flux_temp_smooth

    plt.figure()

    plt.title("Template = %s" % W)

    for j in range(gs200_size):

        id = gs200_id[j]

        if id < 10:

            Q = str(id)

            lambda_keck_gs2000_orig, flux_keck_gs2000_orig, err_keck_gs2000_orig = np.loadtxt(
                '/home/gerwyn/Documents/Physics/Computing Year 4/Data-Analysis/Assignments/Assignment 4/MiniProjectAllData//GS2000/keck_gs2000_0' + Q + '.txt',
                unpack=True)

            flux_keck_gs2000_smooth = smooth(flux_keck_gs2000_orig, window_width)

            # Create the array of knots:
            knots = np.arange(lambda_keck_gs2000_orig[1], lambda_keck_gs2000_orig[-1],
                              (lambda_keck_gs2000_orig[-1] - lambda_keck_gs2000_orig[1]) / np.double(nknots))

            # Use cubic spline:
            spl_keck = cubic_spline(lambda_keck_gs2000_orig, lambda_keck_gs2000_orig, flux_keck_gs2000_smooth, s_data=None,
                         t_data=knots,
                         per_data=0)

            # Continuum subtracted gs200:
            continuum_subtracted_keck = spl_keck - flux_keck_gs2000_smooth

            Chi_squared = np.zeros(len(Vel))

            for i in xrange(len(Vel)):

                # Doppler shift

                doppler_lambda_temp = doppler_shift(lambda_temp_orig, Vel[i])

                # Interpolating template using cubic spline

                continuum_temp_shifted = cubic_spline(lambda_keck_gs2000_orig, doppler_lambda_temp, continuum_subtracted_temp, s_data=None, t_data=None, per_data=0)

                # Reducing the arrays so the may be interpolated correctly through cubic splines and exclude difference in patterns above 6400

                # Keck arrays

                reduced_lambda_keck = lambda_keck_gs2000_orig[
                    (lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

                reduced_continuum_keck = continuum_subtracted_keck[
                    (lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

                reduced_err_keck_gs2000 = err_keck_gs2000_orig[
                    (lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

                # Template arrays

                reduced_continuum_temp_shifted = continuum_temp_shifted[
                    (lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

                # Scaling parameter

                A_parameter = np.sum((reduced_continuum_keck * reduced_continuum_temp_shifted) / (
                        reduced_err_keck_gs2000 ** 2)) / np.sum(
                        (reduced_continuum_temp_shifted ** 2) / (reduced_err_keck_gs2000 ** 2))

                # Chi Squared

                Chi_squared[i] = np.sum((reduced_continuum_keck - reduced_continuum_temp_shifted * A_parameter) ** 2 /
                                        reduced_err_keck_gs2000 ** 2)

            # if id == 1 :
            # else:
            #
            #     continue

            plt.plot(Vel/(1e3), Chi_squared, label=id)

            plt.legend(loc='best')

            plt.xlabel(r"$Velocity \/ shift \/ (Km/s)$")
            plt.ylabel(r"$\chi^2$")

            Chi_min = np.min(Chi_squared)

            radial_vel[k, j] = Vel[Chi_squared == Chi_min]

            Chi_plus_one = Chi_min + 1

            less_than_chi = Chi_squared[Vel < radial_vel[k, j]]

            more_than_chi = Chi_squared[Vel > radial_vel[k, j]]

            less_than_vel = Vel[Vel < radial_vel[k, j]]

            more_than_vel = Vel[Vel > radial_vel[k, j]]

            error_rad_vel_less[k, j] = np.max(less_than_vel[less_than_chi > Chi_plus_one]) - radial_vel[k, j]

            error_rad_vel_more[k, j] = radial_vel[k, j] - np.min(more_than_vel[more_than_chi > Chi_plus_one])


        elif id >= 10:

            Q = str(id)

            lambda_keck_gs2000_orig, flux_keck_gs2000_orig, err_keck_gs2000_orig = np.loadtxt(
                '/home/gerwyn/Documents/Physics/Computing Year 4/Data-Analysis/Assignments/Assignment 4/MiniProjectAllData//GS2000/keck_gs2000_' + Q + '.txt',
                unpack=True)

            flux_keck_gs2000_smooth = smooth(flux_keck_gs2000_orig, window_width)

            # Create the array of knots:
            knots = np.arange(lambda_keck_gs2000_orig[1], lambda_keck_gs2000_orig[-1],
                              (lambda_keck_gs2000_orig[-1] - lambda_keck_gs2000_orig[1]) / np.double(nknots))

            # Use cubic spline:
            spl_keck = cubic_spline(lambda_keck_gs2000_orig, lambda_keck_gs2000_orig, flux_keck_gs2000_smooth, s_data=None,
                         t_data=knots,
                         per_data=0)

            # Continuum subtracted gs200:
            continuum_subtracted_keck = spl_keck - flux_keck_gs2000_smooth

            Chi_squared = np.zeros(len(Vel))

            for i in xrange(len(Vel)):

                # Doppler shift

                doppler_lambda_temp = doppler_shift(lambda_temp_orig, Vel[i])

                # Interpolating template using cubic spline

                continuum_temp_shifted = cubic_spline(lambda_keck_gs2000_orig, doppler_lambda_temp, continuum_subtracted_temp, s_data=None, t_data=None, per_data=0)

                # Reducing the arrays so the may be interpolated correctly through cubic splines and exclude difference in patterns above 6400

                # Keck arrays

                reduced_lambda_keck = lambda_keck_gs2000_orig[
                    (lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

                reduced_continuum_keck = continuum_subtracted_keck[
                    (lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

                reduced_err_keck_gs2000 = err_keck_gs2000_orig[
                    (lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

                # Template arrays

                reduced_continuum_temp_shifted = continuum_temp_shifted[
                    (lambda_keck_gs2000_orig >= Wavelength_min) & (lambda_keck_gs2000_orig <= Wavelength_max)]

                # Scaling parameter

                A_parameter = np.sum((reduced_continuum_keck * reduced_continuum_temp_shifted) / (
                        reduced_err_keck_gs2000 ** 2)) / np.sum(
                        (reduced_continuum_temp_shifted ** 2) / (reduced_err_keck_gs2000 ** 2))

                # Chi Squared

                Chi_squared[i] = np.sum((reduced_continuum_keck - reduced_continuum_temp_shifted * A_parameter) ** 2 /
                                        reduced_err_keck_gs2000 ** 2)


            plt.plot(Vel/(1e3), Chi_squared, label=id)

            plt.legend(loc='best')

            plt.xlabel(r"$Velocity \/ shift \/ (Km/s)$")
            plt.ylabel(r"$\chi^2$")


            Chi_min = np.min(Chi_squared)

            radial_vel[k, j] = Vel[Chi_squared == Chi_min]

            Chi_plus_one = Chi_min + 1

            less_than_chi = Chi_squared[Vel < radial_vel[k, j]]

            more_than_chi = Chi_squared[Vel > radial_vel[k, j]]

            less_than_vel = Vel[Vel < radial_vel[k, j]]

            more_than_vel = Vel[Vel > radial_vel[k, j]]

            error_rad_vel_less[k, j] = np.max(less_than_vel[less_than_chi > Chi_plus_one]) - radial_vel[k, j]

            error_rad_vel_more[k, j] = radial_vel[k, j] - np.min(more_than_vel[more_than_chi > Chi_plus_one])

    ls = 'None'

    fig, ax = plt.subplots()

    # standard error bars
    ax.errorbar(phase, radial_vel[k, :]/(1e3), yerr=[error_rad_vel_less[k, :]/(1e3), error_rad_vel_more[k, :]/(1e3)], marker=".", linestyle=ls, capsize=2, label=W)

    ax.legend(loc='best')

    ax.set_ylabel(r"$Velocity \/ (Km/s)$")
    ax.set_xlabel(r"$binary \/ phase$")

