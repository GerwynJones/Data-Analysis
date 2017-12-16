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
from uncertainties import ufloat
from uncertainties import unumpy


plt.close('all')


def cubic_spline(x, x_data, y_data, s_data, t_data, per_data):

    """ Cubic Spline """

    tck = spi.splrep(x_data, y_data, k=3, s=s_data, t=t_data, per=per_data)

    return spi.splev(x, tck)


def smooth(y_data, window_width):

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


def sinusoidal_curve(gamma, Kx, Ky, phase):

    return gamma + Kx * np.sin(2 * np.pi * phase) + Ky * np.cos(2 * np.pi * phase)


def f_Mx_func(P, K):

    P_d = P*(24*60*60)

    G = 6.67e-11

    return (P_d*(K*1e3)**3)/(2*np.pi*G)


def Mx_func(Mc, i, f_Mx, N):

    M_sun = 1.989e30

    Mx = np.linspace(0.5, 10, N)*M_sun

    Mx_actual = 0

    for j in xrange(len(Mx)):

        q = Mc/Mx[j]

        result = (Mx[j]*(np.sin(i))**3)/(1 + q)**2

        if result >= f_Mx:

            Mx_actual = Mx[j]

            break

        else:

            continue

    return Mx_actual


# template_list = ['g5', 'g9', 'k0', 'k1', 'k2', 'k4', 'k5', 'k7', 'k8', 'm0']

template_list = ['k5']

gs200_id = np.arange(1, 14)

gs200_size = gs200_id.size


v = 1e6

Vel = np.linspace(-v, v, 201)

radial_vel = np.zeros((len(template_list), gs200_size))

error_rad_vel_less = np.zeros((len(template_list), gs200_size))

error_rad_vel_more = np.zeros((len(template_list), gs200_size))


phase = np.array([-0.1405, -0.0583, 0.0325, 0.0998, 0.1740, 0.2310, 0.3079, 0.3699, 0.4388, 0.5008, 0.5698, 0.6371, 0.7276])


Wavelength_min = 5700

Wavelength_max = 6300


# fig, ax = plt.subplots()

for k in range(len(template_list)):

    W = str(template_list[k])

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

    plt.plot(lambda_temp_orig, spl_temp, label='Overall continuum shape')

    plt.plot(lambda_temp_orig, flux_temp_orig, label='keck_k5')

    plt.legend(loc='best')

    plt.figure()

    plt.plot(lambda_temp_orig, continuum_subtracted_temp, label='Template continuum subtracted')

    plt.legend(loc='best')


    fig1, ax1 = plt.subplots()

    ax1.set_title("Template = %s" % W)

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

            A_parameter = np.zeros(len(Vel))
            sigma_A = np.zeros(len(Vel))
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

                A_parameter[i] = np.sum((reduced_continuum_keck * reduced_continuum_temp_shifted) / (
                        reduced_err_keck_gs2000 ** 2)) / np.sum(
                        (reduced_continuum_temp_shifted ** 2) / (reduced_err_keck_gs2000 ** 2))

                sigma_A[i] = 1 / np.sum((reduced_continuum_temp_shifted ** 2) / (reduced_err_keck_gs2000 ** 2))

                # Chi Squared

                Chi_squared[i] = np.sum((reduced_continuum_keck - reduced_continuum_temp_shifted * A_parameter[i]) ** 2 /
                                        reduced_err_keck_gs2000 ** 2)

            if id == 1:

                index = np.argmin(Chi_squared)

                A = unumpy.uarray([A_parameter[index]], [np.abs(sigma_A[index])])

                print A


                fig2, ax2 = plt.subplots()

                ax2.plot(lambda_keck_gs2000_orig, spl_keck, label='Overall continuum shape')

                ax2.plot(lambda_keck_gs2000_orig, flux_keck_gs2000_smooth, label='keck_gs2000_01')

                ax2.legend(loc='best')

                fig4, ax4 = plt.subplots()

                ax4.plot(lambda_keck_gs2000_orig, continuum_subtracted_keck, label='GS2000 continuum subtracted')

                ax4.legend(loc='best')

            else:

                Carry_on = 0


            ax1.plot(Vel/(1e3), Chi_squared, label=id)

            ax1.legend(loc='best')

            ax1.set_xlabel(r"$Velocity \/ shift \/ (Km/s)$")
            ax1.set_ylabel(r"$\chi^2$")

            index = np.argmin(Chi_squared)

            Chi_min = Chi_squared[index]

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

            A_parameter = np.zeros(len(Vel))
            sigma_A = np.zeros(len(Vel))
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

                A_parameter[i] = np.sum((reduced_continuum_keck * reduced_continuum_temp_shifted) / (
                        reduced_err_keck_gs2000 ** 2)) / np.sum(
                        (reduced_continuum_temp_shifted ** 2) / (reduced_err_keck_gs2000 ** 2))

                sigma_A[i] = 1 / np.sum((reduced_continuum_temp_shifted ** 2) / (reduced_err_keck_gs2000 ** 2))

                # Chi Squared

                Chi_squared[i] = np.sum((reduced_continuum_keck - reduced_continuum_temp_shifted * A_parameter[i]) ** 2 /
                                        reduced_err_keck_gs2000 ** 2)


            ax1.plot(Vel/(1e3), Chi_squared, label=id)

            ax1.legend(loc='best')

            ax1.set_xlabel(r"$Velocity \/ shift \/ (Km/s)$")
            ax1.set_ylabel(r"$\chi^2$")


            index = np.argmin(Chi_squared)

            Chi_min = Chi_squared[index]

            radial_vel[k, j] = Vel[Chi_squared == Chi_min]

            Chi_plus_one = Chi_min + 1

            less_than_chi = Chi_squared[Vel < radial_vel[k, j]]

            more_than_chi = Chi_squared[Vel > radial_vel[k, j]]

            less_than_vel = Vel[Vel < radial_vel[k, j]]

            more_than_vel = Vel[Vel > radial_vel[k, j]]

            error_rad_vel_less[k, j] = np.max(less_than_vel[less_than_chi > Chi_plus_one]) - radial_vel[k, j]

            error_rad_vel_more[k, j] = radial_vel[k, j] - np.min(more_than_vel[more_than_chi > Chi_plus_one])


    Average_sig = (error_rad_vel_more + error_rad_vel_less) / 2


    w = 1 / Average_sig ** 2

    W = np.sum(w)

    Wsin = np.sum(w * np.sin(2 * np.pi * phase))

    Wcos = np.sum(w * np.cos(2 * np.pi * phase))

    Wsin2 = np.sum(w * (np.sin(2 * np.pi * phase)) ** 2)

    Wcos2 = np.sum(w * (np.cos(2 * np.pi * phase)) ** 2)

    Wcossin = np.sum(w * np.cos(2 * np.pi * phase) * np.sin(2 * np.pi * phase))


    Hessian_Matrix = np.array([W, Wsin, Wcos, Wsin, Wsin2, Wcossin, Wcos, Wcossin, Wcos2]).reshape(3, 3)

    vel = unumpy.uarray(radial_vel[0, :], np.abs(Average_sig))

    Wv = np.sum(w * vel)

    Wvsin = np.sum(w * np.sin(2 * np.pi * phase) * vel)

    Wvcos = np.sum(w * np.cos(2 * np.pi * phase) * vel)

    Vel_Matrix = np.array([Wv, Wvsin, Wvcos]).reshape(3, 1)

    Inverse_Hessian = np.linalg.inv(Hessian_Matrix)

    Dot_Product = np.dot(Inverse_Hessian, Vel_Matrix)


    gamma = Dot_Product[0] / (1e3)

    Kx = Dot_Product[1] / (1e3)

    Ky = Dot_Product[2] / (1e3)

    phase_trial = np.linspace(np.min(phase), np.max(phase), 100)

    sinusoidal_vel = sinusoidal_curve(gamma, Kx, Ky, phase_trial)

    K = (Kx ** 2 + Ky ** 2) ** (1 / 2)

    ls = 'None'

    fig, ax = plt.subplots()

    # standard error bars
    ax.errorbar(phase, radial_vel[0, :] / (1e3),
                yerr=[error_rad_vel_less[0, :] / (1e3), error_rad_vel_more[0, :] / (1e3)],
                marker=".", linestyle=ls, capsize=2)

    ax.set_ylabel(r"$Velocity \/ (Km/s)$")
    ax.set_xlabel(r"$binary \/ phase$")

    Curve_val = unumpy.nominal_values(sinusoidal_vel)

    ax.plot(phase_trial, Curve_val)


    K = (Kx ** 2 + Ky ** 2) ** (1 / 2)

    P = 0.3440915  # days

    M_sun = 1.989e30


    f_Mx = f_Mx_func(P, K)


    incline = np.deg2rad(90)

    Mc = 0.7*M_sun

    N = 1e6

    Mx = Mx_func(Mc, incline, f_Mx, N)


    print K

    print f_Mx/M_sun

    print Mx/M_sun



plt.show()
