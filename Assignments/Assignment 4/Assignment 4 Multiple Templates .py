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
import uncertainties as u
from uncertainties import unumpy


plt.close('all')


@u.wrap
def roots(n=0, *P):

    """ Adapted the numpy.roots function to take into account the uncertainties in the Coefficients """

    p = np.array(P)
    N = len(p)

    M = np.diag(np.ones((N - 2,), p.dtype), -1)
    M[0, :] = -p[1:] / p[0]

    r = np.linalg.eigvals(M)

    return r[n]


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

        C1 = u.ufloat(c1, sc1)
        C2 = u.ufloat(c2, sc2)
        C3 = u.ufloat(c3, sc3)
        C4 = u.ufloat(c4, sc4)

        Coeff = [C1, C2, C3, C4]


        for m in xrange(len(Coeff) - 1):

            Root = roots(m, *Coeff)

            if Root.nominal_value == Root_very_real:

                Result = Root

            else:

                continue

    return Result


template_list = ['g5', 'g9', 'k0', 'k1', 'k2', 'k4', 'k5', 'k7', 'k8', 'm0']

gs200_id = np.arange(1, 14)

gs200_size = gs200_id.size


v = 1e6

Vel = np.linspace(-v, v, 201)

radial_vel = np.zeros((len(template_list), gs200_size))

error_rad_vel_less = np.zeros((len(template_list), gs200_size))

error_rad_vel_more = np.zeros((len(template_list), gs200_size))


Wavelength_min = 5700

Wavelength_max = 6200


fig, ax = plt.subplots()

for k in range(len(template_list)):

    W = str(template_list[k])

    lambda_temp_orig, flux_temp_orig, err_temp_orig = np.loadtxt(
        '/home/gerwyn/Documents/Physics/Computing Year 4/Data-Analysis/Assignments/Assignment 4/MiniProjectAllData/KeckTemplates/keck_' + W + '.txt',
        unpack=True)

    # Define h:
    window_width = 2

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

                ax.plot(Vel / (1e3), Chi_squared, label=W)

                ax.set_xlabel(r"$Velocity \/ shift \/ (Km/s)$")
                ax.set_ylabel(r"$\chi^2$")

            else:

                Carry_on = 0


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

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    Average_sig = (error_rad_vel_more + error_rad_vel_less) / 2

plt.show()