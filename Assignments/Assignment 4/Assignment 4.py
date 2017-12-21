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

    """ Velocity in terms of ms-1 """

    c = 3e8

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

# setting up the template star

template_list = ['k5']

gs200_id = np.arange(1, 14)

gs200_size = gs200_id.size


# Creating possible trial velocities

v = 1e6

Vel = np.linspace(-v, v, 201)

# Creating arrays for the radial velocities and their errors

radial_vel = np.zeros(gs200_size)

error_rad_vel_less = np.zeros(gs200_size)

error_rad_vel_more = np.zeros(gs200_size)

# Creating an array for the companion star phase

phase = np.array([-0.1405, -0.0583, 0.0325, 0.0998, 0.1740, 0.2310, 0.3079, 0.3699, 0.4388, 0.5008, 0.5698, 0.6371, 0.7276])

# Choosing cut-off wavelengths

Wavelength_min = 5700

Wavelength_max = 6200


linestyle = ['solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'dashdot', 'dashdot', 'dashdot']

# Reading template file

W = str(template_list[0])

lambda_temp_orig, flux_temp_orig, err_temp_orig = np.loadtxt(
    '/home/gerwyn/Documents/Physics/Computing Year 4/Data-Analysis/Assignments/Assignment 4/MiniProjectAllData/KeckTemplates/keck_' + W + '.txt',
    unpack=True)

# Define h:
window_width = 0    # 0 = Not Smoothed

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


# Plotting spectra

plt.figure()

plt.plot(lambda_temp_orig, spl_temp)

plt.plot(lambda_temp_orig, flux_temp_orig)


plt.xlabel(r"$Wavelength \/ (\AA)$")

plt.figure()

plt.plot(lambda_temp_orig, continuum_subtracted_temp)

plt.xlabel(r"$Wavelength \/ (\AA)$")


fig1, ax1 = plt.subplots()

for j in range(gs200_size):

    id = gs200_id[j]

    if id < 10:

        # Reading GS2000 files

        Q = str(id)

        lambda_keck_gs2000_orig, flux_keck_gs2000_orig, err_keck_gs2000_orig = np.loadtxt(
            '/home/gerwyn/Documents/Physics/Computing Year 4/Data-Analysis/Assignments/Assignment 4/MiniProjectAllData//GS2000/keck_gs2000_0' + Q + '.txt',
            unpack=True)

        # Create an array of smoothed data:
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

        # Creating arrays for the parameters and their errors
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

            A = u.ufloat(A_parameter[index], np.abs(sigma_A[index]))

            print 'A parameter : {:.2u}'.format(A)

            # Plotting spectra

            fig2, ax2 = plt.subplots()

            ax2.plot(lambda_keck_gs2000_orig, spl_keck)

            ax2.plot(lambda_keck_gs2000_orig, flux_keck_gs2000_smooth)

            ax2.set_xlabel(r"$Wavelength \/ (\AA)$")

            fig4, ax4 = plt.subplots()

            ax4.plot(lambda_keck_gs2000_orig, continuum_subtracted_keck)

            ax4.set_xlabel(r"$Wavelength \/ (\AA)$")

        else:

            Carry_on = 0

        # Plotting chi-squared

        ax1.plot(Vel/(1e3), Chi_squared, label=phase[j], linestyle=linestyle[j])

        ax1.set_xlabel(r"$Velocity \/ shift \/ (Km/s)$")
        ax1.set_ylabel(r"$\chi^2$")

        index = np.argmin(Chi_squared)

        Chi_min = Chi_squared[index]

        radial_vel[j] = Vel[Chi_squared == Chi_min]

        Chi_plus_one = Chi_min + 1

        less_than_chi = Chi_squared[Vel < radial_vel[j]]

        more_than_chi = Chi_squared[Vel > radial_vel[j]]

        less_than_vel = Vel[Vel < radial_vel[j]]

        more_than_vel = Vel[Vel > radial_vel[j]]

        error_rad_vel_less[j] = np.max(less_than_vel[less_than_chi > Chi_plus_one]) - radial_vel[j]

        error_rad_vel_more[j] = radial_vel[j] - np.min(more_than_vel[more_than_chi > Chi_plus_one])


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

        # Plotting chi-squared

        ax1.plot(Vel/(1e3), Chi_squared, label=phase[j], linestyle=linestyle[j])

        ax1.set_xlabel(r"$Velocity \/ shift \/ (Km/s)$")
        ax1.set_ylabel(r"$\chi^2$")


        index = np.argmin(Chi_squared)

        Chi_min = Chi_squared[index]

        radial_vel[j] = Vel[Chi_squared == Chi_min]

        Chi_plus_one = Chi_min + 1

        less_than_chi = Chi_squared[Vel < radial_vel[j]]

        more_than_chi = Chi_squared[Vel > radial_vel[j]]

        less_than_vel = Vel[Vel < radial_vel[j]]

        more_than_vel = Vel[Vel > radial_vel[j]]

        error_rad_vel_less[j] = np.max(less_than_vel[less_than_chi > Chi_plus_one]) - radial_vel[j]

        error_rad_vel_more[j] = radial_vel[j] - np.min(more_than_vel[more_than_chi > Chi_plus_one])

# Shrink current axis by 20%
box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

Average_sig = (error_rad_vel_more + error_rad_vel_less) / 2

# Creating our Hessian Matrix and solving to find our best fit parameters

w = 1 / Average_sig ** 2   # weight

W = np.sum(w)

Wsin = np.sum(w * np.sin(2 * np.pi * phase))

Wcos = np.sum(w * np.cos(2 * np.pi * phase))

Wsin2 = np.sum(w * (np.sin(2 * np.pi * phase)) ** 2)

Wcos2 = np.sum(w * (np.cos(2 * np.pi * phase)) ** 2)

Wcossin = np.sum(w * np.cos(2 * np.pi * phase) * np.sin(2 * np.pi * phase))


Hessian_Matrix = np.array([W, Wsin, Wcos, Wsin, Wsin2, Wcossin, Wcos, Wcossin, Wcos2]).reshape(3, 3)

vel = unumpy.uarray(radial_vel[:], np.abs(Average_sig))

Wv = np.sum(w * vel)

Wvsin = np.sum(w * np.sin(2 * np.pi * phase) * vel)

Wvcos = np.sum(w * np.cos(2 * np.pi * phase) * vel)

Vel_Matrix = np.array([Wv, Wvsin, Wvcos]).reshape(3, 1)

Inverse_Hessian = np.linalg.inv(Hessian_Matrix)

Dot_Product = np.dot(Inverse_Hessian, Vel_Matrix)


# Determining our curve fit parameters

gamma = Dot_Product[0] / (1e3)

Kx = Dot_Product[1] / (1e3)

Ky = Dot_Product[2] / (1e3)

print "gamma : ", gamma

print "Kx : ", Kx

print "Ky : ", Ky


# Plotting radial velocity and phase

phase_trial = np.linspace(np.min(phase), np.max(phase), 100)

sinusoidal_vel = sinusoidal_curve(gamma, Kx, Ky, phase_trial)

K = (Kx ** 2 + Ky ** 2) ** (1 / 2)

ls = 'None'

fig, ax = plt.subplots()

# standard error bars
ax.errorbar(phase, radial_vel[:] / (1e3),
            yerr=[error_rad_vel_less[:] / (1e3), error_rad_vel_more[:] / (1e3)],
            marker=".", linestyle=ls, capsize=2)

ax.set_ylabel(r"$Velocity \/ (Km/s)$")
ax.set_xlabel(r"$binary \/ phase$")

Curve_val = unumpy.nominal_values(sinusoidal_vel)

ax.plot(phase_trial, Curve_val)


K = (Kx ** 2 + Ky ** 2) ** (1 / 2)

K_vector = u.ufloat(unumpy.nominal_values(K), unumpy.std_devs(K))

# Defining our parameters

P = 0.3440915  # units of days

M_sun = 1.989e30


# Mass function

f_Mx = f_Mx_func(P, K)


# Finding minimum mass using i = 90 degrees and Mc = 0.7

incline = np.deg2rad(90)

Mc = 0.7*M_sun

Mx_i = Mx_func(Mc, incline, f_Mx, 1)

Mx_i_sol_unit = Mx_i/M_sun


print 'K : {:.2f}'.format(K_vector)

f_Mx_sol_unit = f_Mx/M_sun

f_Mx_sol_unit_vector = u.ufloat(unumpy.nominal_values(f_Mx_sol_unit), unumpy.std_devs(f_Mx_sol_unit))

print 'f(Mx) : {:.2f}'.format(f_Mx_sol_unit_vector)

print 'Mx : {:.2f}'.format(Mx_i_sol_unit)


""" Monte-carlo eror propagation """

# Number of monte carlo trials

Ntrial = 1e5

# Creating standard deviations for each parameter

standard_Mc = 0.4*Mc
standard_K = 0.1*K_vector.nominal_value
standard_P = 0.1*P

angle = np.linspace(np.deg2rad(70), np.deg2rad(90), 1e4)

# Creating arrays for paramaters and Mx

incline_trial = np.zeros(np.int(Ntrial))
Mc_trial = np.zeros(np.int(Ntrial))
K_trial = np.zeros(np.int(Ntrial))
P_trial = np.zeros(np.int(Ntrial))

f_Mx_trial = np.zeros(np.int(Ntrial))
Mx_trial = np.zeros(np.int(Ntrial))

for l in xrange(np.int(Ntrial)):

    incline_trial[l] = np.random.choice(angle)

    Mc_trial[l] = np.random.normal(Mc, standard_Mc)

    K_trial[l] = np.random.normal(unumpy.nominal_values(K), standard_K)

    P_trial[l] = np.random.normal(P, standard_P)

    f_Mx_trial[l] = f_Mx_func(P_trial[l], K_trial[l])

    Mx_trial[l] = Mx_func(Mc_trial[l], incline_trial[l], f_Mx_trial[l], 0)


Mx_sol_unit = Mx_trial / M_sun

Mx_mean = np.mean(Mx_sol_unit)

Mx_median = np.median(Mx_sol_unit)

Mx_std = np.std(Mx_sol_unit)

Mx = u.ufloat(Mx_mean, Mx_std)


print 'Monte-Carlo Mx : {:.2f}'.format(Mx)


# Plotting probabilities

hist = np.histogram(Mx_sol_unit, bins=30)
hist_dist = sps.rv_histogram(hist)

X = np.linspace(0.0, 15, 500)

plt.figure()

plt.plot(X, hist_dist.pdf(X), label='PDF')
plt.plot(X, hist_dist.cdf(X), label='CDF')

xcoords = [Mx_mean, Mx_median]

text_coords = [Mx_mean + 0.1, Mx_median - 0.7]

plot_text = ['mean', 'median']

colour = ['red', 'red']

plot_index = 0

for xc in xcoords:

    plt.axvline(x=xc, color=colour[plot_index], linestyle='--')

    plt.text(text_coords[plot_index], 1, plot_text[plot_index], rotation=90)

    plot_index += 1

plt.axvspan((Mx.nominal_value - Mx.std_dev), (Mx.nominal_value + Mx.std_dev), facecolor='g', alpha=0.5)

plt.xlabel(r"$Mass \/ of \/ compact \/ object \/ \/ ( \/ M_\odot \/ )$")
plt.ylabel(r"$Probability$")


print "Probability > 3 Msol : ", 1 - hist_dist.cdf(3)

print "Mean of Mx : ", Mx_mean

print "Median of Mx : ", Mx_median


print "1-sigma : ", (hist_dist.cdf((Mx.nominal_value + Mx.std_dev))) - (hist_dist.cdf((Mx.nominal_value - Mx.std_dev)))

plt.show()