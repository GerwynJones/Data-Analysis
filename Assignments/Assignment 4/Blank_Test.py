# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import uncertainties as u
import numpy as np


@u.wrap
def f(n=0, *P):
    ''' compute the nth root of the polynomial P and the uncertainty of the root'''

    p = np.array(P)
    N = len(p)

    M = np.diag(np.ones((N - 2,), p.dtype), -1)

    M[0, :] = -p[1:] / p[0]

    r = np.linalg.eigvals(M)

    # r.sort()  # there is no telling what order the values come out in

    return r[n]


# our polynomial coefficients and standard errors
c, b, a, d = [-0.99526746, 0.011546, 10.00188999, 22.033]
sc, sb, sa, sd = [2, 1, 0, 1]

A = u.ufloat((a, sa))
B = u.ufloat((b, sb))
C = u.ufloat((c, sc))
D = u.ufloat((d, sd))

Coeff = [A, B, C, D]

Coeff_real = [a, b, c, d]

Root_real = np.roots(Coeff_real)

Root_very_real = Root_real[Root_real.imag == 0].real



for m in xrange(len(Coeff) - 1):

    Root = f(m, *Coeff)

    if Root.nominal_value == Root_very_real:

        print Root

        Real = Root

    else:

        continue

