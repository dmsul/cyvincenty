import numpy as np
cimport numpy as np

import cython

from cython.parallel import prange

from libc.math cimport sin, cos, tan, atan, sqrt, atan2, pi

# WGS 84
cdef double  a = 6378137  # meters
cdef double  f = 1 / 298.257223563
cdef double  b = 6356752.314245  # meters; b = (1 - f)a

cdef int MAX_ITERATIONS = 200
cdef double  CONVERGENCE_THRESHOLD = 1e-12  # .000,000,000,001


@cython.wraparound(False)
@cython.boundscheck(False)
def vincenty_cross(
        np.ndarray[np.float32_t] ax,
        np.ndarray[np.float32_t] ay,
        np.ndarray[np.float32_t] bx,
        np.ndarray[np.float32_t] by):
    """
    Calculates vincenty distance for each pair of longitude/latitude points in
    vectors a (ax, bx) and b (bx, by) using vincenty method. Uses WGS84.

    Arguments must be 1-D numpy arrays (float32).

    Returns numpy array of shape (len(a), len(b)).
    """
    cdef int I = ax.shape[0]
    cdef int J = bx.shape[0]
    cdef np.ndarray[np.float32_t, ndim=2] out = np.zeros((I, J), dtype=np.float32)
    cdef int i, j

    cdef float this_ax, this_ay, this_bx, this_by

    try:
        assert I == ay.shape[0]
        assert J == by.shape[0]
    except AssertionError:
        raise ValueError("Input x/y vectors must be same length.")

    for i in prange(I, nogil=True):
        this_ax = ax[i]
        this_ay = ay[i]
        for j in range(J):
            this_bx = bx[j]
            this_by = by[j]
            out[i, j] = vincenty(this_ax, this_ay, this_bx, this_by)

    return out


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef float vincenty(double ax, double ay, double bx, double by) nogil:
    """
    Calculates the distance in kilometers between longitude/latitude points a
    (ax, bx) and b (bx, by) using vincenty method. Uses WGS84.
    """

    cdef int iteration
    cdef double U1
    cdef double U2
    cdef double L
    cdef double Lambda
    cdef double sinU1
    cdef double sinU2
    cdef double cosU1
    cdef double cosU2

    cdef double sinLambda
    cdef double cosLambda
    cdef double sinSigma
    cdef double cosSigma
    cdef double sigma
    cdef double sinAlpha
    cdef double cosSqAlpha
    cdef double cos2SigmaM
    cdef double C
    cdef double LambdaPrev
    cdef double lambda_diff

    cdef double uSq
    cdef double A
    cdef double B
    cdef double deltaSigma
    cdef double s

    # short-circuit coincident points
    if ay == by and ax == bx:
        return 0.0

    U1 = atan((1 - f) * tan(radians(ay)))
    U2 = atan((1 - f) * tan(radians(by)))
    L = radians(bx - ax)
    Lambda = L

    sinU1 = sin(U1)
    cosU1 = cos(U1)
    sinU2 = sin(U2)
    cosU2 = cos(U2)

    for iteration in range(MAX_ITERATIONS):
        sinLambda = sin(Lambda)
        cosLambda = cos(Lambda)
        sinSigma = sqrt((cosU2 * sinLambda) ** 2 +
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ** 2)
        if sinSigma == 0:
            return 0.0  # coincident points
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
        sigma = atan2(sinSigma, cosSigma)
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
        cosSqAlpha = 1 - sinAlpha ** 2
        if cosSqAlpha == 0:
            cos2SigmaM = 0
        else:
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
        LambdaPrev = Lambda
        Lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma *
                                               (cos2SigmaM + C * cosSigma *
                                                (-1 + 2 * cos2SigmaM ** 2)))
        lambda_diff = Lambda - LambdaPrev
        if lambda_diff < 0:
            lambda_diff = -1 * lambda_diff
        if lambda_diff < CONVERGENCE_THRESHOLD:
            break  # successful convergence
    else:
        return -1  # failure to converge

    uSq = cosSqAlpha * (a ** 2 - b ** 2) / (b ** 2)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *
                 (-1 + 2 * cos2SigmaM ** 2) - B / 6 * cos2SigmaM *
                 (-3 + 4 * sinSigma ** 2) * (-3 + 4 * cos2SigmaM ** 2)))
    s = b * A * (sigma - deltaSigma)

    s /= 1000  # meters to kilometers

    return s


@cython.cdivision(True)
cdef double radians(double x) nogil:
    return x * pi / 180
