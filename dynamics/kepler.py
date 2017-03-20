import numpy as np


"""
Solving Kepler's equation.

M = E - e sin(E)

M : mean anomaly, equal to n(t - τ), where n is the mean motion, t is the time,
    and τ is the time of pericenter passage
E : eccentric anomaly, angle onto the circumscribed circle such that the point
    on the circumscribed circle and the point on the ellipse with true anomaly
    f share the same y-coordinate in the standard frame
"""

def zeroth_function(M, e):
    """Construct the function to zero.

    Parameters
    ----------
    M : mean anomaly
    e : eccentricity

    Returns
    -------
    f : function of E, the eccentric anomaly
    """
    return lambda E: E - e * np.sin(E) - M

def first_function(e):
    """Construct the first derivative of the function to zero.

    Parameters
    ----------
    e : eccentricity

    Returns
    -------
    f : function of E, the eccentric anomaly
    """
    return lambda E: 1 - e * np.cos(E)


def second_function(e):
    """Construct the second derivative of the function to zero.

    Parameters
    ----------
    e : eccentricity

    Returns
    -------
    f : function of E, the eccentric anomaly
    """
    return lambda E: e * np.sin(E)


def third_function(e):
    """Construct the third derivative of the function to zero.

    Parameters
    ----------
    e : eccentricity

    Returns
    -------
    f : function of E, the eccentric anomaly
    """
    return lambda E: -e * np.cos(E)


def initial_guess(M, e, k=0.85):
    """Initial guess for E, with k = 0.85 from Danby (1988).

    Parameters
    ----------
    M : mean anomaly
    e : eccentricity
    k : optional, numeric factor, defaults to 0.85

    Returns
    -------
    E0 : initial guess for E
    """
    return M + np.sign(np.sin(M)) * k * e


def iterate(E, func0, func1, func2, func3):
    """"Produce the next value of E.
    
    Parameters
    ----------
    E : eccentric anomaly
    func0 : function to zero
    funcN : N-th derivative of function
    """
    f0 = func0(E)
    f1 = func1(E)
    f2 = func2(E)
    f3 = func3(E)
    d1 = -f0 / f1
    d2 = -f0 / (f1 + 0.5 * d1 * f2)
    d3 = -f0 / (f1 + 0.5 * d2 * f2 + 1/6 * d2 ** 2 * f3)
    return E + d3


def solve_kepler(M, e, rtol=1e-6):
    """Solve Kepler's equation numerically.
    
    Parameters
    ----------
    M : mean anomaly
    e : eccentricity
    rtol : optional, relative tolerance threshold to end iteration

    Returns
    -------
    E : eccentric anomaly
    """
    rerr = np.inf
    func0 = zeroth_function(M, e)
    func1 = first_function(e)
    func2 = second_function(e)
    func3 = third_function(e)
    E0 = initial_guess(M, e)
    i = 0
    while rerr > rtol:
        E = iterate(E0, func0, func1, func2, func3)
        rerr = np.abs((E - E0) / E0)
        E0 = E
        i += 1
        if i > 100:
            raise RuntimeError("Number of iterations exceeded 100!")
    return E
        
def kepler_orbit(t, a, e, μ=1, cart=False):
    """Solve for the radii and true anomalies of a Keplerian orbit at the specified timesteps.
    t : array of times in years / 2π
    a : semi-major axis in AU
    e : eccentricity
    μ : standard gravitational parameter, μ = G(M1 + M2)
    cart : bool, if true return x, y coords
    """
    P = 1 / np.sqrt(μ) * a ** (3/2)          # period
    n = 1 / P                                # mean motion
    M = n * t                                # mean anomaly
    try:
        E = np.array([solve_kepler(MM, e) for MM in M])
    except TypeError:
        E = solve_kepler(M, e)
    r = a * (1 - e * np.cos(E))
    f = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
    if not cart:
        return r, f
    return r * np.cos(f), r * np.sin(f)
    
