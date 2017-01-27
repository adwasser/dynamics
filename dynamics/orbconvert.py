import numpy as np

'''
All angles in radians.

αβγδεζηθικλμνξοπρστυφχψω
ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ

a : semi-major axis
e : eccentricity
i : inclination
ω : argument of periapsis
Ω : longitude of the ascending node
f : true anomaly

x, y, z : cartesian position
vx, vy, vz : cartesian velocity
'''

π = np.pi

def norm(x):
    """Calculates the euclidean norm of the vector."""
    return np.sqrt(np.dot(x, x))


def Px(φ):
    """Returns the matrix for rotation about the x-axis by angle φ."""
    return np.array([[1, 0, 0],
                     [0, np.cos(φ), -np.sin(φ)],
                     [0, np.sin(φ), np.cos(φ)]])


def Py(φ):
    """Returns the matrix for rotation about the y-axis by angle φ."""
    return np.array([[np.cos(φ), 0, -np.sin(φ)],
                     [0, 1, 0],
                     [np.sin(φ), 0, np.cos(φ)]])


def Pz(φ):
    """Returns the matrix for rotation about the z-axis by angle φ."""    
    return np.array([[np.cos(φ), -np.sin(φ), 0],
                     [np.sin(φ), np.cos(φ), 0],
                     [0, 0, 1]])


def r_ellipse(a, e, f):
    """Calculates the radius frome the foci of an ellipse.

    Parameters
    ----------
    a : semi-major axis
    e : ellipticity
    f : true anomaly

    Returns
    -------
    r : radius between focus and the ellipse at that true anomaly
    """
    return a * (1 - e**2) / (1 + e * np.cos(f))
    

def carttoels(x, y, z, vx, vy, vz, μ=1):
    """ Transforms cartesian coordinates to orbital elements.

    Parameters
    ----------
    x, y, z: cartesian position coordinates, in barycenter frame
    vx, vy, vz: cartesian velocity coordinates, in inertial frame
    μ : optional, standard gravitational parameter, G(m1 + m2), defaults to 1

    Returns
    -------
    a : semi-major axis
    e : eccentricity
    i : inclination
    ω : argument of periapsis
    Ω : longitude of the ascending node
    f : true anomaly    
    """
    r = np.array([x, y, z])
    rdot = np.array([vx, vy, vz])
    R = norm(r)
    v2 = norm(rdot) ** 2
    h = np.cross(r, rdot)
    hx, hy, hz = h
    a = (2 / R - v2 / μ) ** -1
    e = np.sqrt(1 - norm(h) ** 2 / (μ * a))
    i = np.arccos(hz / norm(h))
    Ω = np.arctan2(-hx, hy)
    f = np.arctan2(norm(rdot) * a * (1 - e**2) / (norm(h) * e),
                   (a * (1 - e**2) / R - 1) / e)
    ω = np.arctan2(z / (R * np.sin(i)),
                   np.cos(Ω)**-1 * (x / R + np.sin(Ω) * np.cos(i) * z / (R * np.sin(i)))) - f
    return a, e, i, ω, Ω, f


def elstocart(a, e, i, ω, Ω, f, μ=1):
    """ Transforms orbital elements to cartesian coordinates.

    Parameters
    ----------
    a : semi-major axis
    e : eccentricity
    i : inclination
    ω : argument of periapsis
    Ω : longitude of the ascending node
    f : true anomaly
    μ : optional, standard gravitational parameter, G(m1 + m2), defaults to 1

    Returns
    -------
    x, y, z: cartesian position coordinates, in barycenter frame
    vx, vy, vz: cartesian velocity coordinates, in inertial frame    
    """    
    rot = Pz(Ω) @ Px(i) @ Pz(ω)
    r = r_ellipse(a, e, f) * np.array([np.cos(f), np.sin(f), 0])
    R = norm(r)
    x, y, z = rot @ r
    h = np.sqrt(μ * a * (1 - e**2))
    fdot = h / R ** 2
    Rdot = a * (1 - e ** 2) * e * fdot * np.sin(f) / (1 + e ** 2 * np.cos(f)) ** 2
    vx_plane = Rdot * np.cos(f) - R * fdot * np.sin(f)
    vy_plane = Rdot * np.sin(f) + R * fdot * np.cos(f)
    v_plane = np.array([vx_plane, vy_plane, 0])
    vx, vy, vz = rot @ v_plane
    return x, y, z, vx, vy, vz


def barytohelio(x_bary, y_bary, z_bary, vx_bary, vy_bary, vz_bary, mass_ratio, μ=1):
    """Convert barycentric coordinates to heliocentric coordinates.

    Parameters
    ----------
    x_bary, y_bary, z_bary : barycentric cartesian position coordinates
    vx_bary, vy_bary, vz_bary : barycentric cartesian velocity coordinates
    mass_ration : mass ratio of minor body to Sun
    μ : optional, standard gravitational parameter, G(m1 + m2), defaults to 1
    
    Returns
    -------
    x_helio, y_helio, z_helio : heliocentric cartesian position coordinates
    vx_helio, vy_helio, vz_helio : heliocentric cartesian velocity coordinates
    """
    r2 = np.array([x_bary, y_bary, z_bary])
    v2 = np.array([vx_bary, vy_bary, vz_bary])
    
    R2 = norm(r2)
    R = (mass_ratio + 1) * R2
    R1 = mass_ratio / (mass_ratio + 1) * R
    r_helio = r2 + R1 * r2 / R2

    a2, e, i, ω, Ω, f = carttoels(x_bary, y_bary, z_bary, vx_bary, vy_bary, vz_bary, μ=μ)
    a = (mass_ratio + 1) * a2
    a1 = mass_ratio / (mass_ratio + 1) * a
    x1, y1, z1, vx1, vy1, vz1 = elstocart(a1, e, i, ω + π, Ω, f, μ=μ)
    r1 = np.array([x1, y1, z1])
    v1 = np.array([vx1, vy1, vz1])
    v_helio = v2 - v1
    assert np.all(np.isclose(r_helio, r2 - r1))
    
    x_helio, y_helio, z_helio = r_helio
    vx_helio, vy_helio, vz_helio = v_helio
    return x_helio, y_helio, z_helio, vx_helio, vy_helio, vz_helio


def heliotobary(x_helio, y_helio, z_helio, vx_helio, vy_helio, vz_helio, mass_ratio, μ=1):
    """Convert heliocentric coordinates to barycentric coordinates.

    Parameters
    ----------
    x_helio, y_helio, z_helio : heliocentric cartesian position coordinates
    vx_helio, vy_helio, vz_helio : heliocentric cartesian velocity coordinates
    mass_ration : mass ratio of minor body to Sun
    μ : optional, standard gravitational parameter, G(m1 + m2), defaults to 1
    
    Returns
    -------
    x_bary, y_bary, z_bary : barycentric cartesian position coordinates
    vx_bary, vy_bary, vz_bary : barycentric cartesian velocity coordinates
    """
    pass


