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

def norm(x):
    return np.sqrt(np.dot(x, x))

def Px(φ):
    return np.array([[1, 0, 0],
                     [0, np.cos(φ), -np.sin(φ)],
                     [0, np.sin(φ), np.cos(φ)]])

def Pz(φ):
    return np.array([[np.cos(φ), -np.sin(φ), 0],
                     [np.sin(φ), np.cos(φ), 0],
                     [0, 0, 1]])

def r_ellipse(a, e, f):
    return a * (1 - e**2) / (1 + e * np.cos(f))
    
def carttoels(x, y, z, vx, vy, vz, μ=1):
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

def barytohelio():
    pass

def heliotobary():
    pass


