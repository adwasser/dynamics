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

def Px(phi):
    return np.array([[1, 0, 0],
                     [0, np.cos(phi), -np.sin(phi)],
                     [0, np.sin(phi), np.cos(phi)]])
def Pz(phi):
    return np.array([[np.cos(phi), -np.sin(phi), 0],
                     [np.sin(phi), np.cos(phi), 0],
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
    # Ω = np.arcsin(hx / (norm(h) * np.sin(i)))
    f = np.arctan2(norm(rdot) * a * (1 - e**2) / (norm(h) * e),
                   (a * (1 - e**2) / R - 1) / e)
    # f = np.arcsin(norm(rdot) * a * (1 - e**2) / (norm(h) * e))
    ω = np.arctan2(z / (R * np.sin(i)),
                   np.cos(Ω)**-1 * (x / R + np.sin(Ω) * np.cos(i) * z / (R * np.sin(i)))) - f
    # ω = np.arcsin(z / (R * np.sin(i))) - f
    return a, e, i, ω, Ω, f

def elstocart(a, e, i, ω, Ω, f):
    rot = Pz(Ω) @ Px(i) @ Pz(ω)
    r = r_ellipse(a, e, f) * np.array([np.cos(f), np.sin(f)])
    x, y, z = rot @ r
    pass

def barytohelio():
    pass

def heliotobary():
    pass


