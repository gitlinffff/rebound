import numpy as np

def r2e(rv, mu):
    # Unpack position and velocity vectors
    r = rv[:3]
    v = rv[3:]

    # Calculate specific angular momentum
    h = np.cross(r, v)
    h_norm = np.linalg.norm(h)
    
    # Inclination
    i = np.arccos(h[2] / h_norm)

    # Eccentricity vector
    ev = np.cross(v, h) / mu - r / np.linalg.norm(r)
    e = np.linalg.norm(ev)

    # Semi-major axis
    a = h_norm**2 / (mu * (1 - e**2))
    n = np.sqrt(mu / np.abs(a)**3)

    coe = [a, e, i]

    return coe

# Supporting functions
def f2E(xsit, ecc):
    if xsit >= 0:
        k = int(0.25 * xsit / np.pi)
    else:
        k = int(np.floor(0.25 * xsit / np.pi))
    m = 0.5 * xsit - 2 * k * np.pi

    if m >= 0 and m < 0.5 * np.pi:
        de = 2 * np.arccos(np.sqrt(((ecc + np.cos(xsit)) / (1 + ecc * np.cos(xsit))) + 1) * 0.5) + 4 * k * np.pi
    elif m >= 0.5 * np.pi and m < np.pi:
        de = 2 * np.arccos(-np.sqrt(((ecc + np.cos(xsit)) / (1 + ecc * np.cos(xsit))) + 1) * 0.5) + 4 * k * np.pi
    elif m >= np.pi and m < 1.5 * np.pi:
        de = 4 * np.pi - 2 * np.arccos(-np.sqrt(((ecc + np.cos(xsit)) / (1 + ecc * np.cos(xsit))) + 1) * 0.5) + 4 * k * np.pi
    elif m >= 1.5 * np.pi and m < 2 * np.pi:
        de = 4 * np.pi - 2 * np.arccos(np.sqrt(((ecc + np.cos(xsit)) / (1 + ecc * np.cos(xsit))) + 1) * 0.5) + 4 * k * np.pi

    return de

def E2M(E, e):
    return E - e * np.sin(E)

def M2E(M, e):
    E = M
    fe = 1
    while abs(fe) > 1e-15:
        fe = E - e * np.sin(E) - M
        dfe = 1 - e * np.cos(E)
        E = E - fe / dfe
    return E

def E2f(dde, ecc):
    if dde >= 0:
        k = int(0.25 * dde / np.pi)
    else:
        k = int(np.floor(0.25 * dde / np.pi))
    m = 0.5 * dde - 2 * k * np.pi

    if m >= 0 and m < 0.5 * np.pi:
        true_anomaly = 2 * np.arccos(np.sqrt(((np.cos(dde) - ecc) / (1 - ecc * np.cos(dde))) + 1) * 0.5) + 4 * k * np.pi
    elif m >= 0.5 * np.pi and m < np.pi:
        true_anomaly = 2 * np.arccos(-np.sqrt(((np.cos(dde) - ecc) / (1 - ecc * np.cos(dde))) + 1) * 0.5) + 4 * k * np.pi
    elif m >= np.pi and m < 1.5 * np.pi:
        true_anomaly = 4 * np.pi - 2 * np.arccos(-np.sqrt(((np.cos(dde) - ecc) / (1 - ecc * np.cos(dde))) + 1) * 0.5) + 4 * k * np.pi
    elif m >= 1.5 * np.pi and m < 2 * np.pi:
        true_anomaly = 4 * np.pi - 2 * np.arccos(np.sqrt(((np.cos(dde) - ecc) / (1 - ecc * np.cos(dde))) + 1) * 0.5) + 4 * k * np.pi

    return true_anomaly

