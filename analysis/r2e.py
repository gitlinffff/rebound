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

    # Right ascension of the ascending node (Ω)
    if abs(np.sin(i)) < 1e-12:
        omg = 0
    else:
        omg = np.arctan2(-h[1], h[0])

    if h[1] > 0:
        omg += np.pi

    if omg < 0:
        omg += 2 * np.pi

    # Argument of periapsis (ω) and true anomaly (f)
    if abs(np.sin(i)) < 1e-12 and e > 1e-12:
        w = np.arctan2(ev[1], ev[0])
        u = np.arctan2(r[1], r[0])
        if ev[0] < 0:
            w += np.pi
        if r[0] < 0:
            u += np.pi

    elif abs(np.sin(i)) > 1e-12 and e < 1e-12:
        w = 0
        u = np.arctan2(r[2], (r[1] * np.sin(omg) + r[0] * np.cos(omg)) * np.sin(i))
        if r[2] * u < 0:
            u += np.pi
        elif r[2] < 0 and u < 0:
            u += 2 * np.pi

    elif abs(np.sin(i)) < 1e-12 and e < 1e-12:
        w = 0
        u = np.arctan2(r[1], r[0])
        if r[0] < 0:
            u += np.pi

    elif abs(np.sin(i)) > 1e-12 and e > 1e-12:
        w = np.arctan2(ev[2], (ev[1] * np.sin(omg) + ev[0] * np.cos(omg)) / np.sin(i))
        u = np.arctan2(r[3], (r[1] * np.sin(omg) + r[0] * np.cos(omg)) * np.sin(i))
        if ev[2] * w < 0:
            w += np.pi
        elif ev[2] < 0 and w < 0:
            w += 2 * np.pi

        if r[2] * u < 0:
            u += np.pi
        elif r[2] < 0 and u < 0:
            u += 2 * np.pi

    if e < 1:
        flag = 0
        f = u - w
        f = f % (2 * np.pi)

        coe = [a, e, omg, i, w, f]
    else:
        flag = 1
        f = np.arccos((-a * (e**2 - 1) / np.linalg.norm(r) - 1) / e)
        rp = r + v * 1e-7
        if np.linalg.norm(rp) < np.linalg.norm(r):
            f = -f
        H = 2 * np.arctanh(np.sqrt((e - 1) / (1 + e)) * np.tan(f / 2))
        tp = -(e * np.sinh(H) - H) / n
        coe = [a, e, omg, i, w, f, tp]

    return coe, flag

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

