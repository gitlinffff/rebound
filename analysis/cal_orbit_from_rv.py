import numpy as np
from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit

#h_norm = np.linalg.norm(h)
#r = np.array([-6045, -3490, 2500]) * 1000.
#v = np.array([-3.457, 6.618, 2.533]) * 1000.

#G = 6.6743015e-11
#M_host = Earth.mass.value 

#rv = np.concatenate((r,v))
#mu = G * M_host

#coe = r2e(rv,mu)
#orb = Orbit.from_vectors(Earth, r, v)
#print(dir(Earth))
#print(coe)


# DART impact (2022-Sep-26 23:14:24.1830 UTC)
#r = [1.556582267294774E+11, 1.349129152256939E+10, -8.638156994081538E+09] << u.m
#v = [-7.322785629453647E+03, 3.319798419238497E+04, 9.918308846239352E+02] << u.m / u.s

# 160s after the impact (2022-Sep-26 23:17:04.1830 UTC)
r = [1.556570550147468E+11, 1.349660319404291E+10, -8.637998297279608E+09] << u.m
v = [-7.323648504458615E+03, 3.319790922163674E+04, 9.918791394754933E+02] << u.m / u.s
orbit = Orbit.from_vectors(Sun, r, v)

# Print detailed orbital elements
print(f"Semi-major axis (a): {orbit.a.to(u.AU)}")
print(f"Eccentricity (e): {orbit.ecc}")
print(f"Inclination (i): {orbit.inc.to(u.deg)}")
print(f"Argument of Periapsis (ω): {orbit.argp.to(u.deg)}")
print(f"Longitude of Ascending Node (Ω): {orbit.raan.to(u.deg)}")
print(f"True Anomaly (ν): {orbit.nu.to(u.deg)}")
print(f"Epoch: {orbit.epoch}")

