import numpy as np
from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from r2e_1 import r2e

#h_norm = np.linalg.norm(h)
r = np.array([-6045, -3490, 2500]) * 1000.
v = np.array([-3.457, 6.618, 2.533]) * 1000.

G = 6.6743015e-11
M_host = Earth.mass.value 

rv = np.concatenate((r,v))
mu = G * M_host

coe = r2e(rv,mu)
#orb = Orbit.from_vectors(Earth, r, v)
#print(dir(Earth))
print(coe)



r = [-6045, -3490, 2500] << u.km
v = [-3.457, 6.618, 2.533] << u.km / u.s

orb = Orbit.from_vectors(Earth, r, v)
print(orb)

