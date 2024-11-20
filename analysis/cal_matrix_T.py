# calculate the matrix T that converts a vector from Sun body center frame to the Didymos system barycenter frame
# present this process in the slides
import numpy as np

""" 2022-Sep-26 23:14:24.1830 UTC (moment of the impact)
 *  (https://ssd.jpl.nasa.gov/horizons/app.html#/) 
 *  Coordinate Center: Sun (body center) [500@10] """
r_DSB_t0 = np.array([1.556582267294774E+11, 1.349129152256939E+10, -8.638156994081538E+09])# position of Didymos System Barycenter
v_DSB_t0 = np.array([-7.322785629453647E+03, 3.319798419238497E+04, 9.918308846239352E+02]) # velocity of Didymos System Barycenter
r_Dimor_t0 = np.array([1.556582259038835E+11, 1.349129068894670E+10, -8.638156976265389E+09])
v_Dimor_t0 = np.array([-7.322906365387563E+03, 3.319810316934613E+04, 9.918029911906725E+02])

""" 2022-Sep-26 23:17:04.1830 UTC (160s after the impact)
 *  (https://ssd.jpl.nasa.gov/horizons/app.html#/) 
 *  Coordinate Center: Sun (body center) [500@10] """
r_DSB_t1 = np.array([1.556570550147468E+11, 1.349660319404291E+10, -8.637998297279608E+09]) # position of Didymos System Barycenter
v_DSB_t1 = np.array([-7.323648504458615E+03, 3.319790922163674E+04, 9.918791394754933E+02]) #velocity of Didymos System Barycenter
r_Dimor_t1 = np.array([1.556570541703198E+11, 1.349660237935313E+10, -8.637998283861816E+09])
v_Dimor_t1 = np.array([-7.323764773455447E+03, 3.319802896055199E+04, 9.918516242373627E+02])

# calculate r,v,p vector of Dimorphos relative to Didymos
r_Dimor_rel = r_Dimor_t0 - r_DSB_t0
v_Dimor_rel = v_Dimor_t0 - v_DSB_t0
p_Dimor_rel = np.cross(r_Dimor_rel, v_Dimor_rel)

T11_T12_T13 = (r_Dimor_rel) / np.linalg.norm(r_Dimor_rel)
T31_T32_T33 = -p_Dimor_rel / np.linalg.norm(p_Dimor_rel)
T21_T22_T23 = np.cross(T31_T32_T33, T11_T12_T13)

print(np.linalg.norm(T11_T12_T13))
print(np.linalg.norm(T21_T22_T23))
print(np.linalg.norm(T31_T32_T33))

a = np.dot(T11_T12_T13, T21_T22_T23)
b = np.dot(T11_T12_T13, T31_T32_T33)
c = np.dot(T31_T32_T33, T21_T22_T23)
print(a,b,c)

np.set_printoptions(precision=15)
T = np.vstack([T11_T12_T13,T21_T22_T23,T31_T32_T33])
print('# Matrix T:\n',T)
#print(np.arccos(T31_T32_T33[-1])/np.pi*180+169.3)
