import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import healpy as hp
from ligo.skymap import distance
import numpy as np
from database_methods import *
import csv
from mpl_toolkits.basemap import Basemap
from astropy.io import fits
from ligo.skymap.postprocess import find_greedy_credible_levels
from datetime import datetime
from ligo.skymap import distance
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import scipy.stats
import random

c = 299792.458
H0_min = 0
H0_max = 150
num = 1000
script_start = datetime.now()

### GRAVITATIONAL WAVE
prob, distmu, distsigma, distnorm = hp.read_map("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0",field=range(4))
dist_mean, dist_std, norm = distance.parameters_to_moments(distmu, distsigma)
npix = len(prob)
nside = hp.npix2nside(npix)
gw_bools = np.zeros(npix, dtype=bool)
indexes = np.linspace(0, npix - 1, npix, dtype=int)
credible_levels = find_greedy_credible_levels(prob)

perc_now = 0
start = datetime.now()
print("Loading GW ", str(start.time()))
for i in indexes:
    if i / npix >= perc_now / 100:
        print(str(perc_now) + "%")
        perc_now = perc_now + 10
    if credible_levels[i] <= 0.90:
        gw_bools[i] = True
print("100%")

alphas = 1 - credible_levels[indexes[gw_bools]]
theta, phi = hp.pix2ang(nside, indexes[gw_bools])
gw_ra = [np.rad2deg(x) for x in phi]
gw_dec = [np.rad2deg(0.5 * np.pi - x) for x in theta]

print("Len  GW =", len(gw_ra))
print("Time to Complete: " + str(datetime.now() - start))

# SHOW GW RA AND DEC LIMITS
print("First Blob")
print("RA Limits = [" + str(min([x for x in gw_ra if x <= 20])) + ", " + str(max([x for x in gw_ra if x <= 20])) + "]")
print("Dec Limits = [" + str(min([x for x in gw_dec if x >= -30])) + ", " + str(max([x for x in gw_dec if x >= -30])) + "]")
print("Second Blob")
print("RA Limits = [" + str(min([x for x in gw_ra if x >= 20])) + ", " + str(max([x for x in gw_ra if x >= 20])) + "]")
print("Dec Limits = [" + str(min([x for x in gw_dec if x <= -30])) + ", " + str(max([x for x in gw_dec if x <= -30])) + "]")


### GALAXIES
### PANSTARRS1 - DB
start = datetime.now()
print("Loading PS1 Galaxy " + str(start.time()))
db_query = '''
SELECT raMean, decMean, z_phot, class, objID, ps_score, z_photErr FROM PS1_Galaxy_Final_PS;
'''
PS1 = query_db([db_query])[0]
PS1_ra = [x[0] for x in PS1]
PS1_dec = [x[1] for x in PS1]
PS1_z = [x[2] for x in PS1]
PS1_class = [x[3] for x in PS1]
PS1_objid = [x[4] for x in PS1]
PS1_ps_score = [x[5] for x in PS1]
PS1_z_err = [x[6] for x in PS1]
print("Those with Galaxy Tag: " + str(len([x for x in PS1_class if x == "GALAXY"])/len(PS1_class) * 100) + "%")
print("Those with Quasar Tag: " + str(len([x for x in PS1_class if x == "QSO"]) / len(PS1_class) * 100) + "%")
print("Len PS1 = " + str(len(PS1_ra)))
print("Time to Complete: " + str(datetime.now() - start))


### GLADE
start = datetime.now()
print("Loading GLADE Galaxy " + str(start.time()))
db_query = '''
SELECT RA, _DEC, z, B, z_err FROM GLADE_Galaxy_Final 
WHERE 
    (RA >= 10.0 AND
    RA <= 15.0 AND
    _DEC <= -22.0 AND
    _DEC >= -28.0)
    OR
    (RA >= 20.0 AND
    RA <= 24.0 AND
    _DEC <= -30.0 AND
    _DEC >= -34.0);
'''
GLADE = query_db([db_query])[0]
GLADE_ra = [x[0] for x in GLADE]
GLADE_dec = [x[1] for x in GLADE]
GLADE_z = [x[2] for x in GLADE]
GLADE_B_band = [x[3] for x in GLADE]
GLADE_z_err = [x[4] for x in GLADE]
print("Len GLADE = " + str(len(GLADE_ra)))
print("Time to Complete: " + str(datetime.now() - start))


### Galaxies per GW Pixel/Get Distance and Probs
start = datetime.now()
print("Starting Galaxy per Pixel - " + str(start.time()))
PS1_dist = np.zeros(len(PS1_ra), dtype = float)
PS1_dist_err = np.zeros(len(PS1_ra), dtype = float)
PS1_weight = np.zeros(len(PS1_ra), dtype = float)
GLADE_dist = np.zeros(len(GLADE_ra), dtype = float)
GLADE_dist_err = np.zeros(len(GLADE_ra), dtype = float)
GLADE_weight = np.zeros(len(GLADE_ra), dtype = float)
galaxies_in_pixel = np.zeros(npix, dtype = int)
for i in range(len(PS1_ra)):
    phi = np.deg2rad(PS1_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(PS1_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    PS1_dist[i] = dist_mean[this_pix]
    PS1_dist_err[i] = dist_std[this_pix]
    PS1_weight[i] = 1 - credible_levels[this_pix]
    galaxies_in_pixel[this_pix] = galaxies_in_pixel[this_pix] + 1
for i in range(len(GLADE_ra)):
    phi = np.deg2rad(GLADE_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(GLADE_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    GLADE_dist[i] = dist_mean[this_pix]
    GLADE_dist_err[i] = dist_std[this_pix]
    GLADE_weight[i] = 1 - credible_levels[this_pix]
    galaxies_in_pixel[this_pix] = galaxies_in_pixel[this_pix] + 1
print("Pixel with 1 Galaxy: " + str(len([x for x in galaxies_in_pixel if x == 1])) + ", with 2 Galaxies: " + str(len([x for x in galaxies_in_pixel if x == 2])) + ", with 3+ Galaxies: " + str(len([x for x in galaxies_in_pixel if x >= 3])))
print("Finished Galaxy per Pixel/Get Distance and Probs: " + str(datetime.now() - start))

### H0 Calculation
start = datetime.now()
print("Start H0 Calculations - " + str(start.time()))
H0_dist = [0 for _ in range(len(PS1_ra) + len(GLADE_ra))]
H0_dist_prob = [0 for _ in range(len(H0_dist))]
index_hubble = 0
perc_now = 0
rand_index = random.sample(range(len(H0_dist)), 10)

for i in range(len(PS1_ra)):
    if index_hubble / len(H0_dist) >= perc_now / 100:
        print(str(perc_now) + "%")
        perc_now = perc_now + 10
    phi = np.deg2rad(PS1_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(PS1_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    this_z_dist = np.linspace(PS1_z[i] - 3*PS1_z_err[i], PS1_z[i] + 3*PS1_z_err[i], num)
    this_z_dist_prob = scipy.stats.norm.pdf(x=this_z_dist,loc = PS1_z[i], scale = PS1_z_err[i])
    this_d_dist = np.linspace(PS1_dist[i] - 3*PS1_dist_err[i], PS1_dist[i] + 3*PS1_dist_err[i], num)
    this_d_dist_prob = scipy.stats.norm.pdf(x=this_d_dist,loc = PS1_dist[i], scale = PS1_dist_err[i])
    # Truncate and re-normalize dist
    this_z_dist_prob = [0 if this_z_dist[x] <= 0 else this_z_dist_prob[x] for x in range(len(this_z_dist_prob))]
    this_z_dist_prob = this_z_dist_prob / np.linalg.norm(this_z_dist_prob)
    # integral = np.trapz(this_z_dist_prob, x=this_z_dist)
    # integral = np.trapz([this_z_dist_prob[x] for x in range(len(this_z_dist_prob)) if this_z_dist[x] > 10 ** -4], x=[this_z_dist[x] for x in range(len(this_z_dist)) if this_z_dist[x] > 10 ** -4])
    # if integral != 0:
    #     this_z_dist_prob = this_z_dist_prob/integral
    this_d_dist_prob = [0 if this_d_dist[x] <= 0 else this_d_dist_prob[x] for x in range(len(this_d_dist_prob))]
    this_d_dist_prob = this_d_dist_prob/np.linalg.norm(this_d_dist_prob)
    # integral = np.trapz(this_d_dist_prob, x=this_d_dist)
    # if integral != 0:
    #     this_d_dist_prob = this_d_dist_prob / integral


    # if index_hubble in rand_index or index_hubble == 8608:
    #     plt.figure(index_hubble)
    #     plt.plot(this_z_dist, this_z_dist_prob)
    #     plt.savefig("images/Redshift Sample PDF_"+str(index_hubble)+".png", bbox_inches="tight", dpi=300)
    #     plt.close(index_hubble)


    # Calculate H0 dist
    H0_dist[index_hubble] = c*this_z_dist/this_d_dist
    H0_dist_prob[index_hubble] = (this_z_dist_prob/this_d_dist_prob)*PS1_weight[i]
    index_hubble = index_hubble + 1
for i in range(len(GLADE_ra)):
    if index_hubble / len(H0_dist) >= perc_now / 100:
        print(str(perc_now) + "%")
        perc_now = perc_now + 10
    phi = np.deg2rad(GLADE_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(GLADE_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    this_z_dist = np.linspace(GLADE_z[i] - 3*GLADE_z_err[i], GLADE_z[i] + 3*GLADE_z_err[i], num)
    this_z_dist_prob = scipy.stats.norm.pdf(x=this_z_dist,loc = GLADE_z[i], scale = GLADE_z_err[i])
    this_d_dist = np.linspace(GLADE_dist[i] - 3*GLADE_dist_err[i], GLADE_dist[i] + 3*GLADE_dist_err[i], num)
    this_d_dist_prob = scipy.stats.norm.pdf(x=this_d_dist,loc = GLADE_dist[i], scale = GLADE_dist_err[i])
    # Truncate and re-normalize dist
    this_z_dist_prob = [0 if this_z_dist[x] <= 10**-4 else float(this_z_dist_prob[x]) for x in range(len(this_z_dist_prob))]
    this_z_dist_prob = this_z_dist_prob/np.linalg.norm(this_z_dist_prob)
    # integral = np.trapz([this_z_dist_prob[x] for x in range(len(this_z_dist_prob)) if this_z_dist[x]>10**-4], x=[this_z_dist[x] for x in range(len(this_z_dist)) if this_z_dist[x]>10**-4])
    # integral = np.trapz(this_z_dist_prob, x=this_z_dist)
    # if integral != 0:
    #     this_z_dist_prob = this_z_dist_prob/integral
    this_d_dist_prob = [0 if this_d_dist[x] <= 10**-4 else float(this_d_dist_prob[x]) for x in range(len(this_d_dist_prob))]
    this_d_dist_prob = this_d_dist_prob/np.linalg.norm((this_d_dist_prob))
    # integral = np.trapz(this_d_dist_prob, x=this_d_dist)
    # if integral != 0:
    #     this_d_dist_prob = this_d_dist_prob / integral



    # if index_hubble in rand_index or index_hubble == 8608:
    #     plt.figure(index_hubble)
    #     plt.plot(this_z_dist, this_z_dist_prob)
    #     plt.savefig("images/Redshift Sample PDF_"+str(index_hubble)+".png", bbox_inches="tight", dpi=300)
    #     plt.close(index_hubble)



    # Calculate H0 dist
    H0_dist[index_hubble] = c*this_z_dist/this_d_dist
    H0_dist_prob[index_hubble] = (this_z_dist_prob/this_d_dist_prob)*GLADE_weight[i]
    index_hubble = index_hubble + 1
print("100%")
print("Finished H0 Calculations: " + str(datetime.now() - start))


# Test Some H0
for i in rand_index:
    plt.figure(index_hubble)
    plt.plot(H0_dist[i], H0_dist_prob[i])
    plt.savefig("images/H0 Sample PDF_"+str(i)+".png", bbox_inches="tight", dpi=300)
    plt.close(index_hubble)
# rand_index = random.sample(range(len(H0_dist)), 10)
# plot_i = 0
# for i in rand_index:
#     plt.figure(plot_i)
#     plt.plot(H0_dist[i], H0_dist_prob[i])
#     plt.savefig("images/H0 Sample PDF_"+str(plot_i)+".png", bbox_inches="tight", dpi=300)
#     plt.close(plot_i)
#     plot_i = plot_i + 1





# Add up Hubble Probs
# min_H0 = min([min(H0_dist[x]) for x in range(len(H0_dist))])
# max_H0 = max([max(H0_dist[x]) for x in range(len(H0_dist))])
new_H0_dist = np.linspace(H0_min, H0_max, num)
new_H0_dist_prob = np.zeros(len(new_H0_dist))
for i in range(len(H0_dist)):
    H0_dist_prob[i] = np.interp(new_H0_dist, H0_dist[i], H0_dist_prob[i])
    new_H0_dist_prob = new_H0_dist_prob + H0_dist_prob[i]
new_H0_dist_prob = new_H0_dist_prob/np.trapz(new_H0_dist_prob, x=new_H0_dist)

plt.figure(1)
plt.plot(new_H0_dist,new_H0_dist_prob)
plt.savefig("images/New Dist H0 PDF.png", bbox_inches = "tight", dpi = 300)


print("Finished Script: " + str(datetime.now() - script_start))