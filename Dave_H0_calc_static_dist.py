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
H0_min = 20
H0_max = 150
num = 500
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
PS1_ra = np.array([x[0] for x in PS1])
PS1_dec = np.array([x[1] for x in PS1])
PS1_z = np.array([x[2] for x in PS1])
PS1_class = np.array([x[3] for x in PS1])
PS1_objid = np.array([x[4] for x in PS1])
PS1_ps_score = np.array([x[5] for x in PS1])
PS1_z_err = np.array([x[6] for x in PS1])
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
GLADE_ra = np.array([x[0] for x in GLADE])
GLADE_dec = np.array([x[1] for x in GLADE])
GLADE_z = np.array([x[2] for x in GLADE])
GLADE_B_band = np.array([x[3] for x in GLADE])
# GLADE_z_err = np.array([x[4] for x in GLADE])
GLADE_z_err = np.array([0.001 for _ in range(len(GLADE_ra))])
print("Len GLADE = " + str(len(GLADE_ra)))
print("Time to Complete: " + str(datetime.now() - start))

print("Avg GLADE z" + str(np.mean(GLADE_z_err)))
print("Avg PS1 z" + str(np.mean(PS1_z_err)))
# raise Exception


### Galaxies per GW Pixel/Get Distance and Probs
start = datetime.now()
print("Starting Galaxy per Pixel - " + str(start.time()))
PS1_distance = np.zeros(len(PS1_ra), dtype = float)
PS1_distance_err = np.zeros(len(PS1_ra), dtype = float)
PS1_weight = np.zeros(len(PS1_ra), dtype = float)
GLADE_distance = np.zeros(len(GLADE_ra), dtype = float)
GLADE_distance_err = np.zeros(len(GLADE_ra), dtype = float)
GLADE_weight = np.zeros(len(GLADE_ra), dtype = float)
galaxies_in_pixel = np.zeros(npix, dtype = int)
for i in range(len(PS1_ra)):
    phi = np.deg2rad(PS1_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(PS1_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    PS1_distance[i] = dist_mean[this_pix]
    PS1_distance_err[i] = dist_std[this_pix]
    PS1_weight[i] = (1 - credible_levels[this_pix])
    galaxies_in_pixel[this_pix] = galaxies_in_pixel[this_pix] + 1
for i in range(len(GLADE_ra)):
    phi = np.deg2rad(GLADE_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(GLADE_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    GLADE_distance[i] = dist_mean[this_pix]
    GLADE_distance_err[i] = dist_std[this_pix]
    GLADE_weight[i] = (1 - credible_levels[this_pix])
    galaxies_in_pixel[this_pix] = galaxies_in_pixel[this_pix] + 1
print("Pixel with 1 Galaxy: " + str(len([x for x in galaxies_in_pixel if x == 1])) + ", with 2 Galaxies: " + str(len([x for x in galaxies_in_pixel if x == 2])) + ", with 3+ Galaxies: " + str(len([x for x in galaxies_in_pixel if x >= 3])))
# Divide Weights by Galaxies in Pixel
for i in range(len(PS1_ra)):
    phi = np.deg2rad(PS1_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(PS1_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    PS1_weight[i] = PS1_weight[i]/galaxies_in_pixel[this_pix]
for i in range(len(GLADE_ra)):
    phi = np.deg2rad(GLADE_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(GLADE_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    GLADE_weight[i] = GLADE_weight[i]/galaxies_in_pixel[this_pix]
print("Finished Galaxy per Pixel/Get Distance and Probs: " + str(datetime.now() - start))


print("Min Redshift: " + str(min(np.append(PS1_z, GLADE_z))))
print("Max Redshift: " + str(max(np.append(PS1_z, GLADE_z))))
print("Min URedshift: " + str(min(np.append(PS1_z - 3*PS1_z_err, GLADE_z - 3*GLADE_z_err))))
print("Max URedshift: " + str(max(np.append(PS1_z + 3*PS1_z_err, GLADE_z + 3*GLADE_z_err))))
print("Min Distance: " + str(min(np.unique(np.append(PS1_distance, GLADE_distance)))))
print("Max Distance: " + str(max(np.unique(np.append(PS1_distance, GLADE_distance)))))
print("Min UDistance: " + str(min(np.unique(np.append(PS1_distance - 3*PS1_distance_err, GLADE_distance- 3*GLADE_distance_err)))))
print("Max UDistance: " + str(max(np.unique(np.append(PS1_distance + 3*PS1_distance_err, GLADE_distance+ 3*GLADE_distance_err)))))


redshift_arr = np.linspace(0.005, 0.225, num)
distance_arr = np.linspace(100, 300, num)
H0_arr = (c * redshift_arr) / distance_arr

H0_distributions = []
H0_weights = []

galaxy_z = np.append(PS1_z, GLADE_z)
galaxy_z_err = np.append(PS1_z_err, GLADE_z_err)
galaxy_distance = np.append(PS1_distance, GLADE_distance)
galaxy_distance_err = np.append(PS1_distance_err, GLADE_distance_err)
galaxy_weight = np.append(PS1_weight, GLADE_weight)

print("Redshift Sample:", galaxy_z[408], galaxy_z_err[408])


perc_now = 0
rand_index = random.sample(range(len(galaxy_z)), 4)
start = datetime.now()
print("Start Calculating H0 Distributions - " + str(start.time()))
for i in range(len(galaxy_z)):
    if i/len(galaxy_z) >= perc_now/100:
        print(str(perc_now) + "%")
        perc_now = perc_now + 10

    if i == 408:
        z_dist = scipy.stats.norm.pdf(x=np.linspace(0, galaxy_z[i] + 5*galaxy_z_err[i], num), loc=galaxy_z[i], scale=galaxy_z_err[i])
        d_dist = scipy.stats.norm.pdf(x=np.linspace(galaxy_distance[i] - 5*galaxy_distance_err[i], galaxy_distance[i] + 5*galaxy_distance_err[i], num), loc=galaxy_distance[i], scale=galaxy_distance_err[i])
        H0_arr = (c * np.linspace(0, galaxy_z[i] + 5*galaxy_z_err[i], num)) / np.linspace(galaxy_distance[i]- 5*galaxy_distance_err[i], galaxy_distance[i] + 5*galaxy_distance_err[i], num)
        redshift_arr = np.linspace(0, galaxy_z[i] + 5*galaxy_z_err[i], num)
    else:
        z_dist = scipy.stats.norm.pdf(x = redshift_arr, loc = galaxy_z[i], scale = galaxy_z_err[i])
        d_dist = scipy.stats.norm.pdf(x = distance_arr, loc = galaxy_distance[i], scale = galaxy_distance_err[i])
        H0_arr = (c * redshift_arr) / distance_arr
        redshift_arr = np.linspace(0.005, 0.225, num)

    z_norm = z_dist#/np.trapz(z_dist, x =redshift_arr)
    d_norm = d_dist#/np.trapz(d_dist, x=distance_arr)

    H0_distribution = c*((z_norm/(redshift_arr**3)) / d_norm)
    H0_norm = H0_distribution / np.trapz(H0_distribution, x=H0_arr)
    H0_distributions.append(H0_norm * galaxy_weight[i])


    ### SAMPLE REDSHIFT & DISTANCE PLOTS
    if i in rand_index or i == 408:
        plt.figure(i)
        plt.plot(redshift_arr, z_dist)
        plt.title("Redshift PDF, Galaxy " + str(i))
        plt.xlabel("Redshift")
        plt.ylabel("Probability of Redshift")
        plt.savefig("images/Redshift Sample PDF_" + str(i) + ".png", bbox_inches="tight", dpi=300)
        plt.close(i)
        plt.figure(i)
        plt.plot(distance_arr, d_dist)
        plt.title("Distance PDF, Galaxy " + str(i))
        plt.xlabel("Distance (Mpc)")
        plt.ylabel("Probability of Distance")
        plt.savefig("images/Distance Sample PDF_" + str(i) + ".png", bbox_inches="tight", dpi=300)
        plt.close(i)
# SAMPLE H0 PLOTS
for i in rand_index+[408]:
    plt.figure(i)
    plt.plot(H0_arr, H0_distributions[i])
    plt.title("Sample H0 PDF, Galaxy " + str(i))
    plt.xlabel("Hubble Constant (km/s/Mpc)")
    plt.ylabel("Probability of H0")
    # if i == 408:
    #     plt.xlim(20,150)
    plt.savefig("images/H0 Sample PDF_"+str(i)+".png", bbox_inches="tight", dpi=300)
    plt.close(i)
print("100%")
print("Finished H0 Calculations: " + str(datetime.now() - start))

### H0 PDF of Each Galaxy Superimposed
print("H0 PDF of Each Galaxy Superimposed")
plt.figure(76, figsize = (10,4))
for i in range(len(H0_distributions[:1334])):
    plt.plot(H0_arr, H0_distributions[i])
    if any([x for x in H0_distributions[i] if x > 5]):
        print("Huge H0 PDF Galaxy: " + str(i))
plt.title("All Galaxy H0 PDF - PS1")
plt.xlabel("Hubble Constant (km/s/Mpc)")
plt.ylabel("Probability of H0")
# plt.ylim(0,5)
plt.savefig("images/Dave H0 Superimposed - PS1.png", bbox_inches="tight", dpi=300)
plt.close(76)

plt.figure(76, figsize = (10,4))
for i in range(len(H0_distributions[1334:])):
    plt.plot(H0_arr, H0_distributions[i])
    if any([x for x in H0_distributions[i] if x > 5]):
        print("Huge H0 PDF Galaxy: " + str(i))
plt.title("All Galaxy H0 PDF - GLADE")
plt.xlabel("Hubble Constant (km/s/Mpc)")
plt.ylabel("Probability of H0")
# plt.ylim(0,5)
plt.savefig("images/Dave H0 Superimposed - GLADE.png", bbox_inches="tight", dpi=300)
plt.close(76)


final_h0_dist = np.zeros(len(H0_arr))

for h0_tuple in H0_distributions:
   final_h0_dist += h0_tuple

norm_final_h0 = final_h0_dist/np.trapz(final_h0_dist, x=H0_arr)

plt.figure(99)
plt.plot(H0_arr, norm_final_h0)
plt.title("Final H0 PDF")
plt.xlabel("Hubble Constant (km/s/Mpc)")
plt.ylabel("Probability of H0")
plt.savefig("images/Dave Final H0 PDF - Before Trunk.png", bbox_inches="tight", dpi=300)

norm_final_h0 = np.array([norm_final_h0[x] for x in range(len(H0_arr)) if 20<=H0_arr[x]<=150])
H0_arr = np.array([H0_arr[x] for x in range(len(H0_arr)) if 20<=H0_arr[x]<=150])
norm_final_h0 = norm_final_h0/np.trapz(norm_final_h0,x=H0_arr)

plt.figure(109)
plt.plot(H0_arr, norm_final_h0)
plt.title("Final H0 PDF")
plt.xlabel("Hubble Constant (km/s/Mpc)")
plt.ylabel("Probability of H0")
plt.savefig("images/Dave Final H0 PDF.png", bbox_inches="tight", dpi=300)