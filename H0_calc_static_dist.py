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
num = 10000
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
GLADE_z_err = np.array([x[4] for x in GLADE])
print("Len GLADE = " + str(len(GLADE_ra)))
print("Time to Complete: " + str(datetime.now() - start))

### Redshift Histogram
plt.figure(12)
plt.hist(PS1_z, bins=20, color = "orange",label = "PS1 Photometric Redshift")
plt.hist(GLADE_z, bins=20, color = "aqua", label = "GLADE Spectroscopic Redshift")
plt.title("Redshifts of Each Galaxy")
plt.xlabel("Red Shift")
plt.ylabel("Frequency of Galaxies")
plt.legend()
plt.savefig("images/Redshift Histogram (Super Imposed).png", bbox_inches = "tight", dpi = 300)

plt.figure(13)
counts, bins, bars = plt.hist(PS1_z, bins=20, color = "orange",label = "PS1 Photometric Redshift")
plt.hist(GLADE_z, bins=bins, color = "aqua", label = "GLADE Spectroscopic Redshift", rwidth = 0.8)
plt.title("Redshifts of Each Galaxy")
plt.xlabel("Red Shift")
plt.ylabel("Frequency of Galaxies")
plt.legend()
plt.savefig("images/Redshift Histogram (Combined).png", bbox_inches = "tight", dpi = 300)

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
    PS1_weight[i] = (1 - credible_levels[this_pix])/(PS1_z[i]**3)
    galaxies_in_pixel[this_pix] = galaxies_in_pixel[this_pix] + 1
for i in range(len(GLADE_ra)):
    phi = np.deg2rad(GLADE_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(GLADE_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    GLADE_distance[i] = dist_mean[this_pix]
    GLADE_distance_err[i] = dist_std[this_pix]
    GLADE_weight[i] = (1 - credible_levels[this_pix])/(GLADE_z[i]**3)
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
### Histogram of Distances
plt.figure(14)
plt.hist(np.unique(np.append(PS1_distance, GLADE_distance)), bins=20)
plt.title("Histogram of Unique Distances")
plt.xlabel("Distance (Mpc)")
plt.ylabel("Frequency of Distance")
plt.savefig("images/Distance Histogram.png", bbox_inches = "tight", dpi = 300)


print("Min Redshift: " + str(min(np.append(PS1_z, GLADE_z))))
print("Max Redshift: " + str(max(np.append(PS1_z, GLADE_z))))
print("Min URedshift: " + str(min(np.append(PS1_z - 3*PS1_z_err, GLADE_z - 3*GLADE_z_err))))
print("Max URedshift: " + str(max(np.append(PS1_z + 3*PS1_z_err, GLADE_z + 3*GLADE_z_err))))
print("Min Distance: " + str(min(np.unique(np.append(PS1_distance, GLADE_distance)))))
print("Max Distance: " + str(max(np.unique(np.append(PS1_distance, GLADE_distance)))))
print("Min UDistance: " + str(min(np.unique(np.append(PS1_distance - 3*PS1_distance_err, GLADE_distance- 3*GLADE_distance_err)))))
print("Max UDistance: " + str(max(np.unique(np.append(PS1_distance + 3*PS1_distance_err, GLADE_distance+ 3*GLADE_distance_err)))))


### H0 Calculation
start = datetime.now()
print("Start H0 Calculations - " + str(start.time()))
galaxy_z = np.append(PS1_z, GLADE_z)
galaxy_z_err = np.append(PS1_z_err, GLADE_z_err)
galaxy_distance = np.append(PS1_distance, GLADE_distance)
galaxy_distance_err = np.append(PS1_distance_err, GLADE_distance_err)
galaxy_weight = np.append(PS1_weight, GLADE_weight)

### z-3*z_err = -6, z+3*z_err = 8.5, d-3*d_err = 96, d+3*d_err = 350
z_dist = np.linspace(-6, 8.5, num)
d_dist = np.linspace(96,350, num)
H0_dist = c*z_dist/d_dist
H0_dist_prob = [0 for _ in range(len(galaxy_z))]
perc_now = 0
rand_index = random.sample(range(len(galaxy_z)), 5)

for i in range(len(galaxy_z)):
    this_z_dist_prob = scipy.stats.norm.pdf(x=z_dist,loc = galaxy_z[i], scale = galaxy_z_err[i])
    this_d_dist_prob = scipy.stats.norm.pdf(x=d_dist,loc = galaxy_distance[i], scale = galaxy_distance_err[i])
    # Truncate and re-normalize dist
    this_z_dist_prob = np.array([0.0 if z_dist[x] <= 10**-4 else float(this_z_dist_prob[x]) for x in range(len(this_z_dist_prob))])
    integral = np.trapz(this_z_dist_prob, x=z_dist)
    if integral != 0:
        this_z_dist_prob = this_z_dist_prob/integral
    this_d_dist_prob = np.array([0.0 if d_dist[x] <= 10**-4 else float(this_d_dist_prob[x]) for x in range(len(this_d_dist_prob))])
    integral = np.trapz(this_d_dist_prob, x=d_dist)
    if integral != 0:
        this_d_dist_prob = this_d_dist_prob / integral

    ### SAMPLE REDSHIFT & DISTANCE PLOTS
    if i in rand_index or i == 7248:
        plt.figure(i)
        plt.plot(z_dist, this_z_dist_prob)
        plt.title("Redshift PDF, Galaxy " + str(i))
        plt.xlabel("Redshift")
        plt.ylabel("Probability of Redshift")
        plt.savefig("images/Redshift Sample PDF_"+str(i)+".png", bbox_inches="tight", dpi=300)
        plt.close(i)
    if i in rand_index or i == 7248:
        plt.figure(i)
        plt.plot(d_dist, this_d_dist_prob)
        plt.title("Distance PDF, Galaxy " + str(i))
        plt.xlabel("Distance (Mpc)")
        plt.ylabel("Probability of Distance")
        plt.savefig("images/Distance Sample PDF_"+str(i)+".png", bbox_inches="tight", dpi=300)
        plt.close(i)

    # Calculate H0 dist
    H0_dist_prob[i] = (this_z_dist_prob*this_d_dist_prob)*galaxy_weight[i]

print("100%")
print("Finished H0 Calculations: " + str(datetime.now() - start))


# SAMPLE H0 PLOTS
for i in rand_index+[7248]:
    plt.figure(i)
    plt.plot(H0_dist, H0_dist_prob[i])
    plt.title("Sample H0 PDF, Galaxy " + str(i))
    plt.xlabel("Hubble Constant (km/s/Mpc)")
    plt.ylabel("Probability of H0")
    plt.savefig("images/H0 Sample PDF_"+str(i)+".png", bbox_inches="tight", dpi=300)
    plt.close(i)


# H0 PDF for each Galaxy
print("Begin PDF for each Galaxy")
plt.figure(3)
neg_pdfs = 0
# new_H0_dist = np.array([H0_dist[x] for x in range(len(H0_dist)) if H0_min<=H0_dist[x]<=H0_max])
# for i in range(len(H0_dist_prob)):
#     H0_dist_prob[i] = np.array([0 if H0_dist_prob[i][x] < 10**-4 else float(H0_dist_prob[i][x]) for x in range(len(H0_dist_prob[i])) if H0_min<=H0_dist[x]<=H0_max])
#     integral = np.trapz(H0_dist_prob[i], x=new_H0_dist)
#     if integral != 0:
#         H0_dist_prob[i] = H0_dist_prob[i] / integral
# H0_dist_prob = [x for x in H0_dist_prob if sum(x) >= 10**-4]
for i in range(len(H0_dist_prob)):
    if not any([x for x in H0_dist_prob[i] if x < 0]):
        plt.plot(H0_dist, H0_dist_prob[i])
    else:
        print("Galaxy Index Neg Prob: " +str(i))
        neg_pdfs = neg_pdfs + 1
plt.title("H0 PDF for each Galaxy")
plt.xlabel("Hubble Constant (km/s/Mpc)")
plt.ylabel("Probability of H0")
plt.xlim(0,150)
plt.savefig("images/H0 PDF for each Galaxy.png", bbox_inches="tight", dpi=300)
print("H0 with Neg Probs: " + str(neg_pdfs))





# Add up Hubble Probs
# new_H0_dist = np.linspace(H0_min, H0_max, num)
# new_H0_dist = H0_dist
new_H0_dist_prob = np.zeros(len(H0_dist))
for i in range(len(new_H0_dist_prob)):
    if not any([x for x in H0_dist_prob[i] if x < 0]):
        # new_H0_dist_prob = new_H0_dist_prob + np.interp(new_H0_dist, H0_dist[i], H0_dist_prob[i])
        new_H0_dist_prob = new_H0_dist_prob + new_H0_dist_prob[i]
new_H0_dist = np.array([H0_dist[x] for x in range(len(H0_dist)) if 0<=H0_dist[x]<=140])
new_H0_dist_prob = np.array([new_H0_dist_prob[x] for x in range(len(H0_dist)) if 0<=H0_dist[x]<=140])/np.trapz(np.array([new_H0_dist_prob[x] for x in range(len(H0_dist)) if 0<=H0_dist[x]<=140]), x=new_H0_dist)


print(new_H0_dist)
print(new_H0_dist_prob)
#Get 16%, 50%, 84% confidence interval
interval_16 = -1
interval_50 = -1
interval_84 = -1
running_prob = 0
dx = new_H0_dist[1] - new_H0_dist[0]
for i in range(len(new_H0_dist)):
    running_prob = running_prob + new_H0_dist_prob[i]*dx
    if interval_16 == -1 and running_prob >= 0.16:
        interval_16 = i
    if interval_50 == -1 and running_prob >= 0.50:
        interval_50 = i
    if interval_84 == -1 and running_prob >= 0.84:
        interval_84 = i
frac_measurement = 100*(new_H0_dist[interval_84] - new_H0_dist[interval_16])/(2*new_H0_dist[interval_50])

plt.figure(1)
plt.plot(new_H0_dist,new_H0_dist_prob, color = "blue")
plt.axvline(x = new_H0_dist[interval_16], color = "red", linestyle = "dashed")
plt.axvline(x = new_H0_dist[interval_50], color = "red", label = r"H0 = " + "%0.2f$^{+%0.2f}_{-%0.2f}$ (%0.2f%%)" % (new_H0_dist[interval_50], new_H0_dist[interval_84]-new_H0_dist[interval_50],new_H0_dist[interval_50]-new_H0_dist[interval_16], frac_measurement))
plt.axvline(x = new_H0_dist[interval_84], color = "red", linestyle = "dashed")
plt.title("H0 PDF")
plt.xlabel("Hubble Constant (km/s/Mpc)")
plt.ylabel("Probability of H0")
plt.legend()
plt.savefig("images/New Dist H0 PDF.png", bbox_inches = "tight", dpi = 300)


print("Finished Script: " + str(datetime.now() - script_start))