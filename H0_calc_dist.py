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
counts, bins, bars = plt.hist(np.append(PS1_z,GLADE_z), bins=20, color = "orange",label = "PS1 Photometric Redshift", rwidth=0.9)
plt.hist(GLADE_z, bins=bins, color = "aqua", label = "GLADE Spectroscopic Redshift", rwidth = 0.9)
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
plt.hist(np.append(PS1_distance, GLADE_distance), bins=20, color = "green")
plt.title("Histogram of Distances of Galaxies")
plt.xlabel("Distance (Mpc)")
plt.ylabel("Frequency of Galaxies with Distance")
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
H0_dist = [0 for _ in range(len(PS1_ra) + len(GLADE_ra))]
H0_dist_prob = [0 for _ in range(len(H0_dist))]
index_hubble = 0
perc_now = 0
rand_index = random.sample(range(len(H0_dist)), 5)

for i in range(len(PS1_ra)):
    if index_hubble / len(H0_dist) >= perc_now / 100:
        print(str(perc_now) + "%")
        perc_now = perc_now + 10
    phi = np.deg2rad(PS1_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(PS1_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    this_z_dist = np.linspace(PS1_z[i] - 3*PS1_z_err[i], PS1_z[i] + 3*PS1_z_err[i], num)
    this_z_dist_prob = scipy.stats.norm.pdf(x=this_z_dist,loc = PS1_z[i], scale = PS1_z_err[i])
    this_d_dist = np.linspace(PS1_distance[i] - 3*PS1_distance_err[i], PS1_distance[i] + 3*PS1_distance_err[i], num)
    this_d_dist_prob = scipy.stats.norm.pdf(x=this_d_dist,loc = PS1_distance[i], scale = PS1_distance_err[i])
    # Truncate and re-normalize dist
    this_z_dist_prob = [0 if this_z_dist[x] <= 10**-4 else float(this_z_dist_prob[x]) for x in range(len(this_z_dist_prob))]
    integral = np.trapz(this_z_dist_prob, x=this_z_dist)
    if integral != 0:
        this_z_dist_prob = this_z_dist_prob/integral
    this_d_dist_prob = [0 if this_d_dist[x] <= 10**-4 else float(this_d_dist_prob[x]) for x in range(len(this_d_dist_prob))]
    integral = np.trapz(this_d_dist_prob, x=this_d_dist)
    if integral != 0:
        this_d_dist_prob = this_d_dist_prob / integral



    # ### SAMPLE REDSHIFT & DISTANCE PLOTS
    # if index_hubble in rand_index or index_hubble == 7248:
    #     plt.figure(index_hubble)
    #     plt.plot(this_z_dist, this_z_dist_prob)
    #     plt.title("Redshift PDF, Galaxy " + str(index_hubble))
    #     plt.xlabel("Redshift")
    #     plt.ylabel("Probability of Redshift")
    #     plt.savefig("images/Redshift Sample PDF_"+str(index_hubble)+".png", bbox_inches="tight", dpi=300)
    #     plt.close(index_hubble)
    # if index_hubble in rand_index or index_hubble == 7248:
    #     plt.figure(index_hubble)
    #     plt.plot(this_d_dist, this_d_dist_prob)
    #     plt.title("Distance PDF, Galaxy " + str(index_hubble))
    #     plt.xlabel("Distance (Mpc)")
    #     plt.ylabel("Probability of Distance")
    #     plt.savefig("images/Distance Sample PDF_"+str(index_hubble)+".png", bbox_inches="tight", dpi=300)
    #     plt.close(index_hubble)




    # Calculate H0 dist
    H0_dist[index_hubble] = c * this_z_dist / this_d_dist
    H0_dist_prob[index_hubble] = (this_z_dist_prob * this_d_dist_prob) * PS1_weight[i]
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
    this_d_dist = np.linspace(GLADE_distance[i] - 3*GLADE_distance_err[i], GLADE_distance[i] + 3*GLADE_distance_err[i], num)
    this_d_dist_prob = scipy.stats.norm.pdf(x=this_d_dist,loc = GLADE_distance[i], scale = GLADE_distance_err[i])
    # Truncate and re-normalize dist
    this_z_dist_prob = [0 if this_z_dist[x] <= 10**-4 else float(this_z_dist_prob[x]) for x in range(len(this_z_dist_prob))]
    integral = np.trapz(this_z_dist_prob, x=this_z_dist)
    if integral != 0:
        this_z_dist_prob = this_z_dist_prob/integral
    this_d_dist_prob = [0 if this_d_dist[x] <= 10**-4 else float(this_d_dist_prob[x]) for x in range(len(this_d_dist_prob))]
    integral = np.trapz(this_d_dist_prob, x=this_d_dist)
    if integral != 0:
        this_d_dist_prob = this_d_dist_prob / integral



    # ### SAMPLE REDSHIFT & DISTANCE PLOTS
    # if index_hubble in rand_index or index_hubble == 7248:
    #     plt.figure(index_hubble)
    #     plt.plot(this_z_dist, this_z_dist_prob)
    #     plt.title("Redshift PDF, Galaxy " + str(index_hubble))
    #     plt.xlabel("Redshift")
    #     plt.ylabel("Probability of Redshift")
    #     plt.savefig("images/Redshift Sample PDF_"+str(index_hubble)+".png", bbox_inches="tight", dpi=300)
    #     plt.close(index_hubble)
    # if index_hubble in rand_index or index_hubble == 7248:
    #     plt.figure(index_hubble)
    #     plt.plot(this_d_dist, this_d_dist_prob)
    #     plt.title("Distance PDF, Galaxy " + str(index_hubble))
    #     plt.xlabel("Distance (Mpc)")
    #     plt.ylabel("Probability of Distance")
    #     plt.savefig("images/Distance Sample PDF_"+str(index_hubble)+".png", bbox_inches="tight", dpi=300)
    #     plt.close(index_hubble)



    # Calculate H0 dist
    H0_dist[index_hubble] = c * this_z_dist / this_d_dist
    H0_dist_prob[index_hubble] = (this_z_dist_prob * this_d_dist_prob) * GLADE_weight[i]
    index_hubble = index_hubble + 1
print("100%")
print("Finished H0 Calculations: " + str(datetime.now() - start))


# SAMPLE H0 PLOTS
# for i in rand_index+[7248]:
#     plt.figure(index_hubble)
#     plt.plot(H0_dist[i], H0_dist_prob[i])
#     plt.title("Sample H0 PDF, Galaxy " + str(i))
#     plt.xlabel("Hubble Constant (km/s/Mpc)")
#     plt.ylabel("Probability of H0")
#     plt.savefig("images/H0 Sample PDF_"+str(i)+".png", bbox_inches="tight", dpi=300)
#     plt.close(index_hubble)


# H0 PDF for each Galaxy
plt.figure(3)
neg_pdfs = 0
for i in range(len(H0_dist)):
    if not any([x for x in H0_dist_prob[i] if x < 0]):
        plt.plot(H0_dist[i], H0_dist_prob[i])
    else:
        neg_pdfs = neg_pdfs + 1
plt.title("H0 PDF for each Galaxy")
plt.xlabel("Hubble Constant (km/s/Mpc)")
plt.ylabel("Probability of H0")
plt.xlim(0,150)
plt.savefig("images/H0 PDF for each Galaxy.png", bbox_inches="tight", dpi=300)
print("H0 with Neg Probs: " + str(neg_pdfs))

# Add up Hubble Probs
start = datetime.now()
print("Starting Adding Hubble Probs - " + str(start.time()))
new_H0_dist_bin = np.linspace(H0_min, H0_max, num+1)
new_H0_dist = np.array([new_H0_dist_bin[i]+((new_H0_dist_bin[1]-new_H0_dist_bin[0])/2) for i in range(len(new_H0_dist_bin)-1)])
new_H0_dist_prob = np.zeros(len(new_H0_dist))
perc_now = 0
for i in range(len(H0_dist)):
    if i / len(H0_dist) >= perc_now / 100:
        print(str(perc_now) + "%")
        perc_now = perc_now + 10
    if not any([x for x in H0_dist_prob[i] if x < 0]):
        for ii in [x for x in range(len(H0_dist[i])) if H0_min<=H0_dist[i][x]<=H0_max]:
            for new_index in range(len(new_H0_dist)):
                if new_H0_dist_bin[new_index]<=H0_dist[i][ii]<new_H0_dist_bin[new_index+1]:
                    new_H0_dist_prob[new_index] = new_H0_dist_prob[new_index] + H0_dist_prob[i][ii]
            # new_H0_dist_prob[ii] = new_H0_dist_prob + np.interp(new_H0_dist, H0_dist[i], H0_dist_prob[i])
new_H0_dist_prob = new_H0_dist_prob/np.trapz(new_H0_dist_prob, x=new_H0_dist)
print("Finished Adding Hubble Probs: " + str(datetime.now()-start))

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
plt.axvline(x = new_H0_dist[interval_50], color = "red", label = r"H_{0} = " + "%0.2f$^{+%0.2f}_{-%0.2f}$ km/s/Mpc (%0.2f%%)" % (new_H0_dist[interval_50], new_H0_dist[interval_84]-new_H0_dist[interval_50],new_H0_dist[interval_50]-new_H0_dist[interval_16], frac_measurement))
plt.axvline(x = new_H0_dist[interval_84], color = "red", linestyle = "dashed")
plt.title("H0 PDF")
plt.xlabel("Hubble Constant (km/s/Mpc)")
plt.ylabel("Probability of H0")
plt.legend()
plt.savefig("images/New Dist H0 PDF.png", bbox_inches = "tight", dpi = 300)


print("Finished Script: " + str(datetime.now() - script_start))