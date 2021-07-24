import matplotlib.pyplot as plt
import numpy as np
import csv
from datetime import datetime
from database_methods import *
import healpy as hp
from ligo.skymap.postprocess import find_greedy_credible_levels
from mpl_toolkits.basemap import Basemap
from ligo.skymap import distance

PS1_columns = "ObjID", "uniquePspsOBid", "raStack", "decStack", "raMean", "decMean", "ra", "dec", "ng", "gMeanPSFMag", "gMeanPSFMagErr", "gMeanKronMag", "gMeanKronMagErr", "gMeanApMag", "gMeanApMagErr", "nr", "rMeanPSFMag", "rMeanPSFMagErr", "rMeanKronMag", "rMeanKronMagErr", "rMeanApMag", "rMeanApMagErr", "ni", "iMeanPSFMag", "iMeanPSFMagErr", "iMeanKronMag", "iMeanKronMagErr", "iMeanApMag", "iMeanApMagErr", "nz", "zMeanPSFMag", "zMeanPSFMagErr", "zMeanKronMag", "zMeanKronMagErr", "zMeanApMag", "zMeanApMagErr", "ny", "yMeanPSFMag", "yMeanPSFMagErr", "yMeanKronMag", "yMeanKronMagErr", "yMeanApMag", "yMeanApMagErr", "gQfPerfect", "rQfPerfect", "iQfPerfect", "zQfPerfect", "yQfPerfect", "qualityFlag", "objInfoFlag", "primaryDetection", "bestDetection", "class", "prob_Galaxy", "prob_Star", "prob_QSO", "z_phot", "z_photErr", "z_phot0", "extrapolation_Photoz", "ps_score"
GLADE_columns = "id", "Galaxy_id", "Distance_id", "PGC", "Name_GWGC", "Name_HyperLEDA", "Name_2MASS", "Name_SDSS_DR12", "RA", "_Dec", "Coord", "dist", "dist_err", "z_dist", "z_dist_err", "z", "B", "B_err", "B_abs", "J", "J_err", "H", "H_err", "K", "K_err", "flag1", "flag2", "flag3"

### Load in GW Data
prob, distmu, distsigma, distnorm = hp.read_map("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0",field=range(4))
dist_mean, dist_std, norm = distance.parameters_to_moments(distmu, distsigma)
npix = len(prob)
nside = hp.npix2nside(npix)
gw_90_bools = np.zeros(npix, dtype=bool)
credible_levels = find_greedy_credible_levels(prob)

perc_now = 0
start = datetime.now()
print("Loading GW - ", str(start.time()))
for i in range(npix):
    if i / npix >= perc_now / 100:
        print(str(perc_now) + "%")
        perc_now = perc_now + 10
    if credible_levels[i] <= 0.90:
        gw_90_bools[i] = True
print("100%")

alphas = np.array([1 - credible_levels[x] for x in range(len(credible_levels)) if gw_90_bools[x]])
theta, phi = hp.pix2ang(nside, [x for x in range(npix) if gw_90_bools[x]])
gw_ra = [np.rad2deg(x) for x in phi]
gw_dec = [np.rad2deg(0.5 * np.pi - x) for x in theta]
print(type(gw_ra))
print(type(alphas))


print("Len of 90% Confidence Interval: " + str(len(gw_ra)))
print("Time to Complete: " + str(datetime.now() - start))
### SHOW GW RA AND DEC LIMITS
print("First Blob")
print("RA Limits = [" + str(min([x for x in gw_ra if x <= 20])) + ", " + str(max([x for x in gw_ra if x <= 20])) + "]")
print("Dec Limits = [" + str(min([x for x in gw_dec if x >= -30])) + ", " + str(max([x for x in gw_dec if x >= -30])) + "]")
print("Second Blob")
print("RA Limits = [" + str(min([x for x in gw_ra if x >= 20])) + ", " + str(max([x for x in gw_ra if x >= 20])) + "]")
print("Dec Limits = [" + str(min([x for x in gw_dec if x <= -30])) + ", " + str(max([x for x in gw_dec if x <= -30])) + "]")

### Load in PS1 Data
start = datetime.now()
print("Loading PS1 Galaxy " + str(start.time()))
db_query = '''
SELECT * FROM PS1_Galaxy_v4;
'''
PS1 = query_db([db_query])[0]
PS1_ra = np.array([x[4] for x in PS1])
PS1_dec = np.array([x[5] for x in PS1])
PS1_z = np.array([x[56] for x in PS1])
PS1_b_band = np.array([x[11] for x in PS1])
print("Len PS1 = " + str(len(PS1)))
print("Time to Complete: " + str(datetime.now() - start))

### Load in GLADE Data
start = datetime.now()
print("Loading GLADE Galaxy " + str(start.time()))
db_query = '''
SELECT * FROM GalaxyDistance2 
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
GLADE_ra = np.array([x[8] for x in GLADE])
GLADE_dec = np.array([x[9] for x in GLADE])
GLADE_z = np.array([x[15] for x in GLADE])
GLADE_b_band = np.array([x[16] for x in GLADE])
print("Len GLADE = " + str(len(GLADE)))
print("Time to Complete: " + str(datetime.now() - start))


### Cross Match
start = datetime.now()
print("Starting Cross Match - " + str(start.time()))
nums = len(GLADE)
cross_match = [[] for i in range(nums)]
dist_limit = 1 / 3600 # 1/3600 is 1 arcsecond in degrees
print("Distance Limit: " + str(dist_limit * 3600) + " arcseconds")
for i in range(nums):
    right = GLADE_ra[i] + dist_limit
    left = GLADE_ra[i] - dist_limit
    lower = GLADE_dec[i] - dist_limit
    upper = GLADE_dec[i] + dist_limit
    local_PS1_galaxies = [x for x in range(len(PS1)) if (left <= PS1_ra[x] <= right) and (lower <= PS1_dec[x] <= upper)]
    for ps1_index in local_PS1_galaxies:
        dist = np.sqrt(((GLADE_ra[i] - PS1_ra[ps1_index]) ** 2) + ((GLADE_dec[i] - PS1_dec[ps1_index]) ** 2))
        if dist < dist_limit:
            cross_match[i] = cross_match[i] + [ps1_index]
print("Finished Cross Match: " + str(datetime.now() - start))
print("Number of Cross Matches: " + str(len([x for x in cross_match if len(x) > 0])))
print("Number of Cross Matches above 1: " + str(len([x for x in cross_match if len(x) > 1])))

start = datetime.now()
print("Limiting PS1 Data and Making Histogram of Differences - " + str(start.time()))
GLADE_indexs = [x for x in range(len(cross_match)) if len(cross_match[x]) > 0]
PS1_indexes = [cross_match[x][0] for x in range(len(cross_match)) if len(cross_match[x]) > 0]
band_diff = [abs(PS1_b_band[PS1_indexes[x]] - GLADE_b_band[GLADE_indexs[:,x]]) for x in range(len(GLADE_indexs)) if (type(PS1_b_band[PS1_indexes[x]]) == float and type(GLADE_b_band[GLADE_indexs[x]]) == float)]
z_diff = [abs(PS1_z[PS1_indexes[x]] - GLADE_z[GLADE_indexs[x]]) for x in range(len(GLADE_indexs)) if (type(PS1_z[PS1_indexes[x]]) == float and type(GLADE_z[GLADE_indexs[x]]) == float)]


plt.figure(1)
plt.hist([x for x in band_diff if x < 100], bins=20)
plt.title("Histogram of B-band and G-band difference of GLADE and PS1 Cross Match")
plt.xlabel("Band Difference (mags)")
plt.ylabel("Frequency of difference")
plt.savefig("images/Crossmatch Band diff.png", bbox_inches="tight", dpi=300)

plt.figure(5)
plt.hist([x for x in z_diff if x < 0.25], bins=20)
plt.title("Histogram of Red Shift difference of GLADE and PS1 Cross Match")
plt.xlabel("Red Shift Difference")
plt.ylabel("Frequency of difference")
plt.savefig("images/Red Shift diff.png", bbox_inches="tight", dpi=300)

### Get Rid of Cross Match PS1 Galaxies
PS1_good_indexes = [x for x in range(len(PS1)) if x not in PS1_indexes]
PS1 = [PS1[x] for x in PS1_good_indexes]
print("New PS1 Len: " + str(len(PS1)))
print("Finish Limiting PS1 Data: " + str(datetime.now() - start))

print("Its Working")