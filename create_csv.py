import astropy.cosmology
import matplotlib.pyplot as plt
import numpy as np
import csv
from datetime import datetime
from database_methods import *
import healpy as hp
from ligo.skymap.postprocess import find_greedy_credible_levels
from mpl_toolkits.basemap import Basemap
from ligo.skymap import distance
from astropy.cosmology import *
import astropy.units as u

script_start = datetime.now()

plotting = True
write_csv = True

PS1_columns = ["ObjID", "uniquePspsOBid", "raStack", "decStack", "raMean", "decMean", "ra", "dec", "ng", "gMeanPSFMag", "gMeanPSFMagErr", "gMeanKronMag", "gMeanKronMagErr", "gMeanApMag", "gMeanApMagErr", "nr", "rMeanPSFMag", "rMeanPSFMagErr", "rMeanKronMag", "rMeanKronMagErr", "rMeanApMag", "rMeanApMagErr", "ni", "iMeanPSFMag", "iMeanPSFMagErr", "iMeanKronMag", "iMeanKronMagErr", "iMeanApMag", "iMeanApMagErr", "nz", "zMeanPSFMag", "zMeanPSFMagErr", "zMeanKronMag", "zMeanKronMagErr", "zMeanApMag", "zMeanApMagErr", "ny", "yMeanPSFMag", "yMeanPSFMagErr", "yMeanKronMag", "yMeanKronMagErr", "yMeanApMag", "yMeanApMagErr", "gQfPerfect", "rQfPerfect", "iQfPerfect", "zQfPerfect", "yQfPerfect", "qualityFlag", "objInfoFlag", "primaryDetection", "bestDetection", "class", "prob_Galaxy", "prob_Star", "prob_QSO", "z_phot", "z_photErr", "z_phot0", "extrapolation_Photoz", "ps_score"]
PS1_new_columns = ["ObjID", "raMean", "decMean", "class", "prob_Galaxy", "prob_Star", "prob_QSO", "z_phot", "z_photErr", "ps_score"]
GLADE_columns = ["id", "Galaxy_id", "Distance_id", "PGC", "Name_GWGC", "Name_HyperLEDA", "Name_2MASS", "Name_SDSS_DR12", "RA", "_Dec", "Coord", "dist", "dist_err", "z_dist", "z_dist_err", "z", "B", "B_err", "B_abs", "J", "J_err", "H", "H_err", "K", "K_err", "flag1", "flag2", "flag3"]
GLADE_good_columns = [0, 1, 2, 8, 9, 13, 14, 15, 16]

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

print("Len of 90% Confidence Interval: " + str(len(gw_ra)))
print("Time to Complete: " + str(datetime.now() - start))
# Show GW RA and DEC Limits
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
SELECT * FROM PS1_Galaxy_PS;
'''
PS1 = query_db([db_query])[0]
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
print("Len GLADE = " + str(len(GLADE)))
print("Time to Complete: " + str(datetime.now() - start))


### Only Galaxies in 90% Intervals
PS1_ra = np.array([x[1] for x in PS1])
PS1_dec = np.array([x[2] for x in PS1])

GLADE_ra = np.array([x[8] for x in GLADE])
GLADE_dec = np.array([x[9] for x in GLADE])

start = datetime.now()
print("Start Limit to GW Zone - " + str(start.time()))
PS1_bools = np.ones(len(PS1), dtype=bool)
GLADE_bools = np.ones(len(GLADE), dtype=bool)

for i in range(len(PS1_bools)):
    phi = np.deg2rad(PS1_ra[i])
    theta = 0.5*np.pi - np.deg2rad(PS1_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    if credible_levels[this_pix] > 0.90:
        PS1_bools[i] = False
for i in range(len(GLADE_bools)):
    phi = np.deg2rad(GLADE_ra[i])
    theta = 0.5*np.pi - np.deg2rad(GLADE_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    if credible_levels[this_pix] > 0.90:
        GLADE_bools[i] = False

PS1 = [PS1[x] for x in range(len(PS1)) if PS1_bools[x]]
print("New Len PS1: " + str(len(PS1)))

GLADE = [GLADE[x] for x in range(len(GLADE)) if GLADE_bools[x]]
print("New Len GLADE: " + str(len(GLADE)))
print("Finished Limit to GW Zone: " + str(datetime.now() - start))




### Redshift Limit
PS1_ra = np.array([x[1] for x in PS1])
PS1_dec = np.array([x[2] for x in PS1])
PS1_z = np.array([x[7] for x in PS1])
PS1_z_err = np.array([x[8] for x in PS1])
GLADE_ra = np.array([x[8] for x in GLADE])
GLADE_dec = np.array([x[9] for x in GLADE])
GLADE_z = np.array([x[15] for x in GLADE])
GLADE_z_err = np.full(len(GLADE_z), 10**-4)

start = datetime.now()
print("Start Limit Redshift - " + str(start.time()))
PS1_bools = np.ones(len(PS1), dtype=bool)
GLADE_bools = np.ones(len(GLADE), dtype=bool)
PS1e_bools = np.ones(len(PS1), dtype=bool)
c = 299792.458
cosmo_high = astropy.cosmology.LambdaCDM(H0=20.0, Om0=0.27, Ode0=0.73)
cosmo_low = astropy.cosmology.LambdaCDM(H0=150.0, Om0=0.27, Ode0=0.73)
rel_uncert_max = 0.5
abs_z_min = 0.01
abs_z_max = 0.15

wind_i = 0
rel_i = 0
abs_i = 0

for i in range(len(PS1_bools)):
    phi = np.deg2rad(PS1_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(PS1_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    if abs(PS1_z_err[i]/PS1_z[i]) > rel_uncert_max:
        PS1_bools[i] = False
        PS1e_bools[i] = False
        rel_i = rel_i + 1
    if PS1_bools[i] or PS1e_bools[i]:
        if PS1_z[i] > abs_z_max or PS1_z[i] < abs_z_min:
            PS1_bools[i] = False
            abs_i = abs_i + 1
        if PS1_z[i] > 0.3 or PS1_z[i] < abs_z_min:
            PS1e_bools[i] = False
    if PS1_bools[i]:
        if PS1_z[i]+PS1_z_err[i] < z_at_value(cosmo_high.luminosity_distance, (dist_mean[this_pix] - 2 * dist_std[this_pix])*u.Mpc) or PS1_z[i]-PS1_z_err[i] > z_at_value(cosmo_low.luminosity_distance, (dist_mean[this_pix] + 2 * dist_std[this_pix])*u.Mpc):
            PS1_bools[i] = False
            # PS1e_bools[i] = False
            wind_i = wind_i + 1
for i in range(len(GLADE_bools)):
    phi = np.deg2rad(GLADE_ra[i])
    theta = 0.5 * np.pi - np.deg2rad(GLADE_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    if abs(GLADE_z_err[i]/GLADE_z[i]) > rel_uncert_max:
        GLADE_bools[i] = False
        rel_i = rel_i + 1
    if GLADE_bools[i]:
        if GLADE_z[i] > abs_z_max or GLADE_z[i] < abs_z_min:
            PS1_bools[i] = False
            abs_i = abs_i + 1
    if GLADE_bools[i]:
        if GLADE_z[i]+GLADE_z_err[i] < z_at_value(cosmo_high.luminosity_distance, (dist_mean[this_pix] - 2 * dist_std[this_pix])*u.Mpc) or GLADE_z[i]-GLADE_z_err[i] > z_at_value(cosmo_low.luminosity_distance, (dist_mean[this_pix] + 2 * dist_std[this_pix])*u.Mpc):
            GLADE_bools[i] = False
            wind_i = wind_i + 1



# for i in range(len(PS1_bools)):
#     phi = np.deg2rad(PS1_ra[i])
#     theta = 0.5 * np.pi - np.deg2rad(PS1_dec[i])
#     this_pix = hp.ang2pix(nside, theta, phi)
#     if (PS1_z[i] + PS1_z_err[i] < (20 * (dist_mean[this_pix] - 2 * dist_std[this_pix])) / c) or (PS1_z[i] - PS1_z_err[i] > (150 * (dist_mean[this_pix] + 2 * dist_std[this_pix])) / c):
#         PS1_bools[i] = False
#     if PS1_bools[i]:
#         if abs(PS1_z_err[i]/PS1_z[i]) > rel_uncert_max:
#             PS1_bools[i] = False
#             rel_i = rel_i + 1
#     if PS1_bools[i]:
#         if PS1_z[i] > abs_z_max or PS1_z[i] < abs_z_min:
#             PS1_bools[i] = False
#             abs_i = abs_i + 1
#
# for i in range(len(GLADE_bools)):
#     phi = np.deg2rad(GLADE_ra[i])
#     theta = 0.5*np.pi - np.deg2rad(GLADE_dec[i])
#     this_pix = hp.ang2pix(nside, theta, phi)
#     if (GLADE_z[i]+GLADE_z_err[i] < (20 * (dist_mean[this_pix] - 2*dist_std[this_pix]))/c) or (GLADE_z[i]-GLADE_z_err[i] > (150 * (dist_mean[this_pix] + 2*dist_std[this_pix]))/c):
#         GLADE_bools[i] = False
#     if GLADE_bools[i]:
#         if abs(GLADE_z_err[i]/GLADE_z[i]) > rel_uncert_max:
#             GLADE_bools[i] = False
#             rel_i = rel_i + 1
#     if GLADE_bools[i]:
#         if GLADE_z[i] > abs_z_max or GLADE_z[i] < abs_z_min:
#             PS1_bools[i] = False
#             abs_i = abs_i + 1

print("PS1 Redshift out of bounds: " + str((len([x for x in PS1_bools if not x])/len(PS1_bools))*100) + "%")
print("PS1e Redshift out of bounds: " + str((len([x for x in PS1e_bools if not x])/len(PS1e_bools))*100) + "%")
print("GLADE Redshift out of bounds: " + str((len([x for x in GLADE_bools if not x])/len(GLADE_bools))*100) + "%")
print("Window Cut Off: " + str(wind_i))
print("Relative Cut Off: " + str(rel_i))
print("Abs Cut Off: " + str(abs_i))

print("Redshift Limit PS1e Bools: " + str(len([x for x in PS1e_bools if x])))
PS1e = [PS1[x] for x in range(len(PS1)) if PS1e_bools[x]]
PS1 = [PS1[x] for x in range(len(PS1)) if PS1_bools[x]]
print("New Len PS1: " + str(len(PS1)))
print("New Len PS1e: " + str(len(PS1e)))
print("PS1 Min z = " + str(min([x[7] for x in PS1])) + ", Max z = " + str(max([x[7] for x in PS1])))

GLADE = [GLADE[x] for x in range(len(GLADE)) if GLADE_bools[x]]
print("New Len GLADE: " + str(len(GLADE)))
print("GLADE Min z = " + str(min([x[15] for x in GLADE])) + ", Max z = " + str(max([x[15] for x in GLADE])))
print("Finished Limit Redshift: " + str(datetime.now() - start))

# print("After Redshift Limit PS1e",PS1e)

### Cross Match
PS1_ra = np.array([x[1] for x in PS1])
PS1_dec = np.array([x[2] for x in PS1])
PS1e_ra = np.array([x[1] for x in PS1e])
PS1e_dec = np.array([x[2] for x in PS1e])
GLADE_ra = np.array([x[8] for x in GLADE])
GLADE_dec = np.array([x[9] for x in GLADE])


start = datetime.now()
print("Starting Cross Match - " + str(start.time()))
nums = len(GLADE)
cross_match = [[] for i in range(nums)]
cross_match_e = [[] for i in range(nums)]
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
    local_PS1e_galaxies = [x for x in range(len(PS1e)) if (left <= PS1e_ra[x] <= right) and (lower <= PS1e_dec[x] <= upper)]
    for ps1e_index in local_PS1e_galaxies:
        dist = np.sqrt(((GLADE_ra[i] - PS1e_ra[ps1_index]) ** 2) + ((GLADE_dec[i] - PS1e_dec[ps1_index]) ** 2))
        if dist < dist_limit:
            cross_match_e[i] = cross_match_e[i] + [ps1e_index]
print("Finished Cross Match: " + str(datetime.now() - start))
print("Number of Cross Matches: " + str(len([x for x in cross_match if len(x) > 0])))
print("Number of Cross Matches_e: " + str(len([x for x in cross_match_e if len(x) > 0])))
print("Number of Cross Matches above 1: " + str(len([x for x in cross_match if len(x) > 1])))
print("Number of Cross Matches_e above 1: " + str(len([x for x in cross_match_e if len(x) > 1])))


start = datetime.now()
print("Limiting PS1 Data and Making Histogram of Differences - " + str(start.time()))


# Get Rid of Cross Match PS1 Galaxies
PS1_indexes = [cross_match[x][0] for x in range(len(cross_match)) if len(cross_match[x]) > 0]
PS1e_indexes = [cross_match_e[x][0] for x in range(len(cross_match_e)) if len(cross_match_e[x]) > 0]
PS1_good_indexes = [x for x in range(len(PS1)) if x not in PS1_indexes]
PS1e_good_indexes = [x for x in range(len(PS1e)) if x not in PS1e_indexes]
PS1 = [PS1[x] for x in PS1_good_indexes]
PS1e = [PS1e[x] for x in PS1e_good_indexes]
print("New PS1 Len: " + str(len(PS1)))
print("Finish Limiting PS1 Data: " + str(datetime.now() - start))







# print("After Cross Match PS1e",PS1e)

PS1_z = np.array([x[7] for x in PS1])
PS1e_z = np.array([x[7] for x in PS1e])
GLADE_z = np.array([x[15] for x in GLADE])
print("Max PS1e Redshift", max(PS1e_z))
print("Min Redshift: " + str(min(np.append(GLADE_z, PS1_z))))
print("Min Redshift All: " + str(min(np.append(np.append(GLADE_z, PS1_z), PS1e_z))))

### PLOTTING REDSHIFT HISTOGRAM
print("Plotting Final Redshift Histogram")
fig, ax = plt.subplots()
ax2 = ax.twinx()
counts, bins, bars = ax.hist(PS1e_z[np.where(PS1e_z>=0.15)], bins=np.linspace(0,0.3,21), color = "red",label = "PS1 Photometric Redshift (Extended)", rwidth=0.9)
ax.hist(np.append(GLADE_z,PS1_z), bins=bins, color = "aqua",label = "GLADE Spectroscopic Redshift", rwidth=0.9)
ax.hist(PS1_z, bins=bins, color = "orange", label = "PS1 Photometric Redshift", rwidth = 0.9)
ax2.plot(bins, bins**3, color = "green", label = "Redshift Cubed")
ax2.set_ylabel("Redshift Cubed")
plt.title("Redshifts of Each Galaxy")
ax.set_xlabel("Redshift")
ax.set_ylabel("Frequency of Galaxies")
ax.legend(loc="upper left")
plt.savefig("images/Redshift Histogram CSV.png", bbox_inches = "tight", dpi = 300)
print("Finished Plotting Final Redshift Histogram")


### PLOTTING
if plotting:
    PS1_ra = np.array([x[1] for x in PS1])
    PS1_dec = np.array([x[2] for x in PS1])
    PS1_z = np.array([x[7] for x in PS1])
    PS1_ps_score = np.array([x[9] for x in PS1])
    GLADE_ra = np.array([x[8] for x in GLADE])
    GLADE_dec = np.array([x[9] for x in GLADE])
    GLADE_z = np.array([x[15] for x in GLADE])

    start = datetime.now()
    print("Start Plotting " + str(start.time()))
    plt.figure(2)
    map = Basemap(width=3*(10**6),height=3*(10**6)*0.75,projection='lcc', resolution='c',lat_0=np.mean(gw_dec),lon_0=np.mean(gw_ra))

    p_x, p_y = map(PS1_ra,PS1_dec)
    g_x, g_y = map(GLADE_ra,GLADE_dec)
    gw_x, gw_y = map(gw_ra,gw_dec)

    parallels = np.arange(-90,90,2)
    map.drawparallels(parallels,labels=[True,False,False,False], labelstyle="+/-")
    meridians = np.arange(-180,180,5)
    map.drawmeridians(meridians,labels=[False,False,False,True], labelstyle="+/-")

    ### Position Graph
    dot_alpha = 0.025

    map.scatter(gw_x, gw_y, marker='.', c = alphas, zorder = 10, alpha = 1)
    cbar = plt.colorbar()
    map.scatter(p_x, p_y, marker='.', color = 'orange', zorder = 10, alpha = dot_alpha, label = "PanSTARRS1\n" + "{:,}".format(len(p_x)) + " Galaxies")
    map.scatter(g_x, g_y, marker='.', color='aqua', zorder=10, alpha=dot_alpha, label = "GLADE\n" + "{:,}".format(len(g_x)) + " Galaxies")
    plt.xlabel("Right Ascension", labelpad=20)
    plt.ylabel("Declination", labelpad=30)
    # cbar = plt.colorbar()
    cbar.set_label("Probability")
    leg = plt.legend(loc=2, prop={'size': 6})
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    plt.title("PanSTARRS1, GLADE, & GW190814 Positions")
    # plt.savefig("images/CSV_TESTING_Zoomed Map_inProgress.png", bbox_inches = "tight", dpi = 300)

    ## Red shift Histogram - PS1
    plt.figure(3)
    # plt.hist([x for x in PS1_z if 0<x<=1], bins=20)
    plt.hist(PS1_z, bins=20)
    plt.title("Histogram of Redshifts from PanSTARRS1")
    plt.xlabel("Photometric Red Shift")
    plt.ylabel("Frequency of Galaxies")
    plt.savefig("images/CSV_TESTING_Z Hist PS1_inProgress.png", bbox_inches="tight", dpi=300)

    ## Red shift Histogram - GLADE
    plt.figure(4)
    # plt.hist([x for x in GLADE_z if 0 < x <= 1], bins=20)
    plt.hist(GLADE_z, bins=20)
    plt.title("Histogram of Redshifts from GLADE")
    plt.xlabel("Spectroscopic Red Shift")
    plt.ylabel("Frequency of Galaxies")
    plt.savefig("images/CSV_TESTING_Z Hist GLADE_inProgress.png", bbox_inches="tight", dpi=300)

    ## PS1 ps_score
    plt.figure(6)
    plt.hist(PS1_ps_score, bins=20)
    plt.title("Histogram of ps_score From PanSTARRS1")
    plt.xlabel("ps_score (0 is extended source, 1 is point source)")
    plt.ylabel("Frequency of ps_score")
    # plt.savefig("images/CSV_TESTING_ps_score Histogram.png", bbox_inches="tight", dpi=300)


    print("Finished Plotting " + str(datetime.now() - start))

if write_csv:
    start = datetime.now()
    print("Starting Write to CSV - " + str(start.time()))
    with open("local_data/PS1_new_limit.csv", mode='w') as PS1_file:
        PS1_csv = csv.writer(PS1_file, delimiter = ',')
        PS1_csv.writerow(PS1_new_columns)
        for PS1_row in PS1:
            PS1_csv.writerow(PS1_row)

    with open("local_data/GLADE_new_limit.csv", mode='w') as GLADE_file:
        GLADE_csv = csv.writer(GLADE_file, delimiter = ',')
        GLADE_csv.writerow([GLADE_columns[0], GLADE_columns[1], GLADE_columns[2], GLADE_columns[8], GLADE_columns[9], GLADE_columns[13], GLADE_columns[14], GLADE_columns[15], "z_err", GLADE_columns[16]])
        # GLADE_csv.writerow([GLADE_columns[x] for x in GLADE_good_columns])
        for GLADE_row in GLADE:
            GLADE_csv.writerow([GLADE_row[0], GLADE_row[1], GLADE_row[2], GLADE_row[8], GLADE_row[9], GLADE_row[13], GLADE_row[14], GLADE_row[15], float(10**-4), GLADE_row[16]])
            # GLADE_csv.writerow([GLADE_row[x] for x in GLADE_good_columns])

    print("Finished Writting to CSV: " + str(datetime.now() - start))

print("Finished Entire Script: " + str(datetime.now() - script_start))