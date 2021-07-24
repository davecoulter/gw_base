import matplotlib
matplotlib.use("Agg")

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

load_GLADE_db = True
load_PanSTARRS1_db = True
cross_match = True
plot = True
load_GW = True
gw_zone = True

### GRAVITATIONAL WAVE
if load_GW:
    # hpx = hp.read_map("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0")
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
            print(perc_now, "%")
            perc_now = perc_now + 10
        if credible_levels[i] <= 0.90:
            gw_bools[i] = True

    alphas = 1 - credible_levels[indexes[gw_bools]]
    theta, phi = hp.pix2ang(nside, indexes[gw_bools])
    gw_ra = [np.rad2deg(x) for x in phi]
    gw_dec = [np.rad2deg(0.5 * np.pi - x) for x in theta]
    dist_mean = dist_mean[gw_bools]
    dist_std = dist_std[gw_bools]

    print("Len  GW =", len(gw_ra))
    print("GW Percentage in 90% Confidence of Full File =", 100 * (len(gw_ra) / npix), "%")
    print("Time to Complete: " + str(datetime.now() - start))

    ### SHOW GW RA AND DEC LIMITS
    print("First Blob")
    print("RA Limits = [" + str(min([x for x in gw_ra if x <= 20])) + ", " + str(max([x for x in gw_ra if x <= 20])) + "]")
    print("Dec Limits = [" + str(min([x for x in gw_dec if x >= -30])) + ", " + str(max([x for x in gw_dec if x >= -30])) + "]")
    print("Second Blob")
    print("RA Limits = [" + str(min([x for x in gw_ra if x >= 20])) + ", " + str(max([x for x in gw_ra if x >= 20])) + "]")
    print("Dec Limits = [" + str(min([x for x in gw_dec if x <= -30])) + ", " + str(max([x for x in gw_dec if x <= -30])) + "]")

### GALAXIES

### PANSTARRS1 - DB
if load_PanSTARRS1_db:
    start = datetime.now()
    print("Loading PS1 Galaxy " + str(start.time()))
    db_query = '''
    SELECT ra,PS1_Galaxy_v4.dec, z_phot, class, objID, gMeanKronMag, ps_score FROM PS1_Galaxy_v4;
    '''
    PS1 = query_db([db_query])[0]
    PS1_ra = [x[0] for x in PS1]
    PS1_dec = [x[1] for x in PS1]
    PS1_z = [x[2] for x in PS1]
    PS1_class = [x[3] for x in PS1]
    PS1_objid = [x[4] for x in PS1]
    PS1_G_band = [x[5] for x in PS1]
    PS1_ps_score = [x[6] for x in PS1]
    print("Those with Galaxy Tag: " + str(len([x for x in PS1_class if x == "GALAXY"])/len(PS1_class) * 100) + "%")
    print("Len PS1 = " + str(len(PS1_ra)))
    print("Time to Complete: " + str(datetime.now() - start))

### GLADE
if load_GLADE_db:
    start = datetime.now()
    print("Loading GLADE Galaxy " + str(start.time()))
    db_query = '''
    SELECT RA, _DEC, z, B FROM GalaxyDistance2 
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
    print("Len GLADE = " + str(len(GLADE_ra)))
    print("Time to Complete: " + str(datetime.now() - start))

### Crossmatch of GLADE & PS1
if cross_match:
    start = datetime.now()
    print("Starting Cross Match - " + str(start.time()))
    nums = len(GLADE_ra)
    cross_match = [[] for i in range(nums)]
    dist_limit = 1 / 3600 # 1/3600 is 1 arcsecond in degrees
    print("Distance Limit: " + str(dist_limit * 3600) + " arcseconds")
    for i in range(nums):
        right = GLADE_ra[i] + dist_limit
        left = GLADE_ra[i] - dist_limit
        lower = GLADE_dec[i] - dist_limit
        upper = GLADE_dec[i] + dist_limit
        local_PS1_galaxies = [x for x in range(len(PS1_objid)) if PS1_ra[x] <= right and PS1_ra[x] >= left and PS1_dec[x] >= lower and PS1_dec[x] <= upper]
        for ps1_index in local_PS1_galaxies:
            dist = np.sqrt(((GLADE_ra[i] - PS1_ra[ps1_index]) ** 2) + ((GLADE_dec[i] - PS1_dec[ps1_index]) ** 2))
            if dist < dist_limit:
                cross_match[i] = cross_match[i] + [PS1_objid[ps1_index]]
    print("Finished Cross Match: " + str(datetime.now() - start))
    print("Number of Cross Matches: " + str(len([x for x in cross_match if len(x) > 0])))
    print("Number of Cross Matches above 1: " + str(len([x for x in cross_match if len(x) > 1])))

    start = datetime.now()
    print("Limiting PS1 Data and Making Histogram of Differences - " + str(start.time()))
    GLADE_indexs = [x for x in range(len(cross_match)) if len(cross_match[x]) > 0]
    only_cross = [x for x in cross_match if len(x) > 0]
    PS1_index = [[x for x in range(len(PS1_objid)) if PS1_objid[x] == only_cross[y][0]][0] for y in range(len(only_cross))]
    band_diff = [abs(PS1_G_band[PS1_index[x]] - GLADE_B_band[GLADE_indexs[x]]) for x in range(len(GLADE_indexs)) if (type(PS1_G_band[PS1_index[x]]) == float and type(GLADE_B_band[GLADE_indexs[x]]) == float)]
    z_diff = [abs(PS1_z[PS1_index[x]] - GLADE_z[GLADE_indexs[x]]) for x in range(len(GLADE_indexs)) if (type(PS1_z[PS1_index[x]]) == float and type(GLADE_z[GLADE_indexs[x]]) == float)]
    cross_GLADE_z = [GLADE_z[GLADE_indexs[x]] for x in range(len(GLADE_indexs)) if type(GLADE_z[GLADE_indexs[x]]) == float]
    cross_PS1_z = [PS1_z[PS1_index[x]] for x in range(len(PS1_index)) if type(PS1_z[PS1_index[x]]) == float and PS1_z[PS1_index[x]] >=0 and PS1_z[PS1_index[x]] <= 1]

    plt.figure(1)
    plt.hist([x for x in band_diff if x < 100], bins=20)
    plt.title("Histogram of B-band and G-band difference of GLADE and PS1 Cross Match")
    plt.xlabel("Band Difference (mags)")
    plt.ylabel("Frequency of difference")
    plt.savefig("images/Crossmatch Band diff.png", bbox_inches="tight", dpi=300)

    plt.figure(5)
    print("GLADE cross Z Mean: " + str(np.mean(cross_GLADE_z)))
    print("PS1 cross Z Mean: " + str(np.mean(cross_PS1_z)))
    plt.hist([x for x in z_diff if x < 0.25], bins=20)
    plt.title("Histogram of Red Shift difference of GLADE and PS1 Cross Match")
    plt.xlabel("Red Shift Difference")
    plt.ylabel("Frequency of difference")
    plt.savefig("images/Red Shift diff.png", bbox_inches="tight", dpi=300)

    # Get Rid of Cross Match PS1 Galaxies
    PS1_ra = [PS1_ra[x] for x in range(len(PS1_ra)) if x not in PS1_index]
    PS1_dec = [PS1_dec[x] for x in range(len(PS1_dec)) if x not in PS1_index]
    PS1_z = [PS1_z[x] for x in range(len(PS1_z)) if x not in PS1_index]
    PS1_class = [PS1_class[x] for x in range(len(PS1_class)) if x not in PS1_index]
    PS1_objid = [PS1_objid[x] for x in range(len(PS1_objid)) if x not in PS1_index]
    PS1_G_band = [PS1_G_band[x] for x in range(len(PS1_G_band)) if x not in PS1_index]
    PS1_ps_score = [PS1_ps_score[x] for x in range(len(PS1_ps_score)) if x not in PS1_index]
    print("New PS1 Len: " + str(len(PS1_ra)))
    print("Finish Limiting PS1 Data: " + str(datetime.now() - start))

### ONLY IN GW REGION
if gw_zone:
    start = datetime.now()
    print("Start Limit to GW Zone - " + str(start.time()))
    PS1_bools = np.ones(len(PS1_ra), dtype=bool)
    GLADE_bools = np.ones(len(GLADE_ra), dtype=bool)

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

    PS1_ra = [PS1_ra[x] for x in range(len(PS1_ra)) if PS1_bools[x]]
    PS1_dec = [PS1_dec[x] for x in range(len(PS1_dec)) if PS1_bools[x]]
    PS1_z = [PS1_z[x] for x in range(len(PS1_z)) if PS1_bools[x]]
    PS1_class = [PS1_class[x] for x in range(len(PS1_class)) if PS1_bools[x]]
    PS1_objid = [PS1_objid[x] for x in range(len(PS1_objid)) if PS1_bools[x]]
    PS1_G_band = [PS1_G_band[x] for x in range(len(PS1_G_band)) if PS1_bools[x]]
    PS1_ps_score = [PS1_ps_score[x] for x in range(len(PS1_ps_score)) if PS1_bools[x]]
    print("New Len PS1: " + str(len(PS1_ra)))

    GLADE_ra = [GLADE_ra[x] for x in range(len(GLADE_ra)) if GLADE_bools[x]]
    GLADE_dec = [GLADE_dec[x] for x in range(len(GLADE_dec)) if GLADE_bools[x]]
    GLADE_z = [GLADE_z[x] for x in range(len(GLADE_z)) if GLADE_bools[x]]
    GLADE_B_band = [GLADE_B_band[x] for x in range(len(GLADE_B_band)) if GLADE_bools[x]]
    print("New Len GLADE: " + str(len(GLADE_ra)))
    print("Finished Limit to GW Zone: " + str(datetime.now() - start))

### SKY MAP/HISTOGRAMS
if plot:
    start = datetime.now()
    print("Start Plotting " + str(start.time()))
    plt.figure(2)
    map = Basemap(width=3*(10**6),height=3*(10**6)*0.75,projection='lcc',
                resolution='c',lat_0=np.mean(gw_dec),lon_0=np.mean(gw_ra))
    # map = Basemap(projection='lcc',
    #               resolution='c', lat_0=np.mean(gw_dec), lon_0=np.mean(gw_ra))
    # map.drawmapboundary(fill_color='aqua')
    # map.fillcontinents(color='coral',lake_color='aqua')
    # map.drawcoastlines()

    p_x, p_y = map(PS1_ra,PS1_dec)
    # d_x, d_y = map(d_ra,d_dec)
    g_x, g_y = map(GLADE_ra,GLADE_dec)
    gw_x, gw_y = map(gw_ra,gw_dec)
    ### map.scatter(longitude, latitude)
    ### longitude = ra, latitude = dec

    parallels = np.arange(-90,90,2)
    map.drawparallels(parallels,labels=[True,False,False,False], labelstyle="+/-")
    meridians = np.arange(-180,180,5)
    map.drawmeridians(meridians,labels=[False,False,False,True], labelstyle="+/-")

    ### Position Graph
    dot_alpha = 0.05
    map.scatter(p_x, p_y, marker='.', color = 'm', zorder = 10, alpha = dot_alpha, label = "PanSTARRS1\n" + "{:,}".format(len(p_x)) + " Galaxies")
    # map.scatter(d_x, d_y, marker='.', color='r', zorder=10, alpha=dot_alpha, label = "DES\n" + "{:,}".format(len(p_x)) + " Galaxies")
    map.scatter(g_x, g_y, marker='.', color='aqua', zorder=10, alpha=dot_alpha, label = "GLADE\n" + "{:,}".format(len(g_x)) + " Galaxies")
    map.scatter(gw_x, gw_y, marker='.', c = alphas, zorder = 10, alpha = 1)
    plt.xlabel("Right Ascension", labelpad=20)
    plt.ylabel("Declination", labelpad=30)
    cbar = plt.colorbar()
    cbar.set_label("Probability")
    leg = plt.legend(loc=2, prop={'size': 6})
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    plt.title("PanSTARRS1, GLADE, & GW190814 Positions")
    plt.savefig("images/Zoomed Map_inProgress.png", bbox_inches = "tight", dpi = 300)


    ## PS1 ps_score
    plt.figure(6)
    plt.hist(PS1_ps_score, bins=20)
    plt.title("Histogram of ps_score From PanSTARRS1 with Galaxy Restrictions")
    plt.xlabel("ps_score (0 is galaxy, 1 is star)")
    plt.ylabel("Frequency of ps_Score")
    plt.savefig("images/ps_score Histogram.png", bbox_inches="tight", dpi=300)

    ## Red shift Histogram - PS1
    plt.figure(3)
    # plt.hist([x for x in PS1_z if 0<x<=1], bins=20)
    plt.hist([x for x in PS1_z if -0.5<x<=1], bins=20)
    plt.title("Histogram of Redshifts from PanSTARRS1")
    plt.xlabel("Photometric Red Shift")
    plt.ylabel("Frequency of Galaxies")
    plt.savefig("images/Z Hist PS1_inProgress.png", bbox_inches = "tight", dpi = 300)

    ## Red shift Histogram - PS1
    plt.figure(4)
    # plt.hist([x for x in GLADE_z if 0 < x <= 1], bins=20)
    plt.hist([x for x in GLADE_z if 0 < x <= 1.25], bins=20)
    plt.title("Histogram of Redshifts from GLADE")
    plt.xlabel("Spectroscopic Red Shift")
    plt.ylabel("Frequency of Galaxies")
    plt.savefig("images/Z Hist GLADE_inProgress.png", bbox_inches="tight", dpi=300)

    print("Finished Plotting " + str(datetime.now() - start))