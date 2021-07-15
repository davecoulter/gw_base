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


### DRAW GLOBE
if False:
    # set up orthographic map projection with
    # perspective of satellite looking down at 50N, 100W.
    # use low resolution coastlines.
    map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='coral',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0,360,30))
    map.drawparallels(np.arange(-90,90,30))
    # make up some data on a regular lat/lon grid.
    nlats = 73; nlons = 145; delta = 2.*np.pi/(nlons-1)
    lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])
    lons = (delta*np.indices((nlats,nlons))[1,:,:])
    wave = 0.75*(np.sin(2.*lats)**8*np.cos(4.*lons))
    mean = 0.5*np.cos(2.*lats)*((np.sin(2.*lats))**2 + 2.)
    # compute native map projection coordinates of lat/lon grid.
    x, y = map(lons*180./np.pi, lats*180./np.pi)
    # contour data over the map.
    # cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
    map.scatter([0],[0])
    plt.title('Galaxy Map')
    plt.savefig("map_test.png")
    plt.show()


### DRAW ENTIRE SKY MAP
if False:
    map = Basemap(projection='hammer',
                  lat_0=0, lon_0=0)
    map.drawmapboundary(fill_color='aqua')
    map.fillcontinents(color='coral',lake_color='aqua')
    map.drawcoastlines()
    x, y = map(gw_ra,gw_dec)
    map.scatter(x, y, marker='.',color='m', zorder = 10, alpha = 0.05)
    plt.savefig("map_test_2.png")



### OPENING HEALPIX AS FITS
if False:
    gw_190814 = fits.open("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0")
    gw_190814.info()
    gw_data = gw_190814[1].data
    print(gw_190814[1].columns)
    gw_190814.close()
    for i in np.linspace(0, 12000000, 5, dtype=int):
        print(gw_data[i][0],gw_data[i][1],gw_data[i][2],gw_data[i][3])
    print("Type =", type(gw_data[0]))
    print("Len =",len(gw_data[0]))

### GRAVITATIONAL WAVE
if True:
    # hpx = hp.read_map("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0")
    prob, distmu, distsigma, distnorm = hp.read_map("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0",field=range(4))
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
            perc_now = perc_now + 5
        if credible_levels[i] <= 0.90:
            gw_bools[i] = True

    alphas = 1 - credible_levels[indexes[gw_bools]]
    theta, phi = hp.pix2ang(nside, indexes[gw_bools])
    gw_ra = [np.rad2deg(x) for x in phi]
    gw_dec = [np.rad2deg(0.5 * np.pi - x) for x in theta]

    print("Len  GW =", len(gw_ra))
    print("GW Percentage in 90% Confidence of Full File =", 100 * (len(gw_ra) / npix), "%")
    print("Time to Complete: " + str(datetime.now() - start))

    ### SHOW GW RA AND DEC LIMITS
    print("RA Limits = [" + str(min(gw_ra)) + ", " + str(max(gw_ra)) + "]")
    print("Dec Limits = [" + str(min(gw_dec)) + ", " + str(max(gw_dec)) + "]")

### GALAXY
### PANSTARRS - CSV
if False:
    total_file = 529269
    num_in_array = int(total_file*0.01)
    array_index = np.linspace(0, total_file, num_in_array, dtype=int)
    with open("local_data/PanSTARRS_allData_v1_pjquinonez.csv", mode='r') as csv_file:
        galaxy_code = csv.DictReader(csv_file)
        p_ra = np.zeros(num_in_array)
        p_dec = np.zeros(num_in_array)
        p_z_phot = np.zeros(num_in_array)
        index = 0
        i = 0
        perc = 5
        perc_now = perc
        print("Loading PanSTARRS Galaxy")
        now = datetime.now()
        for row in galaxy_code:
            if i/total_file >= perc_now/100:
                print(perc_now, "%", datetime.now() - now)
                perc_now = perc_now + perc
            if i in array_index:
                p_ra[index] = float(row["ra"])
                p_dec[index] = float(row["dec"])
                p_z_phot[index] = float(row["z_phot"])
                index = index + 1
            i = i + 1

    p_ra = p_ra[:index]
    p_dec = p_dec[:index]
    p_z_phot = p_z_phot[:index]
    print("Len  PanSTARRS =",len(p_ra))
    print("PanSTARRS Percentage of Full File =",100*(len(p_ra)/total_file),"%")
    print("PanSTARRS Redshift Percentage of Full File =", 100*(len(p_z_phot) / total_file),"%")
    print("Time to Complete: " + str(datetime.now() - now))

### PANSTARRS - DB
if True:
    start = datetime.now()
    print("Loading PS1 Galaxy " + str(start.time()))
    db_query = '''
    SELECT ra,PS1_Galaxy_v3.dec, z_phot FROM PS1_Galaxy_v3;
    '''
    PS1 = query_db([db_query])[0]
    PS1_ra = [x[0] for x in PS1]
    PS1_dec = [x[1] for x in PS1]
    PS1_z = [x[2] for x in PS1]
    print("Len PS1 = " + str(len(PS1_ra)))
    print("Time to Complete: " + str(datetime.now() - start))

### DES
if False:
    total_file = 1361811
    num_in_array = int(total_file * 0.01)
    array_index = np.linspace(0, total_file, num_in_array, dtype=int)
    with open("local_data/DES_allBands.csv", mode='r') as csv_file:
        galaxy_code = csv.DictReader(csv_file)
        d_ra = np.zeros(num_in_array)
        d_dec = np.zeros(num_in_array)
        index = 0
        # d_z_phot = np.zeros(num_in_array)
        i = 0
        perc = 5
        perc_now = perc
        print("Loading DES Galaxy")
        now = datetime.now()
        for row in galaxy_code:
            if i/total_file >= perc_now/100:
                print(perc_now, "%", datetime.now() - now)
                perc_now = perc_now + perc
            if i in array_index:
                d_ra[index] = float(row["ALPHAWIN_J2000"])
                d_dec[index] = float(row["DELTAWIN_J2000"])
                # d_z_phot[index] = float(row["z_phot"])
                index = index + 1
            i = i + 1

    d_ra = d_ra[:index]
    d_dec = d_dec[:index]
    # d_z_phot = d_z_phot[:index]
    print("Len  DES =",len(d_ra))
    print("DES Percentage of Full File =",100*(len(d_ra)/total_file),"%")
    # print("DES Redshift Percentage of Full File =", 100*(len(z_phot) / total_file),"%")
    print("Time to Complete: " + str(datetime.now() - now))


### GLADE
if True:
    start = datetime.now()
    print("Loading GLADE Galaxy " + str(start.time()))
    db_query = '''
    SELECT RA, _DEC, z FROM GalaxyDistance2 
    WHERE 
        RA >= 10.0 AND
        RA <= 25.0 AND
        _DEC <= -20.0 AND
        _DEC >= -35.0;
    '''
    GLADE = query_db([db_query])[0]
    GLADE_ra = [x[0] for x in GLADE]
    GLADE_dec = [x[1] for x in GLADE]
    GLADE_z = [x[2] for x in GLADE]
    print("Len GLADE = " + str(len(GLADE_ra)))
    print("Time to Complete: " + str(datetime.now() - start))


### SKY MAP/HISTOGRAM
if True:
    start = datetime.now()
    print("Start Plotting " + str(start.time()))
    map = Basemap(width=(5*10**6),height=(5*10**6)*0.75,projection='lcc',
                resolution='c',lat_0=np.mean(gw_dec),lon_0=np.mean(gw_ra))
    # map.drawmapboundary(fill_color='aqua')
    # map.fillcontinents(color='coral',lake_color='aqua')
    # map.drawcoastlines()

    p_x, p_y = map(PS1_ra,PS1_dec)
    # d_x, d_y = map(d_ra,d_dec)
    g_x, g_y = map(GLADE_ra,GLADE_dec)
    gw_x, gw_y = map(gw_ra,gw_dec)
    ### map.scatter(longitude, latitude)
    ### longitude = ra, latitude = dec

    parallels = np.arange(-90,90,5)
    map.drawparallels(parallels,labels=[True,False,False,False], labelstyle="+/-")
    meridians = np.arange(-180,180,10)
    map.drawmeridians(meridians,labels=[False,False,False,True], labelstyle="+/-")

    ### Position Graph
    map.scatter(p_x, p_y, marker='.', color = 'm', zorder = 10, alpha = 0.1, label = "PanSTARRS1\n" + str(len(p_x))+ " Galaxies")
    # map.scatter(d_x, d_y, marker='.', color='r', zorder=10, alpha=0.1, label = "DES\n" + str(len(d_x)) + " Galaxies")
    map.scatter(g_x, g_y, marker='.', color='aqua', zorder=10, alpha=0.1, label = "GLADE\n" + str(len(g_x)) + " Galaxies")
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

    ## Red shift Histogram - PS1
    plt.figure(5)
    plt.hist([x for x in PS1_z if 0<x<=1], bins=20)
    plt.title("Histogram of Redshifts from PanSTARRS1")
    plt.xlabel("Photometric Red Shift")
    plt.ylabel("Frequency of Galaxies")
    plt.savefig("images/Z Hist PS1_inProgress.png", bbox_inches = "tight", dpi = 300)

    ## Red shift Histogram - PS1
    plt.figure(6)
    plt.hist([x for x in GLADE_z if 0 < x <= 1], bins=20)
    plt.title("Histogram of Redshifts from GLADE")
    plt.xlabel("Spectroscopic Red Shift")
    plt.ylabel("Frequency of Galaxies")
    plt.savefig("images/Z Hist GLADE_inProgress.png", bbox_inches="tight", dpi=300)

    print("Finished Plotting " + str(datetime.now() - start))


# if False:
#     gw_190814 = fits.open("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0")
#     gw_190814.info()
#     gw_190814.close()
#     gw_data = gw_190814[1].data
#     # gw_ra = [x[1] for x in gw_data]
#     # gw_dec = [x[2] for x in gw_data]
#     gw_ra = []
#     gw_dec = []
#     for i in np.linspace(0, 12000000, 5000, dtype=int):
#         gw_ra = gw_ra + [gw_data[i][1]]
#         gw_dec = gw_dec + [[gw_data[i][2]]]
#     gw_ra = np.array(gw_ra)
#     gw_dec = np.array(gw_dec)