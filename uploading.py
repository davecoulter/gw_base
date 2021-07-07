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



if False:
    map = Basemap(projection='hammer',
                  lat_0=0, lon_0=0)
    map.drawmapboundary(fill_color='aqua')
    map.fillcontinents(color='coral',lake_color='aqua')
    map.drawcoastlines()
    x, y = map(gw_ra,gw_dec)
    map.scatter(x, y, marker='.',color='m', zorder = 10, alpha = 0.05)
    plt.savefig("map_test_2.png")




if False:
    gw_190814 = fits.open("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0")
    gw_190814.info()
    gw_data = gw_190814[1].data
    gw_190814.close()
    for i in np.linspace(0, 12000000, 500, dtype=int):
        print(gw_data[i][0],gw_data[i][1],gw_data[i][2],gw_data[i][3])
    print("Type =", type(gw_data[0]))
    print("Len =",len(gw_data[0]))

### GRAVITATIONAL WAVE
if True:
    hpx = hp.read_map("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0")
    npix = len(hpx)
    print(npix)
    nside = hp.npix2nside(npix)
    theta = []
    phi = []
    credible_levels = find_greedy_credible_levels(hpx)
    alphas = []
    print("Loading GW")
    perc = 5
    perc_now = perc
    num_in_array = npix
    now = datetime.now()
    for i in np.linspace(0, npix-1, num_in_array, dtype=int):
        if i/npix >= perc_now/100:
            print(perc_now,"%")
            perc_now = perc_now + perc
        if credible_levels[i] <= 0.90:
            alphas = alphas + [1 - credible_levels[i]]
            theta_0, phi_0 = hp.pix2ang(nside, i)
            theta = theta + [theta_0]
            phi = phi + [phi_0]
    gw_ra = [np.rad2deg(x) for x in phi]
    gw_dec = [np.rad2deg(0.5*np.pi - x) for x in theta]

    rgba_colors = np.zeros((len(gw_ra), 4))
    rgba_colors[:, 2] = 1.0
    rgba_colors[:, 3] = alphas

    print("Len  GW =", len(gw_ra))
    print("GW Percentage of Full File =", 100*(num_in_array / npix),"%")
    print("Time to Complete: " + str(datetime.now() - now))



    # map = Basemap(projection='hammer',
    #               lat_0=0, lon_0=0)
    # map.drawmapboundary(fill_color='aqua')
    # map.fillcontinents(color='coral', lake_color='aqua')
    # map.drawcoastlines()
    # x, y = map(gw_ra, gw_dec)
    # map.scatter(x, y, marker='.', zorder=10, color = 'm', alpha = 1)
    # plt.savefig("map_test_4.png")

### GALAXY
if True:
    total_file = 529270
    num_in_array = int(total_file*0.01)
    array_index = np.linspace(0, total_file, num_in_array, dtype=int)
    with open("GW190814_allData_v1_pjquinonez.csv", mode='r') as csv_file:
        galaxy_code = csv.DictReader(csv_file)
        ra = []
        dec = []
        z_phot = []
        i = 0
        perc = 5
        perc_now = perc
        print("Loading Galaxy")
        now = datetime.now()
        for row in galaxy_code:
            if i/total_file >= perc_now/100:
                print(perc_now, "%")
                perc_now = perc_now + perc
            if i in array_index:
                ra = ra + [float(row["ra"])]
                dec = dec + [float(row["dec"])]
                if -0 < float(row["z_phot"]) <= 1.0:
                    z_phot = z_phot + [float(row["z_phot"])]

            i = i + 1
    ra = np.array(ra)
    dec = np.array(dec)
    print("Len  Galaxies =",len(ra))
    print("Galaxy Percentage of Full File =",100*(len(ra)/total_file),"%")
    print("Redshift Percentage of Full File =", 100*(len(z_phot) / total_file),"%")
    print("Time to Complete: " + str(datetime.now() - now))

### SKY MAP/HISTOGRAM
if True:
    map2 = Basemap(width=12000000,height=9000000,projection='lcc',
                resolution='c',lat_0=np.mean(dec),lon_0=np.mean(ra))
    map2.drawmapboundary(fill_color='aqua')
    map2.fillcontinents(color='coral',lake_color='aqua')
    map2.drawcoastlines()
    x, y = map2(ra,dec)
    x_2, y_2 = map2(gw_ra,gw_dec)
    # map.scatter(longitude, latitude)
    # latitude = dec, longitude = ra

    parallels = np.arange(-90,90,10.)
    map2.drawparallels(parallels,labels=[True,False,False,False])
    meridians = np.arange(-180,180,20.)
    map2.drawmeridians(meridians,labels=[False,False,False,True])

    map2.scatter(x, y, marker='.', color = 'm', zorder = 10, alpha = 0.05)
    map2.scatter(x_2, y_2, marker='.',color=rgba_colors, zorder = 10, alpha = 1)
    plt.title("Pan-STARRS1 Galaxy & GW190814 Positions")
    plt.savefig("Zoomed Map 2.png")

    plt.figure(5)
    plt.hist(z_phot, bins=20)
    plt.title("Histogram of Redshifts from Pan-STARRS1")
    plt.xlabel("Photometric Red Shift")
    plt.ylabel("Frequency of Galaxies")
    plt.savefig("Z Hist.png", bbox_inches = "tight")
    print("It worked")

if False:
    gw_190814 = fits.open("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0")
    gw_190814.info()
    gw_190814.close()
    gw_data = gw_190814[1].data
    # gw_ra = [x[1] for x in gw_data]
    # gw_dec = [x[2] for x in gw_data]
    gw_ra = []
    gw_dec = []
    for i in np.linspace(0, 12000000, 5000, dtype=int):
        gw_ra = gw_ra + [gw_data[i][1]]
        gw_dec = gw_dec + [[gw_data[i][2]]]
    gw_ra = np.array(gw_ra)
    gw_dec = np.array(gw_dec)