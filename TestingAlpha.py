import astropy.cosmology
import matplotlib.pyplot as plt
import numpy as np
import csv
from datetime import datetime
from database_methods import *
import healpy as hp
from ligo.skymap.postprocess import find_greedy_credible_levels
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
from ligo.skymap import distance
import scipy.stats
import random
from astropy.cosmology import *
import astropy.units as u
import scipy

import os

#region Old Maps, Fits, CSV upload, DES
## DRAW GLOBE
# if False:
#     # set up orthographic map projection with
#     # perspective of satellite looking down at 50N, 100W.
#     # use low resolution coastlines.
#     map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
#     # draw coastlines, country boundaries, fill continents.
#     map.drawcoastlines(linewidth=0.25)
#     map.drawcountries(linewidth=0.25)
#     map.fillcontinents(color='coral',lake_color='aqua')
#     # draw the edge of the map projection region (the projection limb)
#     map.drawmapboundary(fill_color='aqua')
#     # draw lat/lon grid lines every 30 degrees.
#     map.drawmeridians(np.arange(0,360,30))
#     map.drawparallels(np.arange(-90,90,30))
#     # make up some data on a regular lat/lon grid.
#     nlats = 73; nlons = 145; delta = 2.*np.pi/(nlons-1)
#     lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])
#     lons = (delta*np.indices((nlats,nlons))[1,:,:])
#     wave = 0.75*(np.sin(2.*lats)**8*np.cos(4.*lons))
#     mean = 0.5*np.cos(2.*lats)*((np.sin(2.*lats))**2 + 2.)
#     # compute native map projection coordinates of lat/lon grid.
#     x, y = map(lons*180./np.pi, lats*180./np.pi)
#     # contour data over the map.
#     # cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
#     map.scatter([0],[0])
#     plt.title('Galaxy Map')
#     plt.savefig("map_test.png")
#     plt.show()
#
# ### DRAW ENTIRE SKY MAP
# if False:
#     map = Basemap(projection='hammer',
#                   lat_0=0, lon_0=0)
#     map.drawmapboundary(fill_color='aqua')
#     map.fillcontinents(color='coral',lake_color='aqua')
#     map.drawcoastlines()
#     x, y = map(gw_ra,gw_dec)
#     map.scatter(x, y, marker='.',color='m', zorder = 10, alpha = 0.05)
#     plt.savefig("map_test_2.png")
#
#
# ### OPENING HEALPIX AS FITS
# if False:
#     gw_190814 = fits.open("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0")
#     gw_190814.info()
#     gw_data = gw_190814[1].data
#     print(gw_190814[1].columns)
#     gw_190814.close()
#     for i in np.linspace(0, 12000000, 5, dtype=int):
#         print(gw_data[i][0],gw_data[i][1],gw_data[i][2],gw_data[i][3])
#     print("Type =", type(gw_data[0]))
#     print("Len =",len(gw_data[0]))
#
# ### PS1 AS CSV
# if False:
#     total_file = 529269
#     num_in_array = int(total_file*0.01)
#     array_index = np.linspace(0, total_file, num_in_array, dtype=int)
#     with open("local_data/PanSTARRS_allData_v1_pjquinonez.csv", mode='r') as csv_file:
#         galaxy_code = csv.DictReader(csv_file)
#         p_ra = np.zeros(num_in_array)
#         p_dec = np.zeros(num_in_array)
#         p_z_phot = np.zeros(num_in_array)
#         index = 0
#         i = 0
#         perc = 5
#         perc_now = perc
#         print("Loading PanSTARRS Galaxy")
#         now = datetime.now()
#         for row in galaxy_code:
#             if i/total_file >= perc_now/100:
#                 print(perc_now, "%", datetime.now() - now)
#                 perc_now = perc_now + perc
#             if i in array_index:
#                 p_ra[index] = float(row["ra"])
#                 p_dec[index] = float(row["dec"])
#                 p_z_phot[index] = float(row["z_phot"])
#                 index = index + 1
#             i = i + 1
#
#     p_ra = p_ra[:index]
#     p_dec = p_dec[:index]
#     p_z_phot = p_z_phot[:index]
#     print("Len  PanSTARRS =",len(p_ra))
#     print("PanSTARRS Percentage of Full File =",100*(len(p_ra)/total_file),"%")
#     print("PanSTARRS Redshift Percentage of Full File =", 100*(len(p_z_phot) / total_file),"%")
#     print("Time to Complete: " + str(datetime.now() - now))
#
# ### DES AS CSV
# if False:
#     total_file = 1361811
#     num_in_array = int(total_file * 0.01)
#     array_index = np.linspace(0, total_file, num_in_array, dtype=int)
#     with open("local_data/DES_allBands.csv", mode='r') as csv_file:
#         galaxy_code = csv.DictReader(csv_file)
#         d_ra = np.zeros(num_in_array)
#         d_dec = np.zeros(num_in_array)
#         index = 0
#         # d_z_phot = np.zeros(num_in_array)
#         i = 0
#         perc = 5
#         perc_now = perc
#         print("Loading DES Galaxy")
#         now = datetime.now()
#         for row in galaxy_code:
#             if i/total_file >= perc_now/100:
#                 print(perc_now, "%", datetime.now() - now)
#                 perc_now = perc_now + perc
#             if i in array_index:
#                 d_ra[index] = float(row["ALPHAWIN_J2000"])
#                 d_dec[index] = float(row["DELTAWIN_J2000"])
#                 # d_z_phot[index] = float(row["z_phot"])
#                 index = index + 1
#             i = i + 1
#
#     d_ra = d_ra[:index]
#     d_dec = d_dec[:index]
#     # d_z_phot = d_z_phot[:index]
#     print("Len  DES =",len(d_ra))
#     print("DES Percentage of Full File =",100*(len(d_ra)/total_file),"%")
#     # print("DES Redshift Percentage of Full File =", 100*(len(z_phot) / total_file),"%")
#     print("Time to Complete: " + str(datetime.now() - now))
#
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
#endregion

# PS1_columns = "ObjID", "uniquePspsOBid", "raStack", "decStack", "raMean", "decMean", "ra", "dec", "ng", "gMeanPSFMag", "gMeanPSFMagErr", "gMeanKronMag", "gMeanKronMagErr", "gMeanApMag", "gMeanApMagErr", "nr", "rMeanPSFMag", "rMeanPSFMagErr", "rMeanKronMag", "rMeanKronMagErr", "rMeanApMag", "rMeanApMagErr", "ni", "iMeanPSFMag", "iMeanPSFMagErr", "iMeanKronMag", "iMeanKronMagErr", "iMeanApMag", "iMeanApMagErr", "nz", "zMeanPSFMag", "zMeanPSFMagErr", "zMeanKronMag", "zMeanKronMagErr", "zMeanApMag", "zMeanApMagErr", "ny", "yMeanPSFMag", "yMeanPSFMagErr", "yMeanKronMag", "yMeanKronMagErr", "yMeanApMag", "yMeanApMagErr", "gQfPerfect", "rQfPerfect", "iQfPerfect", "zQfPerfect", "yQfPerfect", "qualityFlag", "objInfoFlag", "primaryDetection", "bestDetection", "class", "prob_Galaxy", "prob_Star", "prob_QSO", "z_phot", "z_photErr", "z_phot0", "extrapolation_Photoz", "ps_score"
# PS1_new_columns = ["ObjID", "raMean", "decMean", "class", "prob_Galaxy", "prob_Star", "prob_QSO", "z_phot", "z_photErr", "ps_score"]
# GLADE_columns = "id", "Galaxy_id", "Distance_id", "PGC", "Name_GWGC", "Name_HyperLEDA", "Name_2MASS", "Name_SDSS_DR12", "RA", "_Dec", "Coord", "dist", "dist_err", "z_dist", "z_dist_err", "z", "B", "B_err", "B_abs", "J", "J_err", "H", "H_err", "K", "K_err", "flag1", "flag2", "flag3"
# print("PS1 Columns")
# for i in range(len(PS1_new_columns)):
#     print(i,PS1_new_columns[i])
# print("\nGLADE Columns")
# for i in range(len(GLADE_columns)):
#     print(i,GLADE_columns[i])

# x = np.linspace(0,10,100)
# y = x**2
# x2 = np.linspace(0,10,200)
# y2 = np.interp(x2, x, y)
# print(x2)
# print(y2)
#
# # h0 = np.array([1,1,1,1,2,3,4])
# # probs = np.array([0.25, 0.25,0.25,0.25,3, 0.6, 0.01])
# # plt.figure(1)
# # plt.hist(h0, weights=probs/max(probs), density=False)
# # plt.savefig("images/aaaaaaaaaH0 PDF hist TEST.png", bbox_inches = "tight", dpi = 300)
#
# x = scipy.stats.norm(loc = 5, scale = 1.5).pdf(np.linspace(0,10,1000))
# # print(np.trapz(y2, x=x2))
# print(type(x))
# print(x)


# x = np.linspace(0,10,10)
# y = scipy.stats.norm.pdf(x= x, loc=5, scale = 1.5)*9
# # y = np.append(np.array([0,0,0,0,0,0]),y)
# print(y)
# print("Int:",np.trapz(y,x=x))
# y = y/np.linalg.norm(y)
# print(y)
# print("Int:",np.trapz(y,x=x))
# y = y/np.linalg.norm(y)
# print(y)
# print("Int:",np.trapz(y,x=x))
# # for _ in range(5):
# #     y = y/np.linalg.norm(y)
# # print(y)
# # print("Int:",np.trapz(y,x=x))






# ### PLOTTING REDSHIFT HISTOGRAM
# print("Plotting Final Redshift Histogram")
# fig, ax = plt.subplots()
# ax2 = ax.twinx()
#
# counts, bins, bars = ax.hist(PS1e_z, bins=20, color = "red",label = "PS1 Photometric Redshift (Extended)", rwidth=0.9)
#
# print(len(bins))
# ax.hist(np.append(PS1_z,GLADE_z), bins=bins, color = "orange",label = "PS1 Photometric Redshift", rwidth=0.9)
# ax.hist(GLADE_z, bins=bins, color = "aqua", label = "GLADE Spectroscopic Redshift", rwidth = 0.9)
#
# ax2.plot(x,y)
# ax2.set_ylabel("z^3")
# plt.title("Redshifts of Each Galaxy")
# ax.set_xlabel("Red Shift")
# ax.set_ylabel("Frequency of Galaxies")
# ax.legend(loc="upper left")
# plt.savefig("images/TESTING Redshift Histogram CSV.png", bbox_inches = "tight", dpi = 300)

x = np.linspace(0,10,1000)
y = 3 + 5*x ** 3

(a, b), pcov = scipy.optimize.curve_fit(f=lambda x, a, b: a*(x**3) + b, xdata=x, ydata=y)
print(str(a) + "x^3 + " + str(b))