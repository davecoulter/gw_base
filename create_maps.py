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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import scipy.stats
c = 299792.458

script_start = datetime.now()

def gaussian_func(x, norm, mean, sigma):
    return norm * ((1.0/(sigma*np.sqrt(2.0*np.pi))) * (np.e**(-1.0*((x - mean)**2)/(2.0*sigma**2))))
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

thing = lambda x: gaussian_func(x, 2, 3, 4)
thing2 = lambda x: 2*gaussian(x=x,mu=3, sig=4)
print(thing(np.array([1,2,3,4])))
print(thing2(np.array([1,2,3,4])))
x = np.arange(20,81,1)
y = gaussian_func(x,3,50,25)
plt.figure(99)
plt.plot(x,y)
plt.savefig("images/TESTING Gauss Func.png", bbox_inches = "tight", dpi = 300)

### GRAVITATIONAL WAVE
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
print("100%")

alphas = 1 - credible_levels[indexes[gw_bools]]
theta, phi = hp.pix2ang(nside, indexes[gw_bools])
gw_ra = [np.rad2deg(x) for x in phi]
gw_dec = [np.rad2deg(0.5 * np.pi - x) for x in theta]
# dist_mean = dist_mean[gw_bools]
# dist_std = dist_std[gw_bools]

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


# ### Galaxies per GW Pixel/Get Distance and Probs
# start = datetime.now()
# print("Starting Galaxy per Pixel - " + str(start.time()))
# PS1_dist = np.zeros(len(PS1_ra), dtype = float)
# PS1_dist_err = np.zeros(len(PS1_ra), dtype = float)
# PS1_probs = np.zeros(len(PS1_ra), dtype = float)
# GLADE_dist = np.zeros(len(GLADE_ra), dtype = float)
# GLADE_dist_err = np.zeros(len(GLADE_ra), dtype = float)
# GLADE_probs = np.zeros(len(GLADE_ra), dtype = float)
# galaxies_in_pixel = np.zeros(npix, dtype = int)
# for i in range(len(PS1_ra)):
#     phi = np.deg2rad(PS1_ra[i])
#     theta = 0.5 * np.pi - np.deg2rad(PS1_dec[i])
#     this_pix = hp.ang2pix(nside, theta, phi)
#     PS1_dist[i] = dist_mean[this_pix]
#     PS1_dist_err[i] = dist_std[this_pix]
#     PS1_probs[i] = 1 - credible_levels[this_pix]
#     galaxies_in_pixel[this_pix] = galaxies_in_pixel[this_pix] + 1
# for i in range(len(GLADE_ra)):
#     phi = np.deg2rad(GLADE_ra[i])
#     theta = 0.5 * np.pi - np.deg2rad(GLADE_dec[i])
#     this_pix = hp.ang2pix(nside, theta, phi)
#     GLADE_dist[i] = dist_mean[this_pix]
#     GLADE_dist_err[i] = dist_std[this_pix]
#     GLADE_probs[i] = 1 - credible_levels[this_pix]
#     galaxies_in_pixel[this_pix] = galaxies_in_pixel[this_pix] + 1
# print("Pixel with 1 Galaxy: " + str(len([x for x in galaxies_in_pixel if x == 1])) + ", with 2 Galaxies: " + str(len([x for x in galaxies_in_pixel if x == 2])) + ", with 3+ Galaxies: " + str(len([x for x in galaxies_in_pixel if x >= 3])))
# print("Finished Galaxy per Pixel/Get Distance and Probs: " + str(datetime.now() - start))
#
# #region### H0 Calculations - Mine
# hubble_const = np.zeros(len(PS1_ra) + len(GLADE_ra), dtype=float)
# hubble_const_probs = np.zeros(len(hubble_const), dtype = float)
# hubble_const_err = np.zeros(len(hubble_const), dtype = float)
# index_hubble = 0
# for i in range(len(PS1_ra)):
#     phi = np.deg2rad(PS1_ra[i])
#     theta = 0.5 * np.pi - np.deg2rad(PS1_dec[i])
#     this_pix = hp.ang2pix(nside, theta, phi)
#     hubble_const[index_hubble] = c*PS1_z[i]/PS1_dist[i]
#     hubble_const_err[index_hubble] = hubble_const[index_hubble] * np.sqrt(((PS1_z[i]/PS1_z_err[i])**2) + ((PS1_dist[i]/PS1_dist_err[i])**2))
#     hubble_const_probs[index_hubble] = (PS1_probs[i]/galaxies_in_pixel[this_pix])/(PS1_z[i]**3)
#     index_hubble = index_hubble + 1
# for i in range(len(GLADE_ra)):
#     phi = np.deg2rad(GLADE_ra[i])
#     theta = 0.5 * np.pi - np.deg2rad(GLADE_dec[i])
#     this_pix = hp.ang2pix(nside, theta, phi)
#     hubble_const[index_hubble] = c*GLADE_z[i]/GLADE_dist[i]
#     hubble_const_err[index_hubble] = hubble_const[index_hubble] * np.sqrt(((GLADE_z[i] / GLADE_z_err[i]) ** 2) + ((GLADE_dist[i] / GLADE_dist_err[i]) ** 2))
#     hubble_const_probs[index_hubble] = (GLADE_probs[i]/galaxies_in_pixel[this_pix])/(GLADE_z[i]**3)
#     index_hubble = index_hubble + 1
# # Calculate Hubble Const
# # H0_min = 20
# # H0_max = 300
# # # H0_gauss_funcs = [lambda x: hubble_const_probs[i]* gaussian(x=x, mu=hubble_const[i], sig=hubble_const_err[i]) for i in range(len(hubble_const))]
# # H0_gauss_funcs = [lambda x: gaussian_func(x=x, norm=hubble_const_probs[i], mean=hubble_const[i], sigma=hubble_const_err[i]) for i in range(len(hubble_const)) if H0_min <= hubble_const[i] <= H0_max]
# # print("Len Gauss Funcs: " + str(len(H0_gauss_funcs)))
# # print(type(H0_gauss_funcs[0]))
# # H0_input = np.arange(start=H0_min, stop=H0_max+1, step=5, dtype=float)
# # H0_input_probs = np.zeros(len(H0_input))
# # for i in range(len(H0_input)):
# #     H0_input_probs[i] = sum([y(H0_input[i]) for y in H0_gauss_funcs])
# # integrate_sum = 0
# # for i in range(len(H0_input_probs)-1):
# #     integrate_sum = integrate_sum + ((H0_input_probs[i] + H0_input_probs[i+1])/2)*(H0_input[i+1] - H0_input[i])
# # H0_input_probs = H0_input_probs/integrate_sum
# # H0_input_probs = np.array([sum([y(x) for y in H0_gauss_funcs]) for x in H0_input])
#
#
# ### Dave H0 Calculation
# print("\n\nDoing Dave H0 Calculation")
# H0_dist = [scipy.stats.norm(loc = hubble_const[x], scale = hubble_const_err[x]) for x in range(len(hubble_const)) if -0 <= hubble_const[x] <= 140]
# limit_hubble_const_probs = [hubble_const_probs[x] for x in range(len(hubble_const_probs))if -0 <= hubble_const[x] <= 140]
# H0_input_dave = np.linspace(-0, 140, 500)
# H0_phot_sum = np.zeros(len(H0_input_dave))
# # for i in range(len(H0_dist)):
# #     dist = H0_dist[i]
# #     prob = hubble_const_probs[i]
# #     result = prob * dist.pdf(H0_input_dave)
# #     # H0_phot_sum[i] = result
# #     if i <= 5:
# #         print(len(result))
# # # print(H0_phot_sum)
#
# for i in range(len(H0_input_dave)):
#     # print(len(H0_dist))
#     # print(len(hubble_const_probs))
#     # print(len(H0_phot_sum))
#     # print(len(H0_input_dave))
#     print(i)
#     result = [limit_hubble_const_probs[x] * H0_dist[x].pdf(H0_input_dave[i]) for x in range(len(H0_dist))]
#     result = [float(x) if x > 10**-10 else 0.0 for x in result]
#     # print("Result:",result)
#     H0_phot_sum[i] = sum(result)
#     # print("H0 Prob Sum:", H0_phot_sum[i])
#
# # Normalize
# integrate_sum = 0
# for i in range(len(H0_phot_sum)-1):
#     integrate_sum = integrate_sum + ((H0_phot_sum[i] + H0_phot_sum[i+1])/2)*(H0_input_dave[i+1] - H0_input_dave[i])
# H0_phot_sum = H0_phot_sum/integrate_sum
# # H0_mean = sum([H0_input_dave[x]*H0_phot_sum[x] for x in range(len(H0_input_dave))])/sum(H0_phot_sum)
#
# #Get 16%, 50%, 84% confidence interval
# interval_16 = -1
# interval_50 = -1
# interval_84 = -1
# running_prob = 0
# dx = H0_input_dave[1] - H0_input_dave[0]
# for i in range(len(H0_input_dave)):
#     running_prob = running_prob + H0_phot_sum[i]*dx
#     if interval_16 == -1 and running_prob >= 0.16:
#         interval_16 = i
#     if interval_50 == -1 and running_prob >= 0.50:
#         interval_50 = i
#     if interval_84 == -1 and running_prob >= 0.84:
#         interval_84 = i
# frac_measurement = 100*(H0_input_dave[interval_84] - H0_input_dave[interval_16])/(2*H0_input_dave[interval_50])
#
# plt.figure(23)
# plt.plot(H0_input_dave,H0_phot_sum, color = "blue")
# plt.xlabel("Hubble Constant (km/s/Mpc)")
# plt.ylabel("Probability of H0")
# plt.title("H0 PDF")
# plt.axvline(x = H0_input_dave[interval_16], color = "red", linestyle = "dashed")
# plt.axvline(x = H0_input_dave[interval_50], color = "red", label = r"H0 = " + "%0.2f$^{+%0.2f}_{-%0.2f}$ (%0.2f%%)" % (H0_input_dave[interval_50], H0_input_dave[interval_84]-H0_input_dave[interval_50],H0_input_dave[interval_50]-H0_input_dave[interval_16], frac_measurement))
# # plt.axvline(x = H0_input_dave[interval_50], color = "red", label = r"H0 = " + "${}^+{}_-{}$".format(H0_input_dave[interval_50], H0_input_dave[interval_84]-H0_input_dave[interval_50],H0_input_dave[interval_50]-H0_input_dave[interval_16]))
# plt.axvline(x = H0_input_dave[interval_84], color = "red", linestyle = "dashed")
# plt.legend()
# plt.savefig("images/TESTINNG DAVE H0 PDF.png", bbox_inches = "tight", dpi = 300)
# #endregion




### SKY MAP/HISTOGRAMS
start = datetime.now()
print("Start Plotting " + str(start.time()))
fig_size = 4
fig = plt.figure(2, figsize=(fig_size,fig_size*(3/4)))
ax = fig.add_subplot(111)
# map = Basemap(width=2.5*(10**6),height=2.5*(10**6)*0.75,projection='lcc', resolution='c',lat_0=-28,lon_0=17.5)
ra_coords = [9,26]
dec_coords = [-35,-21]
map = Basemap(llcrnrlon=ra_coords[1],llcrnrlat=dec_coords[0], urcrnrlon=ra_coords[0], urcrnrlat=dec_coords[1],projection='lcc', resolution='c',lat_0=(dec_coords[1] + dec_coords[0])/2,lon_0=(ra_coords[1] + ra_coords[0])/2)
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
dot_alpha = 0.00
# # map.scatter(d_x, d_y, marker='.', color='r', zorder=10, alpha=dot_alpha, label = "DES\n" + "{:,}".format(len(p_x)) + " Galaxies")
# map.scatter(p_x, p_y, marker='.', color = 'hotpink', zorder = 10, alpha = dot_alpha, label = "PanSTARRS1\n" + "{:,}".format(len(p_x)) + " Galaxies")
# map.scatter(g_x, g_y, marker='.', color='aqua', zorder=10, alpha=dot_alpha, label = "GLADE\n" + "{:,}".format(len(g_x)) + " Galaxies")
map.scatter(gw_x, gw_y, marker='.', c = alphas, zorder = 10, alpha = 1)
cbar = plt.colorbar()
cbar.set_label("Probability")
cbar.set_ticks(np.linspace(0.0,1.0,11))
# map.scatter(p_x, p_y, marker = 's', facecolors='none', edgecolors = 'hotpink', zorder = 10, alpha = dot_alpha, label = "PanSTARRS1\n" + "{:,}".format(len(p_x)) + " Galaxies")
# map.scatter(g_x, g_y, marker = '^', facecolors='none', edgecolors ='hotpink', zorder=10, alpha=dot_alpha, label = "GLADE\n" + "{:,}".format(len(g_x)) + " Galaxies")
plt.xlabel("Right Ascension", labelpad=20)
plt.ylabel("Declination", labelpad=30)
# leg = plt.legend(loc=2, prop={'size': 6})
# for lh in leg.legendHandles:
#     lh.set_alpha(1)
# # plt.xlim(10,10)
# plt.ylim(-34,-22)
plt.title("GW190814 Positions")
plt.savefig("images/Zoomed Map_inProgress.png", bbox_inches = "tight", dpi = 300)





## PS1 ps_score
plt.figure(6)
plt.hist(PS1_ps_score, bins=20)
plt.title("Histogram of ps_score From PanSTARRS1 with Galaxy Restrictions")
plt.xlabel("ps_score (0 is galaxy, 1 is star)")
plt.ylabel("Frequency of ps_Score")
# plt.savefig("images/ps_score Histogram.png", bbox_inches="tight", dpi=300)

## Red shift Histogram - PS1
plt.figure(3)
# plt.hist([x for x in PS1_z if 0<x<=1], bins=20)
plt.hist(PS1_z, bins=20)
plt.title("Histogram of Redshifts from PanSTARRS1")
plt.xlabel("Photometric Red Shift")
plt.ylabel("Frequency of Galaxies")
# plt.savefig("images/Z Hist PS1_inProgress.png", bbox_inches = "tight", dpi = 300)

## Red shift Histogram - GLADE
plt.figure(4)
# plt.hist([x for x in GLADE_z if 0 < x <= 1], bins=20)
plt.hist(GLADE_z, bins=20)
plt.title("Histogram of Redshifts from GLADE")
plt.xlabel("Spectroscopic Red Shift")
plt.ylabel("Frequency of Galaxies")
# plt.savefig("images/Z Hist GLADE_inProgress.png", bbox_inches="tight", dpi=300)

## H0 PDF
# plt.figure(8)
# # plt.hist(hubble_const,bins = 20, weights=hubble_const_probs, density=False)
# plt.hist([hubble_const[x] for x in range(len(hubble_const)) if 20 <= hubble_const[x] <= 150],bins = 20, weights=[hubble_const_probs[x] for x in range(len(hubble_const)) if 20 <= hubble_const[x] <= 150], density=True, log = True)
# # plt.yscale('log', nonposy='clip')
# plt.title("Histogram of Hubble Constant")
# plt.xlabel("Hubble Constant (km/s/Mpc)")
# plt.ylabel("Probability of H0")
# # plt.savefig("images/HO PDF_inProgress.png", bbox_inches="tight", dpi=300)
# plt.savefig("images/HO PDF_inProgress.png", bbox_inches="tight", dpi=300)

# plt.figure(8)
# print(H0_input_probs)
# plt.plot(H0_input, H0_input_probs)
# # plt.yscale('log')
# plt.title("H0 PDF")
# plt.xlabel("Hubble Constant (km/s/Mpc)")
# plt.ylabel("Probability of H0")
# plt.savefig("images/H0 PDF_inProgress.png", bbox_inches="tight", dpi=300)

print("Finished Plotting " + str(datetime.now() - start))


print("Finished Script: " + str(datetime.now() - script_start))