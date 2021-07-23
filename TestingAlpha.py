import matplotlib.pyplot as plt
import numpy as np
import csv
from datetime import datetime
from database_methods import *
import healpy as hp
from ligo.skymap.postprocess import find_greedy_credible_levels
from mpl_toolkits.basemap import Basemap


### Test out database methods
# detector_select = '''
# SELECT RA, _DEC FROM GalaxyDistance2
# WHERE
#     RA >= 10.0 AND
#     RA <= 25.0 AND
#     _DEC <= -20.0 AND
#     _DEC >= -35.0
# LIMIT 5;
# '''
# result = query_db([detector_select])[0]
# for r in result:
#     print("RA: " + str(r[0]) + ", DEC: " + str(r[1]))

### Testing Distance Healpy
# if False:
#     prob, distmu, distsigma, distnorm = hp.read_map("Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0", field = range(4))
#     print(len(prob), len(distmu), len(distsigma), len(distnorm))
#     for i in range(0):
#         print(distmu[i], distsigma[i], distnorm[i])

# ### PS1_Scores
# filename = "local_data/PanSTARRS1_scores_local_data_pjquinonez.csv"
# data = np.loadtxt(filename, delimiter=',', skiprows=1)
# objid = data[:,0]
# scores = data[:,1]
# print(len(scores))
#
# plt.figure(6)
# plt.hist([x for x in GLADE_z if 0 < x <= 1], bins=20)
# plt.title("Histogram of Redshifts from GLADE")
# plt.xlabel("Spectroscopic Red Shift")
# plt.ylabel("Frequency of Galaxies")
# plt.savefig("images/Z Hist GLADE_inProgress.png", bbox_inches="tight", dpi=300)

#region Cross Match
# ### Load in PS1
# print("Loading PS1 Galaxy ")
# db_query = '''
# SELECT ra,PS1_Galaxy_vNewLimits.dec, z_phot, class, prob_Galaxy, objID, gMeanKronMag FROM PS1_Galaxy_vNewLimits;
# '''
# PS1 = query_db([db_query])[0]
# PS1_ra = [x[0] for x in PS1]
# PS1_dec = [x[1] for x in PS1]
# PS1_z = [x[2] for x in PS1]
# PS1_class = [x[3] for x in PS1]
# PS1_prob_Galaxy = [x[4] for x in PS1]
# PS1_objid = [x[5] for x in PS1]
# PS1_G_band = [x[6] for x in PS1]
# print("Unique Objects: " + str(len(np.unique(PS1_objid))))
# print("Len PS1: " + str(len(PS1_objid)))
#
# print()
#
# ### Load in GLADE
# print("Loading GLADE Galaxy")
# db_query = '''
#     SELECT RA, _DEC, z, B FROM GalaxyDistance2
#     WHERE
#
#         RA >= 10.0 AND
#         RA <= 15.0 AND
#         _DEC <= -22.0 AND
#         _DEC >= -28.0;
#     '''
# GLADE = query_db([db_query])[0]
# GLADE_ra = [x[0] for x in GLADE]
# GLADE_dec = [x[1] for x in GLADE]
# GLADE_z = [x[2] for x in GLADE]
# GLADE_B_band = [x[3] for x in GLADE]
# print("Len GLADE = " + str(len(GLADE_ra)))
#
# print()
#
# start = datetime.now()
# print("Starting Cross Match - " + str(start.time()))
# nums = len(GLADE_ra)
# cross_match = [[] for i in range(nums)]
# dist_limit = 1/3600
# print("Distance Limit: " + str(dist_limit*3600) + " arcseconds")
# for i in range(nums):
#     right = GLADE_ra[i] + dist_limit
#     left = GLADE_ra[i] - dist_limit
#     lower = GLADE_dec[i] - dist_limit
#     upper = GLADE_dec[i] + dist_limit
#     local_PS1_galaxies = [x for x in range(len(PS1_objid)) if PS1_ra[x] <= right and PS1_ra[x] >= left and PS1_dec[x] >= lower and PS1_dec[x] <= upper]
#     for ps1_index in local_PS1_galaxies:
#         dist = np.sqrt(((GLADE_ra[i] - PS1_ra[ps1_index])**2) + ((GLADE_dec[i] - PS1_dec[ps1_index])**2))
#         if dist < dist_limit:
#             cross_match[i] = cross_match[i] + [PS1_objid[ps1_index]]
# print("Finished: " + str(datetime.now()-start))
# print("Time for Full GLADE Catalog: " + str((datetime.now()-start)*(len(GLADE_ra)/nums)))
#
# print("Number of Cross Matches: " + str(len([x for x in cross_match if len(x) > 0])))
# print("Number of Cross Matches above 1: " + str(len([x for x in cross_match if len(x) > 1])))
#
# GLADE_ids = [x for x in range(len(cross_match)) if len(cross_match[x]) > 0]
# only_cross = [x for x in cross_match if len(x) > 0]
# PS1_ids = [[x for x in range(len(PS1_objid)) if PS1_objid[x] == only_cross[y][0]][0] for y in range(len(only_cross))]
# band_diff = [abs(PS1_G_band[PS1_ids[x]] - GLADE_B_band[GLADE_ids[x]]) for x in range(len(GLADE_ids)) if (type(PS1_G_band[PS1_ids[x]]) == float and type(GLADE_B_band[GLADE_ids[x]]) == float)]
# good_b_bands = [GLADE_B_band[GLADE_ids[x]] for x in range(len(GLADE_ids)) if (type(PS1_G_band[PS1_ids[x]]) == float and type(GLADE_B_band[GLADE_ids[x]]) == float)]
#
#
# plt.figure(6)
# print(len([x for x in band_diff if x < 100])/len(band_diff))
# print(np.mean([x for x in good_b_bands if type(x) == float]))
# plt.hist([x for x in band_diff if x < 100], bins=20)
# plt.title("Histogram of B-band and G-band difference of GLADE and PS1 Cross Match")
# plt.xlabel("Band Difference (mags)")
# plt.ylabel("Frequency of difference")
# plt.savefig("images/Crossmatch Band diff.png", bbox_inches="tight", dpi=300)
#
# # for i in range(5):
# #     print(GLADE_ra[GLADE_ids[i]],GLADE_dec[GLADE_ids[i]])
# #     # ps1_index = np.where(PS1_objid == only_cross[i][0])[0]
# #     # ps1_index = [x for x in range(len(PS1_objid)) if PS1_objid[x] == only_cross[i][0]][0]
# #     print(PS1_ra[PS1_ids[i]], PS1_dec[PS1_ids[i]])
# #     print()
#endregion


PS1_columns = "ObjID", "uniquePspsOBid", "raStack", "decStack", "raMean", "decMean", "ra", "dec", "ng", "gMeanPSFMag", "gMeanPSFMagErr", "gMeanKronMag", "gMeanKronMagErr", "gMeanApMag", "gMeanApMagErr", "nr", "rMeanPSFMag", "rMeanPSFMagErr", "rMeanKronMag", "rMeanKronMagErr", "rMeanApMag", "rMeanApMagErr", "ni", "iMeanPSFMag", "iMeanPSFMagErr", "iMeanKronMag", "iMeanKronMagErr", "iMeanApMag", "iMeanApMagErr", "nz", "zMeanPSFMag", "zMeanPSFMagErr", "zMeanKronMag", "zMeanKronMagErr", "zMeanApMag", "zMeanApMagErr", "ny", "yMeanPSFMag", "yMeanPSFMagErr", "yMeanKronMag", "yMeanKronMagErr", "yMeanApMag", "yMeanApMagErr", "gQfPerfect", "rQfPerfect", "iQfPerfect", "zQfPerfect", "yQfPerfect", "qualityFlag", "objInfoFlag", "primaryDetection", "bestDetection", "class", "prob_Galaxy", "prob_Star", "prob_QSO", "z_phot", "z_photErr", "z_phot0", "extrapolation_Photoz", "ps_score"
GLADE_columns = "id", "Galaxy_id", "Distance_id", "PGC", "Name_GWGC", "Name_HyperLEDA", "Name_2MASS", "Name_SDSS_DR12", "RA", "_Dec", "Coord", "dist", "dist_err", "z_dist", "z_dist_err", "z", "B", "B_err", "B_abs", "J", "J_err", "H", "H_err", "K", "K_err", "flag1", "flag2", "flag3"
# for i in range(len(PS1_columns)):
#     print(i,PS1_columns[i])
for i in range(len(GLADE_columns)):
    print(i,GLADE_columns[i])

# print(np.array([1,23,3])["1st"])










