import matplotlib.pyplot as plt
import numpy as np
import csv
from datetime import datetime
from database_methods import *

# if False:
#     c = 3*10**8
#     x = np.linspace(0,c*0.9,10000)
#     y = np.sqrt((c+x)/(c-x)) - 1
#     cap = np.linspace(0,1, 10000)
#
#     plt.figure(1)
#     plt.scatter(x,y, c = cap, cmap = "viridis",vmin=0, vmax=1)
#     plt.xlabel("Recessional Velocity [m/s]")
#     plt.ylabel("Redshift")
#
#     cbar = plt.colorbar()
#     cbar.set_label("hello")
#     plt.savefig("Redshift vs Velocity.png")
#
# # x = np.zeros(10)
# # index = 0
# # for i in range(6):
# #     x[index] = 1
# #     index = index + 1
# # print(x)
# # print(x[:index])
# with open("local_data/DES_allBands.csv", mode='r') as csv_file:
#     galaxy_code = csv.DictReader(csv_file)
#     index = 0
#     i = 0
#     perc = 5
#     perc_now = perc
#     total_file = 529270
#     print("Loading Galaxy")
#     now = datetime.now()
#     for row in galaxy_code:
#         index = index + 1
#         if i/total_file >= perc_now/100:
#             print(perc_now, "%", datetime.now() - now)
#             perc_now = perc_now + perc
#         i = i + 1
# print(index)
# print(i)
str = "objID	uniquePspsOBid	raStack	decStack	raMean	decMean	ra	dec	ng	gMeanPSFMag	gMeanPSFMagErr	gMeanKronMag	gMeanKronMagErr	gMeanApMag	gMeanApMagErr	nr	rMeanPSFMag	rMeanPSFMagErr	rMeanKronMag	rMeanKronMagErr	rMeanApMag	rMeanApMagErr	ni	iMeanPSFMag	iMeanPSFMagErr	iMeanKronMag	iMeanKronMagErr	iMeanApMag	iMeanApMagErr	nz	zMeanPSFMag	zMeanPSFMagErr	zMeanKronMag	zMeanKronMagErr	zMeanApMag	zMeanApMagErr	ny	yMeanPSFMag	yMeanPSFMagErr	yMeanKronMag	yMeanKronMagErr	yMeanApMag	yMeanApMagErr	gQfPerfect	rQfPerfect	iQfPerfect	zQfPerfect	yQfPerfect	qualityFlag	objInfoFlag	gpetRadius	rpetRadius	ipetRadius	zpetRadius	ypetRadius	gpetR50	rpetR50	ipetR50	zpetR50	ypetR50	primaryDetection	bestDetection	gKronFlux	gKronFluxErr	rKronFlux	rKronFluxErr	iKronFlux	iKronFluxErr	zKronFlux	zKronFluxErr	yKronFlux	yKronFluxErr	class	prob_Galaxy	prob_Star	prob_QSO	z_phot	z_photErr	z_phot0	extrapolation_Photoz".replace("	",", ")
print(str)
