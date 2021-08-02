import numpy as np
import csv
from datetime import datetime
from database_methods import *
import healpy as hp
from ligo.skymap.postprocess import find_greedy_credible_levels
from ligo.skymap import distance

script_start = datetime.now()

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
print("Loading PS1 Galaxy - " + str(start.time()))
db_query = '''
SELECT * FROM PS1_Galaxy_v4;
'''
PS1 = query_db([db_query])[0]
print("Len PS1 = " + str(len(PS1)))
print("Time to Complete: " + str(datetime.now() - start))

### Load in PS1 ps_score CSV
start = datetime.now()
print("Starting Load Point Source Score CSV - " + str(start.time()))
total_file = 23800
num_in_array = int(total_file*1)
array_index = np.linspace(0, total_file, num_in_array, dtype=int)
with open("local_data/PS1_PS_BoxLimitJoined_pjquinonez.csv", mode='r') as csv_file:
    galaxy_code = csv.DictReader(csv_file)
    PS1_ps_objid = np.zeros(num_in_array)
    PS1_ps = np.zeros(num_in_array)
    PS1_ps_ra = np.zeros(num_in_array)
    PS1_ps_dec = np.zeros(num_in_array)
    index = 0
    i = 0
    perc_now = 0
    print("Loading PS1 Point Source Scores")
    for row in galaxy_code:
        if i/total_file >= perc_now/100:
            print(perc_now, "%", datetime.now() - start)
            perc_now = perc_now + 5
        if i in array_index:
            PS1_ps_objid[index] = int(row["objid"])
            PS1_ps[index] = float(row["ps_score"])
            PS1_ps_ra[index] = float(row["raMEAN"])
            PS1_ps_dec[index] = float(row["decMEAN"])
            index = index + 1
        i = i + 1

PS1_ps_objid = PS1_ps_objid[:index]
PS1_ps = PS1_ps[:index]
PS1_ps_ra = PS1_ps_ra[:index]
PS1_ps_dec = PS1_ps_dec[:index]
print("Len  PS1 Point Source Scores =",len(PS1_ps_objid))
print("Time to Complete: " + str(datetime.now() - start))


### Only Galaxies in 90% Intervals
PS1_ra = np.array([x[4] for x in PS1])
PS1_dec = np.array([x[5] for x in PS1])

start = datetime.now()
print("Start Limit to GW Zone - " + str(start.time()))
PS1_bools = np.ones(len(PS1), dtype=bool)
PS_bools = np.ones(len(PS1_ps), dtype=bool)

for i in range(len(PS1_bools)):
    phi = np.deg2rad(PS1_ra[i])
    theta = 0.5*np.pi - np.deg2rad(PS1_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    if credible_levels[this_pix] > 0.90:
        PS1_bools[i] = False
for i in range(len(PS_bools)):
    phi = np.deg2rad(PS1_ps_ra[i])
    theta = 0.5*np.pi - np.deg2rad(PS1_ps_dec[i])
    this_pix = hp.ang2pix(nside, theta, phi)
    if credible_levels[this_pix] > 0.90:
        PS_bools[i] = False

PS1 = [PS1[x] for x in range(len(PS1)) if PS1_bools[x]]
print("New Len PS1: " + str(len(PS1)))

PS1_ps_objid = PS1_ps_objid[PS_bools]
PS1_ps_ra = PS1_ps_ra[PS_bools]
PS1_ps_dec = PS1_ps_dec[PS_bools]
PS1_ps = PS1_ps[PS_bools]
print("New Len Point Source: " + str(len(PS1_ps)))
print("Finished Limit to GW Zone: " + str(datetime.now() - start))



print("Time to Complete Script: " + str(datetime.now() - script_start))