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









