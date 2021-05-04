import matplotlib
matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import healpy as hp
from ligo.skymap import distance
import numpy as np
from database_methods import *



# print("*** DEBUG ***\n")
#
# # test out database methods
# detector_select = '''
# SELECT id, Name FROM Detector;
# '''
# result = query_db([detector_select])[0]
#
# for r in result:
#     print(r)
#
# print("\n*** DEBUG ***\n")



# DC: Tutorial code and text taken from: https://healpy.readthedocs.io/en/latest/tutorial.html

'''
NSIDE and ordering
Maps are simply numpy arrays, where each array element refers to a location in the sky as defined by the Healpix 
pixelization schemes (see the healpix website).

Note: Running the code below in a regular Python session will not display the maps; it’s recommended to use an IPython 
shell or a Jupyter notebook.

The resolution of the map is defined by the NSIDE parameter, which is generally a power of 2.
'''
NSIDE = 32
print(
    "Approximate resolution at NSIDE {} is {:.2} deg".format(
        NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60
    )
)

'''
The function healpy.pixelfunc.nside2npix gives the number of pixels NPIX of the map:
'''
NPIX = hp.nside2npix(NSIDE)
print(NPIX)
m = np.arange(NPIX)

'''
The same pixels in the map can be ordered in 2 ways, either RING, where they are numbered in the array in horizontal 
rings starting from the North pole:
'''

# DC: note that using Docker, you need to serialize your plots so you can't use the default way of creating a model for
# this...
# Also, I am just showing the syntax for how to create subplots with Healpy's built-in visualization functions
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(224)

# Set the current "active" axis
plt.axes(ax1)
hp.mollview(m, title="Mollview image RING w/ Graticule", hold=True)
hp.graticule()

plt.axes(ax2)
hp.mollview(m, title="Mollview image RING no Graiticule", hold=True)
fig.savefig("mollview_test1.png")
plt.close('all')


'''
The standard coordinates are the colatitude 𝜃, 0 at the North Pole, 𝜋/2 at the equator and 𝜋 at the South Pole and the 
longitude 𝜙 between 0 and 2𝜋 eastward, in a Mollview projection, 𝜙=0 is at the center and increases eastward toward the 
left of the map.

We can also use vectors to represent coordinates, for example vec is the normalized vector that points to 𝜃=𝜋/2,𝜙=3/4𝜋:
'''
vec = hp.ang2vec(np.pi / 2, np.pi * 3 / 4)
print(vec)

'''
We can find the indices of all the pixels within 10 degrees of that point and then change the value of the map at those 
indices:
'''
ipix_disc = hp.query_disc(nside=32, vec=vec, radius=np.radians(10))

m = np.arange(NPIX)
m[ipix_disc] = m.max()

fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(111)
plt.axes(ax1)
hp.mollview(m, title="Mollview image RING", hold=True)
fig.savefig("mollview_test2.png")
plt.close('all')

'''
We can retrieve colatitude and longitude of each pixel using pix2ang, in this case we notice that the first 4 pixels 
cover the North Pole with pixel centers just ~1.5 degrees South of the Pole all at the same latitude. The fifth pixel 
is already part of another ring of pixels.
'''
theta, phi = np.degrees(hp.pix2ang(nside=32, ipix=[0, 1, 2, 3, 4]))
print(theta)
print(phi)

'''
The RING ordering is necessary for the Spherical Harmonics transforms, the other option is NESTED ordering which is very 
efficient for map domain operations because scaling up and down maps is achieved just multiplying and rounding pixel 
indices. See below how pixel are ordered in the NESTED scheme, notice the structure of the 12 HEALPix base pixels 
(NSIDE 1):
'''
m = np.arange(NPIX)

fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(111)
plt.axes(ax1)
hp.mollview(m, nest=True, title="Mollview image NESTED", hold=True)
fig.savefig("mollview_test3.png")
plt.close('all')

'''
All healpy routines assume RING ordering, in fact as soon as you read a map with read_map, even if it was stored as 
NESTED, it is transformed to RING. However, you can work in NESTED ordering passing the nest=True argument to most 
healpy routines.
'''

# The rest of the tutorial is on the webpage -- take a look and get comfortable with Healpy! Next up will be MySQL...