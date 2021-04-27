import matplotlib.pyplot as plt


import numpy as np
import healpy as hp

##### Stuff from tutorial, but kinda usless because it won't show up in pycharm
# # The resolution of the map is defined by the NSIDE parameter, which is generally a power of 2.
# NSIDE = 32
# print("Approximate resolution at NSIDE {} is {:.2} deg".format(NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60))
#
# # The function healpy.pixelfunc.nside2npix gives the number of pixels NPIX of the map:
# NPIX = hp.nside2npix(NSIDE)
# print("NPIX = " + str(NPIX))
#
# # The same pixels in the map can be ordered in 2 ways, either RING, where they are numbered in the array in horizontal rings starting from the North pole:
# m = np.arange(NPIX)
# hp.mollview(m, title="Mollview image RING")
# hp.graticule()
#####

# Map won't show because not in jupyter, seek help
print("Reading GW190425 Map")
wmap_map_I = hp.read_map("GW190425_PublicationSamples_flattened.fits.gz")
# By default, read_map loads the first column, for reading other columns you can specify the field keyword.

# write_map writes a map to disk in FITS format, if the input map is a list of 3 maps, they are written to a single file as I,Q,U polarization components:
# Aparently, overwrite is no longer a parameter of write_map (??)
#hp.write_map("my_map.fits", wmap_map_I)
# AttributeError: module 'astropy.io.fits' has no attribute 'new_table'
# I apparently need to update astropy, should ask how to do that maybe on docker

#Visualization
plt.figure(1, figsize=(10,8))
hp.mollview(
    wmap_map_I,
    coord=["G", "E"],
    title="Histogram equalized Ecliptic",
    unit="mK",
    norm="hist",
    min=-1,
    max=1,
)
hp.graticule()
plt.savefig("pedros_figure_1.png")