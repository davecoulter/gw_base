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
hpx_path = "Events/S190814bv/GW190814_PublicationSamples_flattened.fits.gz,0"
wmap_map_I = hp.read_map(hpx_path)
orig_prob, orig_distmu, orig_distsigma, orig_distnorm, header_gen = hp.read_map(hpx_path, field=(0, 1, 2, 3), h=True)
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

mask = wmap_map_I.astype(np.bool)
wmap_map_I_masked = hp.ma(wmap_map_I)
wmap_map_I_masked.mask = np.logical_not(mask)

plt.figure(4, figsize=(10,8))
plt.hist(wmap_map_I_masked.compressed(), bins=1000)
plt.savefig("pedros_figure_4.png")

plt.figure(2, figsize=(10,8))
hp.gnomview(wmap_map_I, rot=[0, 0.3], title="GnomView", unit="mK", format="%.2g")
plt.savefig("pedros_figure_2.png")

plt.figure(3, figsize=(10,8))
hp.mollview(wmap_map_I_masked.filled())
plt.savefig("pedros_figure_3.png")
plt.close(3)


