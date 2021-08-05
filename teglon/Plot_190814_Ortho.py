import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.pyplot import cm
from matplotlib.patches import CirclePolygon
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from database_methods import *

import sys

sys.path.append('../')

import os
import optparse

from configparser import RawConfigParser
import multiprocessing as mp
import mysql.connector

import mysql.connector as test

print(test.__version__)

from mysql.connector.constants import ClientFlag
from mysql.connector import Error
import csv
import time
import pickle
from collections import OrderedDict

import numpy as np
from scipy.special import erf
from scipy.optimize import minimize, minimize_scalar
import scipy.stats as st
from scipy.integrate import simps
from scipy.interpolate import interp2d

import astropy as aa
from astropy import cosmology
from astropy.cosmology import WMAP5, WMAP7, LambdaCDM
from astropy.coordinates import Distance
from astropy.coordinates.angles import Angle
from astropy.cosmology import z_at_value
from astropy import units as u
import astropy.coordinates as coord
from dustmaps.config import config
from dustmaps.sfd import SFDQuery

import shapely as ss
from shapely.ops import transform as shapely_transform
from shapely.geometry import Point
from shapely.ops import linemerge, unary_union, polygonize, split
from shapely import geometry
from shapely.geometry import JOIN_STYLE

import healpy as hp
from ligo.skymap import distance

from HEALPix_Helpers import *
from Tile import *
from SQL_Polygon import *
from Pixel_Element import *
from Completeness_Objects import *

import psutil
import shutil
import urllib.request
import requests
from bs4 import BeautifulSoup
from dateutil.parser import parse

import glob
import gc
import json


isDEBUG = False

class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--gw_id', default="S190814bv", type="str",
                          help='LIGO superevent name, e.g. `S190814bv` ')

        parser.add_option('--healpix_dir', default='../Events/{GWID}', type="str",
                          help='Directory for where to look for the healpix file.')

        parser.add_option('--healpix_file', default="GW190814_PublicationSamples_flattened.fits.gz,0",
                          type="str", help='Healpix filename.')

        return (parser)

    def main(self):

        # Load map pixels from file...
        # Note: I am only downscaling the map for performance reasons... this isn't technically part of the plotting
        # code
        nside256 = 256
        nside1024 = 1024

        hpx_path = "{}/{}".format(self.options.healpix_dir, self.options.healpix_file).replace('{GWID}', self.options.gw_id)
        orig_prob, orig_distmu, orig_distsigma, orig_distnorm, header_gen = hp.read_map(hpx_path, field=(0, 1, 2, 3),
                                                                                        h=True)
        orig_npix = len(orig_prob)

        rescaled_prob, rescaled_distmu, rescaled_distsigma, rescaled_distnorm = hp.ud_grade([orig_prob,
                                                                                             orig_distmu,
                                                                                             orig_distsigma,
                                                                                             orig_distnorm],
                                                                                            nside_out=nside256,
                                                                                            order_in="RING",
                                                                                            order_out="RING")

        rescaled_npix = len(rescaled_prob)
        print("\tRescaled number of pix in '%s': %s" % (self.options.healpix_file, rescaled_npix))

        rescaled_map_nside = hp.npix2nside(rescaled_npix)
        print("\tRescaled resolution (nside) of '%s': %s\n" % (self.options.healpix_file, rescaled_map_nside))

        original_pix_per_rescaled_pix = orig_npix / rescaled_npix
        print("Original pix per rescaled pix for %s" % original_pix_per_rescaled_pix)

        print("Renormalizing and initializing rescaled map...")
        rescaled_prob = rescaled_prob * original_pix_per_rescaled_pix

        map_pix = []
        for i, p in enumerate(rescaled_prob):
            map_pix.append(Pixel_Element(i, nside256, p))

        # Original Resolution (SLOW)
        # map_pix = []
        # for i, p in enumerate(orig_prob):
        #     map_pix.append(Pixel_Element(i, nside1024, p))


        map_pix_sorted = sorted(map_pix, key=lambda x: x.prob, reverse=True)
        print("...done")

        cutoff_50th = 0.5
        cutoff_90th = 0.9
        index_50th = 0
        index_90th = 0

        print("Find index for 50th...")
        cum_prob = 0.0
        for i in range(len(map_pix_sorted)):
            cum_prob += map_pix_sorted[i].prob
            index_50th = i

            if (cum_prob >= cutoff_50th):
                break

        print("... %s" % index_50th)

        print("Find index for 90th...")
        cum_prob = 0.0
        for i in range(len(map_pix_sorted)):
            cum_prob += map_pix_sorted[i].prob
            index_90th = i

            if (cum_prob >= cutoff_90th):
                break
        print("... %s" % index_90th)


        print("Build multipolygons...")
        net_50_polygon = []
        for p in map_pix_sorted[0:index_50th]:
            net_50_polygon += p.query_polygon
        joined_50_poly = unary_union(net_50_polygon)

        # Fix any seams
        eps = 0.00001
        merged_50_poly = []
        smoothed_50_poly = joined_50_poly.buffer(eps, 1, join_style=JOIN_STYLE.mitre).buffer(-eps, 1, join_style=JOIN_STYLE.mitre)

        try:
            test_iter = iter(smoothed_50_poly)
            merged_50_poly = smoothed_50_poly
        except TypeError as te:
            merged_50_poly.append(smoothed_50_poly)

        print("Number of sub-polygons in `merged_50_poly`: %s" % len(merged_50_poly))
        sql_50_poly = SQL_Polygon(merged_50_poly)

        net_90_polygon = []
        for p in map_pix_sorted[0:index_90th]:
            net_90_polygon += p.query_polygon
        joined_90_poly = unary_union(net_90_polygon)

        # Fix any seams
        merged_90_poly = []
        smoothed_90_poly = joined_90_poly.buffer(eps, 1, join_style=JOIN_STYLE.mitre).buffer(-eps, 1, join_style=JOIN_STYLE.mitre)

        try:
            test_iter = iter(smoothed_90_poly)
            merged_90_poly = smoothed_90_poly
        except TypeError as te:
            merged_90_poly.append(smoothed_90_poly)

        print("Number of sub-polygons in `merged_90_poly`: %s" % len(merged_90_poly))
        sql_90_poly = SQL_Polygon(merged_90_poly)
        print("... done.")

        fig = plt.figure(figsize=(10, 10), dpi=600)
        ax = fig.add_subplot(111)

        m = Basemap(projection='stere',
                    lon_0=15.0,
                    lat_0=-20.0,
                    llcrnrlat=-35.0,
                    urcrnrlat=-19.5,
                    llcrnrlon=8.0,
                    urcrnrlon=24.5)

        # Scale colormap
        pix_90 = map_pix_sorted[0:index_90th]
        pixel_probs = [p.prob for p in pix_90]
        min_prob = np.min(pixel_probs)
        max_prob = np.max(pixel_probs)
        print("min prob: %s" % min_prob)
        print("max prob: %s" % max_prob)

        norm = colors.Normalize(min_prob, max_prob)

        print("Plotting (%s) `pixels`..." % len(pix_90))
        for i, p in enumerate(pix_90):
            p.plot(m, ax, facecolor=plt.cm.Greys(norm(p.prob)), edgecolor='None', linewidth=0.5,
                   alpha=norm(p.prob) * 0.8)

        print("Plotting SQL Multipolygons")

        sql_50_poly.plot(m, ax, edgecolor='black', linewidth=1.0, facecolor='None')
        sql_90_poly.plot(m, ax, edgecolor='gray', linewidth=0.75, facecolor='None')
        sql_50_poly.plot(m, ax, edgecolor='black', linewidth=2.0, facecolor='None')
        sql_90_poly.plot(m, ax, edgecolor='black', linestyle="-", linewidth=1.5, facecolor='None')

        # draw parallels.
        sm_label_size = 18
        label_size = 32
        title_size = 36

        _90_x1 = 0.77
        _90_y1 = 0.558

        _90_x2 = 0.77
        _90_y2 = 0.40

        _90_text_y = 0.37
        _90_text_x = 0.32

        _50_x1 = 0.60
        _50_y1 = 0.51

        _50_x2 = 0.48
        _50_y2 = 0.40

        _50_text_y = 0.37
        _50_text_x = 0.64

        parallels = np.arange(-90., 90., 10.)
        dec_ticks = m.drawparallels(parallels, labels=[0, 1, 0, 0])
        for i, tick_obj in enumerate(dec_ticks):
            a = coord.Angle(tick_obj, unit=u.deg)

            for text_obj in dec_ticks[tick_obj][1]:
                direction = '+' if a.dms[0] > 0.0 else '-'
                text_obj.set_text(r'${0}{1:0g}^{{\degree}}$'.format(direction, np.abs(a.dms[0])))
                text_obj.set_size(sm_label_size)
                x = text_obj.get_position()[0]

                new_x = x * (1.0 + 0.08)
                text_obj.set_x(new_x)

        # draw meridians
        meridians = np.arange(0., 360., 7.5)
        ra_ticks = m.drawmeridians(meridians, labels=[0, 0, 0, 1])

        RA_label_dict = {
            7.5: r'$00^{\mathrm{h}}30^{\mathrm{m}}$',
            15.0: r'$01^{\mathrm{h}}00^{\mathrm{m}}$',
            22.5: r'$01^{\mathrm{h}}30^{\mathrm{m}}$',

        }

        for i, tick_obj in enumerate(ra_ticks):
            for text_obj in ra_ticks[tick_obj][1]:
                if tick_obj in RA_label_dict:
                    text_obj.set_text(RA_label_dict[tick_obj])
                    text_obj.set_size(sm_label_size)

        sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.Greys)
        sm.set_array([])  # can be an empty list

        tks = np.linspace(min_prob, max_prob, 6)
        tks_strings = []

        for t in tks:
            tks_strings.append('%0.2f' % (t * 100))

        cb = fig.colorbar(sm, ax=ax, ticks=tks, orientation='vertical', fraction=0.04875, pad=0.02,
                          alpha=0.80)  # 0.08951
        cb.ax.set_yticklabels(tks_strings, fontsize=16)
        cb.set_label("2D Pixel Probability", fontsize=label_size, labelpad=9.0)

        cb.ax.tick_params(width=2.0, length=6.0)

        cb.outline.set_linewidth(2.0)

        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(2.0)

        ax.invert_xaxis()

        plt.ylabel(r'$\mathrm{Declination}$', fontsize=label_size, labelpad=36)
        plt.xlabel(r'$\mathrm{Right\;Ascension}$', fontsize=label_size, labelpad=30)

        fig.savefig("{}/Plots/GW190814_Contours.png".format(
            self.options.healpix_dir.replace('{GWID}', self.options.gw_id)), bbox_inches='tight')
        plt.close('all')
        print("... Done.")


if __name__ == "__main__":
    useagestring = """python Plot_190814_Ortho.py [options]

Example with healpix_dir defaulted to 'Events/<gwid>':
python TileCandidatesMap.py --gw_id <gwid> --healpix_file <filename>
"""

    start = time.time()

    teglon = Teglon()
    parser = teglon.add_options(usage=useagestring)
    options, args = parser.parse_args()
    teglon.options = options

    teglon.main()

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `Plot_190814_Ortho` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
