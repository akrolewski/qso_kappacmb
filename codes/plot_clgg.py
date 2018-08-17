import numpy as np
import healpy as hp
from compute_clgg import *
from parameters import *
import matplotlib.pyplot as plt

ls, Clgg, fsky, shotnoise = Clgg(histogram_file,False,nside,1500)

lcenterbin, binnedgg, sigmavecth = BinClgg(ls, Clgg, fsky, 30,1200, 8)

plt.ion()
plt.figure()
plt.errorbar(lcenterbin,binnedgg-shotnoise,yerr=sigmavecth,fmt = 'o')