# Initialization
import numpy as np
import matplotlib.pyplot as pl
from scipy import interpolate
from scipy import stats
import healpy as hp
import sys
from scipy.special import legendre

def Clgg(gal_map, gal_mask, N_side, lmax, normalized=False):
	# CMASS map, number counts. Should convert to delta_n/n 
	gal = hp.read_map(gal_map)
	if gal_mask != False:
		gal_mask = hp.read_map(gal_mask)
		gal_mask = hp.ud_grade(gal_mask, N_side)    # Bring to same resolution
	else:
		gal_mask = np.ones_like(gal)

	# The mask is product of PLANCK lensing and CMASS masks
	mask = gal_mask
	fsky = np.sum(mask) * 1. / len(mask)

	# Apply the mask. Same for all, for now
	masked_pl = mask
	masked_gal = gal * mask
	
	# Compute shot noise
	mean_gal = np.sum(masked_gal)/np.sum(mask)
	A_per_healpix = (4*np.pi)/(12*N_side**2)
	shot_noise = (mean_gal/A_per_healpix) ** (-1.)

	# Normalized tells you if the map has already been normalized or is a raw map of counts
	if normalized == False:
		mean_gal = np.sum(masked_gal) / np.sum(mask)
		masked_gal_dn = masked_gal/ mean_gal - 1.
		masked_gal_dn = mask * masked_gal_dn 
	else:
		masked_gal_dn = masked_gal

	Clgg = hp.anafast(masked_gal_dn, lmax = lmax)
	ls = np.arange(len(Clgg))
	pixwinf = hp.pixwin(N_side)[0:len(Clgg)]

	Clgg = Clgg / (pixwinf **2)
	Clgg = Clgg / fsky

	return ls, Clgg, fsky, shot_noise

'''def w_theta(ls, Clgg):
	theta = np.linspace(0,0.1,100) # In radians
	w = np.zeros_like(theta)
	legl = map(lambda ell: legendre(ell))
	all_legendre = legl(ls)
	
	for i in range(len(theta)):
		w[i] = np.sum((2*ls+1)/(4*np.pi)*Clgg*legl'''

def BinClgg(ls, Clgg, fsky, lmin, lmax, Nbins):
	# Bins and plots.  Requires the auto and cross spectra.
	# Requires lmin and lmax and Nbins for the binning
	# lmin may not go lower than 8 (see readme in planck lensing package)
	# Computes error bars using smooth model spectra
	# For Clkk, uses the Planck noise + power spectrum, nlkk.dat
	# For Clgg, the options are to smooth Clgg or to leave it unsmoothed
	# For ClKg, the options are 'zero' for no signal (appropriate for null tests), 'smooth' for
	# smoothed Clkg, or 'none' for no smoothing at all	
	if lmin < 8:
		print 'Your lmin is too low. Please give lmin >= 8'
		return None
	print 'Binning...'
	bins = np.round(np.linspace(lmin, lmax, Nbins+1))   # Bin edges
	bins = bins.astype(int)

	lcenterbin = np.zeros(len(bins)-1)
	sigmavecth = np.zeros(len(bins)-1)
	binnedgg = np.zeros(len(bins)-1)
	for k in range(0, len(bins)-1):  
		lmaxvec = np.arange(bins[k], bins[k+1], 1)
		lcenterbin[k] = np.round(0.5 * (bins[k] + bins[k+1]))   # bin center
		wt_gg = 0
		for l in lmaxvec:
			sigmavecth[k] += fsky * (2. * l + 1.) / (2*Clgg[l]**2)
			binnedgg[k] += fsky * (2. * l + 1.) / (2*Clgg[l]**2) * Clgg[l]
			wt_gg += fsky * (2. * l + 1.) / (2*Clgg[l]**2)
		sigmavecth[k] = 1. / sigmavecth[k]
		binnedgg[k] = binnedgg[k] / wt_gg
	sigmavecth = np.sqrt(sigmavecth)
	
	return lcenterbin, binnedgg, sigmavecth