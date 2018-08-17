import numpy as np
import healpy as hp
from read_pksc import read_pksc
from parameters import *

xpos, ypos, zpos, _, _, _, _, _, _, _ = read_pksc(halo_file)
print 'reading done'

ipix = hp.pixelfunc.vec2pix(nside,xpos,ypos,zpos)
ipix.tofile(healpix_file)