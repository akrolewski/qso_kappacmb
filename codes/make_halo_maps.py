import numpy as np
from read_pksc import read_pksc
from mpi4py import MPI
from astropy.io import fits
from parameters import *

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size() # Fiducial size is 100

# Get number of halos
pkfile           = open(halo_file,"rb")
Nhalo            = np.fromfile(pkfile, dtype=np.int32, count=1)[0]
pkfile.close()

# Use this to define size of the parallelized chunks that I compute the histograms over
# Make each chunk as large as the last chunk would be (which has the roundoff remainder)
chunk_size = Nhalo - (size-1)*(Nhalo/size)
# Populate each chunk with -10s so that I can just histogram and get the right answer,
# given that I have restricted bins to being positive
recv_data = -10.*np.ones(chunk_size,dtype=np.float64)

# Initialize variables in all MPI ranks
ipix = None
counts = None
dspls = None

if rank == 0:
	ipix = np.fromfile(healpix_file,dtype=np.int64)
	ipix = ipix.astype(np.float64)
	
	# Create tuple of "counts" telling me how long each chunk is
	counts = (Nhalo/size)*np.ones(size)
	counts[-1] = chunk_size
	counts = tuple(counts)
	
	# Create displacements telling me where to read each chunk
	dspls = tuple(np.arange(size)*Nhalo/size)

# Scatter to each MPI rank	
comm.Scatterv([ipix,counts,dspls,MPI.DOUBLE],recv_data,root=0)

# Run the histogram
hist = np.histogram(recv_data,bins=np.arange(12*nside**2+1.)-0.5)[0]

del ipix
del recv_data

# Sum over all histograms
allhist = np.zeros_like(hist)
comm.Reduce(hist,allhist,root=0,op=MPI.SUM)

if rank == 0:
	# Write the histogram to fits file
	hdr = fits.Header()
	hdr['NSIDE'] = nside

	hdu = fits.BinTableHDU.from_columns(
	[fits.Column(name='hist',format='E',array=allhist)],header=hdr)

	hdr = fits.Header()
	primary_hdu = fits.PrimaryHDU(header=hdr)

	hdul = fits.HDUList([primary_hdu, hdu])
	hdul.writeto(histogram_file,overwrite=True)

