import numpy as np


# Read in .pksc file
def read_pksc(pkfile_name):
	pkfile           = open(pkfile_name,"rb")
	Nhalo            = np.fromfile(pkfile, dtype=np.int32, count=1)[0]
	RTHMAXin         = np.fromfile(pkfile, dtype=np.float32, count=1)
	redshiftbox      = np.fromfile(pkfile, dtype=np.float32, count=1)
	print "\nNhalo = ", Nhalo

	nfloats_perhalo = 10
	npkdata         = nfloats_perhalo*Nhalo
	peakdata        = np.fromfile(pkfile, dtype=np.float32, count=npkdata)
	peakdata        = np.reshape(peakdata,(Nhalo,nfloats_perhalo))

	xpos      = peakdata[:,0] # Units of position: Mpc
	ypos     = peakdata[:,1]
	zpos     = peakdata[:,2]
	vxpos   = peakdata[:,3] # There is something wrong with the velocities: from Martin
	# By the way, I found that the velocities in the DarkSky mock were hard to understand.  
	# Normally to go into redshift space you add the line-of-sight velocity, in units of aH, to the position.  
	# I found when I did this I didn't get the anisotropy in the clustering that I was expecting.  
	#I had to boost the velocities at high redshift.  I tried increasing the velocities by 1/sqrt{a} and by 1/a.  
	#The latter worked better than the former (though neither was as good as I'd hoped).  
	#So I am guessing they have some sort of comoving velocity thing going on -- though what that 
	#quantity means physically I have no idea.  I spent an hour or so flicking through their paper 
	#and web site, but I didn't see anything explaining this to me.

	#So you might want to take a look and see if you can improve upon that ad hoc correction, or at least confirm that you like or dislike the velocities they give you!
	vypos   = peakdata[:,4]
	vzpos   = peakdata[:,5]
	Rth     = peakdata[:,6] # Lagrangian radius = (3*mvir/(4*pi*rho_crit(z)))**(1./3.)
	M500c   = peakdata[:,7] # Msun.  These are halos defined with respect to
	# the radius at which the density is 500 times critical, 200 times critical, and 200
	# times baryon density
	M200c       = peakdata[:,8]
	M200b       = peakdata[:,9]
	return xpos, ypos, zpos, vxpos, vypos, vzpos, Rth, M500c, M200c, M200b