This code computes the angular clustering (Cls) for a mock catalog.

Set the filepaths and nside in parameters.py.

Start by computing the healpixel for each halo in compute_healpixels.py

Then run make_halo_maps_edison.slurm (on edison) to generate the map
(histogram of healpix values) in parallel

Then run plot_clgg.py to compute and plot Cl_{gg}
Note that lmax and the binning can be adjusted in plot_clgg.py