A Brief Tour of the Dark Sky Simulations Early Data Release
===========================================================

This material is intended to be useful to a wide audience, so it
starts at a very basic level.  It ends at a very advanced level. Skip
ahead if you need to, or slow down and review the supplemental
material.

Get the Release Metadata
------------------------
First, you need to get the package with the basic data release.  You can use your web
browser download feature to get a .zip or .tar.gz file.  They can be found at
<https://bitbucket.org/darkskysims/data_release/downloads> Or, if you prefer the command line,

    wget https://bitbucket.org/darkskysims/data_release/get/default.tar.gz

The file you retrieve will be named something like darkskysims-data\_release-dba678c57371.tar.gz
and is currently about 4 Mbytes in size. You should unpack it.

    tar xzf darkskysims-data_release-dba678c57371.tar.gz
    cd darkskysims-data_release-dba678c57371/ds14_a

If you are familiar with the mercurial version control system, and it is installed on
your system, you can replace the steps above with

    hg clone https://bitbucket.org/darkskysims/data_release
    cd data_release/ds14_a

Using mercurial has additional advantages, since it can tell you if any of the files
have been changed or corrupted.  Read more about mercurial at
<http://mercurial.selenic.com/>

Data Exploration
----------------
Let's look at some metadata,

    head ds14_a_1.0000.head

should produce the following output on a Unix-like operating system If you don't have
UNIX, use whatever tool is available to view an ASCII file.

    # SDF 1.0
    int header_len =  2528;
    parameter byteorder = 0x78563412;
    int version = 2;
    int version_2HOT = 2;
    int units_2HOT = 2;
    int64_t npart = 1073741824000;
    float particle_mass = 5.6749434;
    int iter = 543;
    int do_periodic = 1;

This is an SDF header.  Our raw data is distributed in the SDF format.  The data\_release
directory just contains the initial part (.head) of this large file.  Otherwise, you
would need to wait to download the 34 Terabytes in the ds14\_a\_1.0000 file You can read
more about SDF at <http://bitbucket.org/JohnSalmon/sdf>.  If you browse further down in the
file (perhaps using the more or less command) you will see,

    double Omega0_m = 0.295037918703847;
    double Omega0_lambda = 0.7048737821671822;
    double Omega0_b = 0.04676431995034128;
    double h_100 = 0.6880620000000001;
    char length_unit[] = "kpc";
    char mass_unit[] = "1e10 Msun";
    char time_unit[] = "Gyr";
    char velocity_unit[] = "kpc/Gyr";
    char compiled_version_nln[] = "2HOT_nln-1.1.0-17-gbb2d669";
    char compiled_date_nln[] = "Apr 19 2014";
    char compiled_time_nln[] = "06:34:59";

That is the metadata in the SDF file which describes the cosmological parameters used,
the physical units, and information describing the exact version of the code.  A bit
further down is the first structure declaration,

    struct {
        unsigned int sha1_len;
        unsigned char sha1[20];
    }[65536];

This structure describes the layout in the data file of the sha1 checksums.  There are
65536 of them, each with the length of data that was checksummed and the sha1 checksum
value.  This allows one to verify the integrity of each segment of the 34 Terabyte file
independently.  Since the sha1 values are contained in the .head file in the
data\_release, their integrity can be verified with the checksum contained within the
mercurial repository.  A change to any one of the 272 trillion bits in the data file can
thereby be detected.

The main purpose of the SDF header is found in the last structure definition,

    struct {
        float x, y, z;              /* position of body */
        float vx, vy, vz;           /* velocity of body */
        int64_t ident;              /* unique identifier */
    }[1073741824000];

This tells us how to read the positions, velocities and identity of the particles.
There are 1073741824000 of them.  Where are they?  We tell you in the .url file.

    cat ds14_a_1.0000.url
	http://darksky.slac.stanford.edu/simulations/ds14_a/ds14_a_1.0000

The full 34 TB data file is on a machine at Stanford University (feel free to explore
the server at http://darksky.slac.stanford.edu/simulations/). You could download it
like any of the millions of other files on the Internet, but we have better ways.

SDF-enabled Exploration
-----------------------

If you like C and the command line, the SDF library is a good choice.  Check out the
development version from our repository and compile it.  It will help for the next step
if you have libcurl installed first.

    cd ..
    git clone http://bitbucket.org/darkskysims/sdf.git SDF
    cd SDF
    make

This will build an executable called SDFcvt.  It can be used to browse SDF files.
cp SDFcvt /usr/local/bin or elsewhere in your path if you like.  Change back to the
ds14\_a subdirectory of the data\_release, and try it out.

    cd ../ds14_a
    ../SDF/SDFcvt ds14_a_1.0000.head Omega0_m Omega0_lambda Omega0_b h_100
	0.29503791870384699 0.70487378216718222 0.046764319950341277 0.68806200000000006

SDFcvt has parsed the cosmological values you specified from the header.  Now try this,

    ../SDF/SDFcvt -n 1 -s 1000000000000 http://dsdata.org/ds14_a_1.0000 x y z vx vy vz
	3222417 3812104 2801807.5 -216.09906 -215.322311 273.964539

That is the position and velocity of the trillionth particle in the data file.
-n 1 specifies that you want to read 1 element from the structure
-s specifies the offset in the array
(Note that dsdata.org just redirects to the Stanford web site listed above.  Its purpose
is to make URLs shorter, so DarkSky names are easy to type, or fit on a line.)

If you see a message like this, you do not have libcurl installed

    SFhdrio: MPMY_Fread returns -1, errno=22

If you don't have libcurl, read on to see how you can use our python library to access
the data files over the Internet.  You can also download any of the smaller data files
and use SDFcvt locally.

Python-based Exploration
------------------------

If you like Python, we have implemented several methods to interact with the
data, both directly using Numpy arrays and through yt.  Let's start with
relatively simple methods of interacting with the data.  To start, you'll need
to install a few packages. The simplest method is to use pip to install some 
of the basic packages.  To install yt, we suggest you follow their
documentation (http://yt-project.org/docs/dev-3.0/installing.html). 

First let's get ThingKing:

    pip install thingking

or install from souce, located at http://bitbucket.org/zeropy/thingking.  
ThingKing exposes data on the WWW to Python in a memory-mapped interface.

If you start up a Python interpreter, try the following:

    import thingking
    ds14_a = thingking.HTTPArray("http://darksky.slac.stanford.edu/simulations/ds14_a/ds14_a_1.0000")

``ds14_a`` is now array-like, and you can do things like examine its size:

    print ds14_a.size

    34359739943392

That's an array with 34 trillion elements! Let's look at the first few.
Since we did not prescribe a type to the HTTPArray, it defaults to being
an array of characters. That means that if you try to print, say, the 
first 10 elements, you'll see each character:

    print ds14_a[:10]

    [('#',) (' ',) ('S',) ('D',) ('F',) (' ',) ('1',) ('.',) ('0',) ('\n',)]

That is not so useful, but if we print the data attribute that hangs off ``ds14_a``,
we get something more readable: 

    print ds14_a[:10].data

    # SDF 1.0

Congratulations, you can now examine any part of a 34 TB file that you'd like.
If you were inclined, you could use the information in the header, like the ``header_len`` value
and the sizes of the variables to access the particle data. However, we've
implemented all of that within yt, and suggest we move on to there to access
the particle data in Python.


DarkSky \& yt
-------------

We are working towards integrating our extensions to yt into the main yt
development repository.  Until then, we are maintaining a separate fork, which
you can get by pulling changes from the repository at http://bitbucket.org/darkskysims/yt-dark
into your local yt repository. You may also wish to just download a separate clone of this
repository.

    hg clone http://bitbucket.org/darkskysims/yt-dark
    cd yt-dark
    python setup.py install  # or python setup.py develop 

At this point, it would also be most useful to download the example scripts at 
http://bitbucket.org/darkskysims/darksky\_tour. Like most things, you can do this through
downloading a tar file, or use mercurial:
 
    hg clone http://bitbucket.org/darkskysims/darksky_tour
    cd darksky_tour 

There are a bunch of examples in here. Let's just walk through the first one,
which will create a nice visualization of all the particles within a 50 Mpc box
centered around the most massive galaxy cluster in ds14\_a.  The following is
the ``splat_viz.py`` example from the darksky\_tour:

    import yt
    import numpy as np
    from enhance import enhance
    from yt.utilities.lib.image_utilities import add_rgba_points_to_image
    from darksky_catalog import darksky
    
    # Define a bounding box of 100 Mpc on a side.
    center = np.array([-2505805.31114929,  -3517306.7572399, -1639170.70554688])
    width = 50.0e3 # 5 Mpc
    bbox = np.array([center-width/2, center+width/2])
    
    ds = darksky['ds14_a'].load(midx=10, bounding_box=bbox)
    
    ad = ds.all_data()
    Npix = 1024
    image = np.zeros([Npix, Npix, 4], dtype='float64')
    
    cbx = yt.visualization.color_maps.mcm.RdBu
    col_field = ad['particle_velocity_z']
    
    # Calculate image coordinates ix and iy based on what your view width is
    ix = (ad['particle_position_x'] - ds.domain_left_edge[0])/ds.domain_width[0]
    iy = (ad['particle_position_y'] - ds.domain_left_edge[1])/ds.domain_width[1]
    
    # Normalize the color field so that it doesn't get maxed out
    col_field = (col_field - col_field.min()) / (col_field.mean() + 4*col_field.std() - col_field.min())
    add_rgba_points_to_image(image, ix.astype('float64'), iy.astype('float64'), cbx(col_field))
      
    # Write out a color-enhanced image
    yt.write_bitmap(enhance(image), 'enhanced.png')
    print 'Splatted %i particles' % ad['particle_position_x'].size

This imports yt, a few extras from the darksky\_tour repository, and the
darksky\_catalog.  After defining a bounding box into the entire dataset, we
load the data using a midx level 10 file.  The yt dataset is returned to the 
``ds`` object, and we then manually splat particles onto a canvas, colored by
their velocity along the line of sight, and output an image.
