#!/usr/bin/env python

'''
This script generates many images for creating simulated flyby images
from AMAREX.  It reads average oval characteristics from the AE files,
rotates into geographic coordinates, adds dayglow, and saves an image to file.
'''

import os
import datetime as dt

from spacepy import coordinates as coord
from spacepy.time import Ticktock

# Create output array
outdir = 'movie_frames/'
if not os.path.exists(outdir):
    os.mkdir(outdir)


# Create time array:
day  = dt.datetime(2018,12,21,0,0,0)
time = [day + dt.timedelta(minutes=int(x)) for x in np.arange(0,24*60, 15)]


>>> from spacepy import coordinates as coord
>>> cvals = coord.Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
>>> cvals.x  # returns all x coordinates
array([1, 1])
>>> from spacepy.time import Ticktock
>>> cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticks
>>> newcoord = cvals.convert('GSM', 'sph')
>>> newcoord
