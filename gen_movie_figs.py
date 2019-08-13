#!/usr/bin/env python

'''
This script generates many images for creating simulated flyby images
from AMAREX.  It reads average oval characteristics from the AE files,
rotates into geographic coordinates, adds dayglow, and saves an image to file.
'''

import os
import re
from glob import glob
import datetime as dt

import numpy as np
from scipy.interpolate import griddata
from spacepy import coordinates as coord
from spacepy.time import Ticktock

from oval_tools import Aurora

###  IMPORTANT PARAMTERS  ###
dLat = 1#.25
dLon = 1#.25
ival = 'ilong'

# Create output directory
outdir = 'movie_frames/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

# Create time array:
day  = dt.datetime(2018,12,21,0,0,0)
time = [day + dt.timedelta(minutes=int(x)) for x in np.arange(0,24*60, 15)]

# Convert to ticktocks:
ticks = Ticktock(time)

# Create lat-lon grid in geo coordinates:
lat = np.arange( -90, 90 +dLat, dLat)
lon = np.arange(-180, 180+dLon, dLon)
nlon, nlat = lon.size, lat.size
lat,  lon  = np.meshgrid( lat, lon )
xyz_geo = np.zeros( [lon.size*lat.size, 3] )
#xy_geo = np.array([lat.ravel(), lon.ravel()]).transpose()

# Loop through AE files:
for f in glob('data/AE*'):
    # Extract AE range from file:
    match = re.search('\d{3}\_\d{3}',f)
    ae = f[match.start():match.end()]
    aedir = outdir+'/AE_{}/'.format(ae)
    
    # Create folder for results:
    if not os.path.exists(aedir):
        os.mkdir(aedir)
    
    # Open oval data:
    data =  Aurora(f)

    # Create vectors of points:
    npts = data['north']['ilong'].size * 2
    half = int(npts/2)
    mlat = np.zeros(npts)
    mlon = np.zeros(npts)
    I    = np.zeros(npts)

    # Populate 'em
    rad = 1.0 + np.zeros( npts )
    mlon_north, mlat_north = np.meshgrid(data['mlt']*15, data['lat'])
    mlon_south, mlat_south = np.meshgrid(data['mlt']*15, data['lat']*-1)
    mlon[:half], mlon[half:] = mlon_north.ravel(), mlon_south.ravel()
    mlat[:half], mlat[half:] = mlat_north.ravel(), mlat_south.ravel()
    I[   :half], I[   half:] = data['north'][ival].ravel(), \
                               data['south'][ival].ravel()

    # Create an array that can be used by Coords:
    #xyz = np.array( [ rad, mlat, mlon ] ).transpose()
    xyz = np.zeros( npts, 3)
    xyz[:,0] = 
    
    # Loop through time:
    for t in time:
        # Create coord object, add time:
        xyz_now = coord.Coords( xyz, 'SM', 'sph')
        xyz_now.ticks = Ticktock(npts*[t])

        # New coordinates!
        xyz_geo = xyz_now.convert('GEO', 'sph').data[:,1:]

        # Now, interpolate onto old grid:
        Igeo = griddata( xyz_geo, I, xy_geo, method='linear').reshape(lat.shape)

        break
    break

