#!/usr/bin/env python

'''
This script uses the oval_tools package to explore the impact of integration
noise and the magnitude of oval shift on AMAREX's ability to correlate the
northern and southern auroral patterns.
'''


from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from oval_tools import Aurora

# Set some constants:
dip_tilt   = 9.6
bkgd_noise = 10 # Rayleighs

# Create items over which to iterate:
AE   = glob('data/AE*.save')       # Activity level
time = np.arange(30, 320, 20)      # Image integration time
offs = np.arange(0, 15/2+.5, .5)   # offset between north/south oval


# Loop over combinations.  Get CCs.
for filename in AE[:1]:
    # Create an array to store correlation Coefficients:
    cc_short = np.zeros( (time.size, offs.size) )
    cc_long  = np.zeros( (time.size, offs.size) )

    # Loop over oval offsets:
    for j, o in enumerate(offs): 
        # Open file:
        data = Aurora(filename)

        # Add tilt:
        data.mutate('n', roll=+o)
        data.mutate('s', roll=-o)

        # Add dayglow:
        data.add_dayglow(pitch=dip_tilt)

        for i, t in enumerate(time):
            # add noise based on time integration:
            data.add_bright_noise('n', bkgd=bkgd_noise, t=t)
            data.add_bright_noise('s', bkgd=bkgd_noise, t=t)

            # Save correlation coeffs:
            cc_short[i,j] = data.corr_hemi('ishort').max()
            cc_long[ i,j] = data.corr_hemi('ilong').max()
            

