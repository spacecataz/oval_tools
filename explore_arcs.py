#!/usr/bin/env python

'''
This script uses the oval_tools package to explore the potential ability to 
identify features from one hemisphere to another.

Specifically, we examine the worst-case scenario: finding an arc feature in
the night side hemisphere and trying to match that feature in the southern 
hemisphere.
'''

import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from oval_tools import Aurora

# Set some constants:
dip_tilt   = 33
bkgd_noise = 10 # Rayleighs

arcs = glob('data/Aurora_PolarCap/polararc_eflux3_time*')

def plot_hemicomp(data):
    '''
    Compare northern and southern hemispheres.
    Do feature finding for both hemis.
    Return the maximum correlation coeff for the comparison.
    '''

    # Calc CC, get location of peak:
    c    = data.corr_feature('ilong')
    ij   = np.unravel_index(np.argmax(c), c.shape)
    x, y = data.phi_mesh[ij[0]], data.colat_mesh[ij[1]]
    
    # Create and adjust figure:
    fig = plt.figure(figsize=[9.34, 4.57])
    fig.subplots_adjust(left=.05,bottom=.11,right=.95,top=.88,wspace=.28)
    
    # Add ovals:
    out1 = data.add_dial_plot('ilong', target=fig, loc=121, feature=True)
    out2 = data.add_dial_plot('ilong', target=fig, loc=122, feature=True,
                              hemi='s')

    # Add loc of matched feature:
    out2[1].plot(x,y,'*',ms=10, c='dodgerblue')

    # Sort some variables for convenience:
    a1, a2, cont = out1[1], out2[1], out1[-2]

    # Add some labels
    a1.set_title('Northern Hemi.', size=16)
    a2.set_title('Southern Hemi.', size=16)
    
    # Add a shared colorbar
    box1 = out1[1].get_position()
    box2 = out2[1].get_position()
    width = 1.75*(box2.x0 - box1.x1)
    x0 = .5*(box1.x1+box2.x0-width)#box1.x1 -.25*width
    ax = fig.add_axes( [x0, .17, width, .02] )
    cb = fig.colorbar(out2[2], cax=ax, orientation='horizontal')
    cb.set_label('LBH Long ($R$)', size=14)

    return fig, c.max()
    
