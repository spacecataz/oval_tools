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
time = np.arange(30, 200, 20)      # Image integration time
offs = np.arange(0, 15/2+.5, .5)   # offset between north/south oval

# Create place to put results:
outdir = 'oval_analysis/'
if not os.path.exists(outdir):
    os.mkdir(outdir)


def plot_cc(time, offsets, cc_short, cc_long, AE_range):
    '''
    Make a plot of CC vs. integration time and oval offset.
    '''

    from matplotlib.ticker import MultipleLocator, FuncFormatter
    
    fig = plt.figure( figsize=(7.3,4) )
    fig.subplots_adjust(left=.12,bottom=.16,right=.84,top=.89,wspace=.07)
    a1, a2 = fig.subplots(1,2)

    levs = np.linspace(0, 1, 25)
    
    cont = a1.contourf(2*offs, time, cc_short, levels=levs)
    cont = a2.contourf(2*offs, time, cc_long,  levels=levs)

    ae = AE_range.split('.')[0].split('_')[1:3]
    fig.suptitle('Interhemis. Correlation: AE={}-{}$nT$'.format(ae[0], ae[1]))
    a1.set_title('LBH Short', size=12)
    a2.set_title('LBH Long',  size=12)

    box = a2.get_position()
    a3 = fig.add_axes( [box.x1+.01, box.y0, .02, box.height] )
    cb=fig.colorbar(cont, cax=a3)
    cb.set_label('2D Corr. Coeff.')
    #cb.ax.yaxis.set_major_locator(MultipleLocator(.2))
    
    for a in (a1,a2):
        a.xaxis.set_major_locator(MultipleLocator(5))
        a.yaxis.set_major_locator(MultipleLocator(60))
        a.xaxis.set_major_formatter(FuncFormatter(
            lambda x,y: '{:.1f}$^{{\\circ}}$'.format(x)))
        a.set_xlabel('Oval Separation')
        a.set_ylabel('Integration Time ($s$)')
    a2.set_ylabel('')
    a2.set_yticklabels([])

    return fig
    
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

        for i, t in enumerate(time):
            # Add dayglow:
            data.add_dayglow(pitch=dip_tilt)

            # add noise based on time integration:
            data.add_bright_noise('n', bkgd=bkgd_noise, t=t)
            data.add_bright_noise('s', bkgd=bkgd_noise, t=t)

            # remove dayglow, leaving dayglow noise:
            data.remove_dayglow()
            
            # Save correlation coeffs:
            cc_short[i,j] = data.corr_hemi('ishort', rect=True).max()
            cc_long[ i,j] = data.corr_hemi('ilong',  rect=True).max()
            
    fig = plot_cc(time, offs, cc_short, cc_long, filename)
    fig.savefig(outdir+filename.split('/')[-1][:-5]+'.png')
