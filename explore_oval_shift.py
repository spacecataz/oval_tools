#!/usr/bin/env python

'''
This script uses the oval_tools package to explore the impact of integration
noise and the magnitude of oval shift on AMAREX's ability to correlate the
northern and southern auroral patterns.
'''

import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from oval_tools import Aurora

# Set some constants:
dip_tilt   = 33
bkgd_noise = 10 # Rayleighs

# Create items over which to iterate:
AE   = glob('data/AE*.save')       # Activity level
time = np.arange(20, 620, 20)      # Image integration time
offs = np.arange(0, 10/2+.5, .5)   # offset between north/south oval

AE.sort()

# Create place to put results:
outdir = 'oval_analysis/'
if not os.path.exists(outdir):
    os.mkdir(outdir)


def plot_cc(snr, offsets, cc_short, cc_long, AE_range, goal=.8):
    '''
    Make a plot of CC vs. SNR and oval offset.
    '''

    from matplotlib.ticker import MultipleLocator, FuncFormatter
    
    fig = plt.figure( figsize=(7.3,4) )
    fig.subplots_adjust(left=.12,bottom=.16,right=.84,top=.92,wspace=.04)
    a1, a2 = fig.subplots(1,2)

    levs = np.linspace(0, 1, 55)

    # Add raw contours:
    cont = a1.contourf(2*offs, snr[0,:], cc_short, levels=levs)
    cont = a2.contourf(2*offs, snr[1,:], cc_long,  levels=levs)

    # Add line at 30s int time
    #a1.hlines(30, 2*offsets.min(),2*offsets.max(), colors='k', linestyle='--')
    #a2.hlines(30, 2*offsets.min(),2*offsets.max(), colors='k', linestyle='--')
    
    # Add cc=goal line:
    kwargs = {'levels':[goal], 'colors':'red', 'linewidths':2.0}
    a1.contour(2*offs, snr[0,:], cc_short, **kwargs) 
    a2.contour(2*offs, snr[1,:], cc_long,  **kwargs) 

    # Add labels:
    a1.text(.92, .08, 'LBH Short', transform=a1.transAxes, size=12,
            bbox={'fc':'lightgray'}, ha='right')
    a2.text(.08, .08, 'LBH Long', transform=a2.transAxes, size=12,
            bbox={'fc':'lightgray'}, ha='left')
    
    ae = AE_range.split('.')[0].split('_')[1:3]
    fig.suptitle('Interhemis. Correlation: AE={}-{}$nT$'.format(ae[0], ae[1]),
                 size=16)

    # Add colorbar:
    box = a2.get_position()
    a3  = fig.add_axes( [box.x1+.01, box.y0, .02, box.height] )
    cb  = fig.colorbar(cont, cax=a3, ticks=np.arange(0, 1.25, .25))
    cb.set_label('2D Corr. Coeff.')
    a3.hlines(goal, 0,1, colors='red')

    fig.text(.5, .015, 'Interhemispheric Oval Offset', ha='center', size=16)
    for a in (a1,a2):
        a.set_ylim([1.5, 3.1])
        a.xaxis.set_major_locator(MultipleLocator(5))
        a.yaxis.set_major_locator(MultipleLocator(.5))
        a.xaxis.set_major_formatter(FuncFormatter(
            lambda x,y: '{:.1f}$^{{\\circ}}$'.format(x)))
        a.set_ylabel('Median SNR')
    a2.set_ylabel('')
    a2.set_yticklabels([])

    # Remove last label on left axes:
    #labels = [tick.get_text() for tick in a1.get_xticklabels()]
    plt.setp(a1.get_xticklabels()[-1], visible=False)
    plt.setp(a1.get_xticklabels()[-1], visible=False)
    
    return fig

def plot_signal(time, data, filename, style='fivethirtyeight'):
    '''
    Examine noise-to-signal ratio for both channels.
    '''

    from matplotlib.ticker import MultipleLocator
    
    # Temporarily use a specific style sheet:
    with plt.style.context(style):
        # Create figure and adjust.
        fig = plt.figure(figsize=[8,6])
        ax  = fig.add_subplot(111)

    # Plot the two channels.
    l1 = ax.plot(time, data[0,:])[0]
    l2 = ax.plot(time, data[1,:])[0]

    # Add labels:
    ax.text(time[-1], data[0,-1], 'LBH Short', va='center', size=16,
            c=l1.get_color())
    ax.text(time[-1], data[1,-1], 'LBH Long',  va='center', size=16,
            c=l2.get_color())

    # Adjust axes:
    ax.set_xlim( [ax.get_xlim()[0], ax.get_xlim()[-1]*1.15] )
    ax.xaxis.set_major_locator(MultipleLocator(60))
    ax.set_xlabel('Integration Time ($s$)')
    ax.set_ylabel('Median SNR')

    title = filename.split('/')[-1][3:-5].replace('_','-')
    ax.set_title('Signal-to-Noise for AE {} $nT$'.format(title), size=18)

    fig.tight_layout()
    
    return fig

def plot_hemicomp(data):
    '''
    Compare northern and southern hemispheres:
    '''

    # Create and adjust figure:
    fig = plt.figure(figsize=[9.34, 4.57])
    fig.subplots_adjust(left=.05,bottom=.11,right=.95,top=.88,wspace=.28)
    
    # Add ovals:
    out1 = data.add_dial_plot('ilong', target=fig, loc=121)
    out2 = data.add_dial_plot('ilong', hemi='south', target=fig, loc=122)

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

    return fig
    
# Loop over combinations.  Get CCs.
for filename in AE:
    # Create an array to store correlation Coefficients:
    cc_short = np.zeros( (time.size, offs.size) )
    cc_long  = np.zeros( (time.size, offs.size) )

    # Create an array to store signal ratios:
    signal = np.zeros( (2,time.size) )
    
    # Loop over oval offsets:
    for j, o in enumerate(offs): 
        # Open file:
        data = Aurora(filename)
        orig = Aurora(filename) # for signal calculations.
        
        # Add tilt:
        data.mutate('n', roll=+o)
        data.mutate('s', roll=-o)
        orig.mutate('n', roll=+o)
        orig.mutate('s', roll=-o)

        # Get weakest signal of original oval:
        Ilong_min = data['north']['ilong'][data['north']['ilong']>0].min()
        Ishrt_min = data['north']['ilong'][data['north']['ilong']>0].min()

        # Save locations of signal>min and signal itself:
        sig_shrt = orig['north']['ishort'][ orig['north']['ishort']>Ishrt_min ]
        sig_long = orig['north']['ilong' ][ orig['north']['ilong'] >Ilong_min ]
        
        for i, t in enumerate(time):
            # Add dayglow:
            data.add_dayglow(pitch=dip_tilt)

            # add noise based on time integration:
            data.add_bright_noise('n', bkgd=bkgd_noise, t=t)
            data.add_bright_noise('s', bkgd=bkgd_noise, t=t)

            # remove dayglow, leaving dayglow noise:
            data.remove_dayglow()

            # Create and save comparison plot:
            fig = plot_hemicomp(data)
            off_str = '{:.1f}'.format(o).replace('.','p')
            fig.savefig(outdir+'ovals_{}_t{:03d}_offs{}.png'.format(
                filename.split('/')[-1][:-5],t,off_str))
            plt.close(fig)
            
            # Save correlation coeffs:
            cc_short[i,j] = data.corr_hemi('ishort', rect=True).max()
            cc_long[ i,j] = data.corr_hemi('ilong',  rect=True).max()

            # Save brightness ratio:
            std_shrt = np.std(data['north']['ishort'] + \
                              np.abs(data['north']['ishort_n']))
            std_long = np.std(data['north']['ilong']  + \
                              np.abs(data['north']['ilong_n']))
            signal[0,i] += np.median( sig_shrt/std_shrt )
            signal[1,i] += np.median( sig_long/std_long )

    signal /= len(offs)
            
    # Create figures relevant to analysis:
    fig = plot_signal(time, signal, filename, 'fivethirtyeight')
    fig.savefig(outdir+'SNR_{}.png'.format(
        filename.split('/')[-1][:-5],off_str))
    plt.close(fig)
            
    fig = plot_cc(signal, offs, cc_short, cc_long, filename)
    fig.savefig(outdir+filename.split('/')[-1][:-5]+'.png')
    fig.savefig(outdir+filename.split('/')[-1][:-5]+'.pdf')
