#!/usr/bin/env python
'''
This is a module for manipulating and exploring auroral oval information for
the purpose of investigating observational requirements.
'''

import numpy as np
import matplotlib.pyplot as plt

# Set a good style sheet because we are not savages.
plt.style.use('seaborn-talk')

# Some default labels for plotting:
labels = {'eflux':'Energy Flux', 'avee':'Avg. Energy'}
units = {'eflux':r'$ergs/cm^{2}/s$', 'avee':'$ergs?$'} #DON'T KNOW AVEE UNITS.

def set_target(target, figsize=None, loc=111, polar=False):
    '''
    This is a helper function for plotting that makes building multi-panel
    plots easier.

    Given a *target* on which to plot a figure, determine if that *target*
    is **None** or a matplotlib figure or axes object.  Based on the type
    of *target*, a figure and/or axes will be either located or generated.
    Both the figure and axes objects are returned to the caller for further
    manipulation.  This is used in nearly all *add_plot*-type methods.

    Parameters
    ==========
    target : object
        The object on which plotting will happen.

    Other Parameters
    ================
    figsize : tuple
        A two-item tuple/list giving the dimensions of the figure, in inches.  
        Defaults to Matplotlib defaults.
    loc : integer 
        The subplot triple that specifies the location of the axes object.  
        Defaults to 111.
    polar : bool
        Set the axes object to polar coodinates.  Defaults to **False**.

    Returns
    =======
    fig : object
      A matplotlib figure object on which to plot.

    ax : object
      A matplotlib subplot object on which to plot.

    Examples
    ========
    >>> import matplotlib.pyplot as plt
    >>> from spacepy.pybats import set_target
    >>> fig = plt.figure()
    >>> fig, ax = set_target(target=fig, loc=211)

    '''
    import matplotlib.pyplot as plt
    
    # Is target a figure?  Make a new axes.
    if type(target) == plt.Figure:
        fig = target
        ax  = fig.add_subplot(loc, polar=polar)
    # Is target an axes?  Make no new items.
    elif issubclass(type(target), plt.Axes):
        ax  = target
        fig = ax.figure
    # Is target something else?  Make new everything.
    else:
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot(loc, polar=polar)

    return fig, ax


class Aurora(dict):
    '''
    This class handles an auroral oval in terms of location, strength, and
    other characteristics.

    Upon instantiation, an IDL save file (indicated with argument *filename*)
    will be loaded into the object for use with the object methods.

    Parameters
    ==========
    filename : string
        A string setting the path of the data file to load into the object.

    To-Do:
       --Create options to read different file types on instantiation.
    '''

    
    def __init__(self, filename, *args, **kwargs):
        '''
        Instantiate object, read file, populate object.
        '''

        # Initialize as empty dict
        super(dict, self).__init__(*args, **kwargs)
        
        # Store filename within object:
        self.filename = filename

        # Initialize basic structure of data object:
        self['north'] = {}
        self['south'] = {}
        
        # Load data.
        self._read_IDL()

    def _read_IDL(self):
        '''
        Load an IDL save file into object.  Should only be called on 
        instantiation.
        '''

        from scipy.io import readsav

        # Load data from IDL save file.
        temp = readsav(self.filename)

        ## COPY DATA FROM SAVE FILE ##
        # Save lat/lons:
        self['lat'] = temp['lats']
        self['mlt'] = temp['mlts']

        # Put fluxes into correct hemispheres:
        self['north']['avee']  = temp['aveenorth' ]
        self['south']['avee']  = temp['aveesouth' ]
        self['north']['eflux'] = temp['efluxnorth']
        self['south']['eflux'] = temp['efluxsouth']
        
        # Calculate some related variables necessary for polar plotting:
        self.colat = 90.-self['lat']  #radial 
        self.phi   = np.pi*15/180*self['mlt'] - np.pi/2.  #Azimuthal+offset.

        # Calculate variables related to creating permutations:
        self.dMlt = self['mlt'][2] - self['mlt'][1]
        self.dLat = self['lat'][2] - self['lat'][1]
        self.dLon = self.dMlt*15 # This is useful.

        return True

    def add_dial_plot(self, var, nlev=101, hemi='north', target=None, loc=111,
                      zlim=None, title=None, add_cbar=False, clabel=None,
                      cmap='hot_r', dolog=False, lat_ticks=15,
                      *args, **kwargs):
        '''
        Add a dial plot of variable **var** using hemisphere **hemi**.  
        The **target** and **loc** syntax sets the default location of the
        resulting plot.  The resulting matplotlib objects (figure, axes, etc.)
        are returned to the caller for further customization.

        If kwarg **target** is None (default), a new figure is 
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object, 
        the plot is placed into that axis.

        Extra arguments and keyword arguments are passed to the 
        *matplotlib.axes._subplots.AxesSubplot.contourf* method.

        Parameters
        ==========
        var : string
            The value to be plotted, e.g., 'eflux'.

        Returns
        =======
        fig : matplotlib figure object
        ax  : matplotlib axes object
        cont : matplotlib contour object
        cbar : matplotlib colorbar object

        Other Parameters:
        =================
        hemi : string
            Set hemisphere to plot.  Defaults to "north".  Single letters can
            be used (e.g., "n" instead of "north".)
        target : Figure or Axes
            If None (default), a new figure is generated from scratch.
            If a matplotlib Figure object, a new axis is created
            to fill that figure.
            If a matplotlib Axes object, the plot is placed
            into that axis.
        loc : int
            Use to specify the subplot placement of the axis
            (e.g. loc=212, etc.) Used if target is a Figure or None.
            Default 111 (single plot).
        nlev : int
            Set the number of contour integers.  Defaults to 101.
        zlim : two-element list or array
            Set the color bar range.  Defaults to max/min of data.
        add_cbar : bool
            Set whether to add a color bar or not.  Defaults to **False**.
        dolog : bool
            If **True**, use a log scale to plot *var*.  Defaults to **False**.
        title : string
            Sets the plot title.  Defaults to 'auto', using the variable label.
        clabel : string
            Set label for the color bar, defaults to *var* and associated units.
        cmap : string color map name
            Set the color map to be used.  Defaults to 'hot_r'.
        lat_ticks : int
            Set the cadence for latitude ticks.  Defaults to 15 degrees.

        Examples:
        =========
        # Quickly examine energy flux:
        >>> import matplotlib.pyplot as plt
        >>> import oval_tools
        >>> plt.ion()
        >>> oval = oval_tools.Aurora('data/AE_350_400.save')
        >>> oval.add_dial_plot('eflux', add_cbar=True, zlim=[0,10])
        # Compare hemispheres on a single plot:
        >>> fig=plt.figure( figsize=(9,4) )
        >>> oval.add_dial_plot('eflux',target=fig,loc=121,hemi='n',zlim=[0,10])
        >>> oval.add_dial_plot('eflux',target=fig,loc=122,hemi='s',zlim=[0,10])
        '''
        from numpy import pi
        import matplotlib as mpl
        from matplotlib.colors import (LogNorm, Normalize)
        from matplotlib.ticker import (LogLocator, LogFormatter, 
                                       LogFormatterMathtext, MultipleLocator)


        ### Plot Prep ###
        # Set hemisphere name:
        if hemi[0].lower()=='n':
            hemi='north'
        else:
            hemi='south'
        
        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10,10), loc=loc, polar=True)

        # Get max/min if none given:
        if zlim is None:
            zlim=[0,0]
            zlim[0]=self[hemi][var].min(); zlim[1]=self[hemi][var].max()
            # If log scale, ensure minimum does not cross zero:
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

        # Set contour levels.  Safest to do this "by hand":
        if dolog:
            # Set levels evenly in log space:
            levs = np.power(10, np.linspace(np.log10(zlim[0]), 
                                            np.log10(zlim[1]), nlev))
            z=np.where(self[hemi][var]>zlim[0], self[hemi][var], 1.01*zlim[0])
            norm=LogNorm()
            ticks=LogLocator()
            fmt=LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            z=self[hemi][var]
            norm=None
            ticks=None
            fmt=None

        ### Create Plot ###
        # Add contour to axes:
        cont = ax.contourf(self.phi, self.colat, z, levs, cmap=cmap, *args,
                           norm=norm, **kwargs)

        # Add cbar as required:
        if add_cbar:
            cbar=plt.colorbar(cont, ax=ax, pad=0.05)
            # Build label:
            if clabel==None: 
                clabel="{} ({})".format(labels[var], units[var])
            cbar.set_label(clabel)
        else:
            cbar=None # Need to return something, even if none.

        ### Customize pot ###
        # Set MLT labels:
        lt_labels = ['06',   '',    '18',   '00'  ]
        xticks    = [   0,   pi/2,   pi,     3*pi/2]
        ax.set_xticks(xticks)
        ax.set_xticklabels(lt_labels)

        # Set better radial range:
        ax.yaxis.set_major_locator(MultipleLocator(lat_ticks))
        ax.set_ylim([0,self.colat.max()])

        # Turn off old tick labels.  They suck.
        ax.set_yticklabels('')

        # Get locations of new ticks:
        yticks = np.floor(self['lat']/lat_ticks)*lat_ticks
        yticks = np.unique(yticks[ yticks>=self['lat'].min()])

        # Use text function to manually add pretty ticks.
        opts = {'size':mpl.rcParams['ytick.labelsize'],
                'rotation':-45, 'ha':'center', 'va':'center'}
        for theta in yticks:
            txt = '{:02.0f}'.format(theta)+r'$^{\circ}$'
            ax.text(pi/4., 90.-theta, txt, color='w', weight='heavy', **opts)
            ax.text(pi/4., 90.-theta, txt, color='k', weight='light', **opts)
        
        # Set title:
        if title: ax.set_title(title)
        
        # Return plot objects:
        return fig, ax, cont, cbar

    def _roll(self, nLat, nMlt, hemi='north'):
        '''
        Roll array by **nLat** and **nMlt** in each respective direction.
        This is a helper method for mutate.
        '''
        from numpy import roll
        
        for x in self[hemi].keys():
            self[hemi][x] = roll(self[hemi][x], [nLat, nMlt], [0,1])

        return True
    
    def mutate(self, hemi='north', rotate=0., expand=0.,
               shift_dir=0, shift_lat=0.0):
        '''
        Mutate an oval picture by rotating, expanding, or shifting.
        All variables (energy flux, ave. energy, etc.) will be affected.
        Only one hemisphere will be mutated at a time.

        Keep in mind that permutations are discretized to the resolution of
        the oval: if the azimuthal resolution is ~5 degrees, a rotation of 
        of 1 degree will do nothing.  The permutations are rounded to the
        nearest integer multiple of the grid spacings.

        The order of the permutations is as follows:

        1. The pattern will be rotated about the poll.
        2. The pattern will be grown/shrunk in latitude.
        3. The directional shift will be applied.

        Multiple calls to this method with one mutation at a time can 
        achieve different orders if desired.

        Parameters
        ==========
        No required params.

        Other Parameters
        ==========
        rotate : float
            Angle (degrees) to rotate the oval picture about the magnetic pole.
        expand : float
            Latitude, in degrees, to grow (or shrink if negative) the oval.
        hemi : string
            Select which hemisphere to change.  Defaults to 'north'.
            Single letters can be used (e.g., "n" for "north").
        '''

        ### Set hemisphere name ###
        if hemi[0].lower()=='n':
            hemi='north'
        else:
            hemi='south'
        
        ### Rotate oval ###
        lon_shift = int( np.round(rotate / self.dLon) )
        self._roll(0, lon_shift, hemi)
        
        ### Expand oval ###
        lat_shift = int( np.round(expand / self.dLat) )
        self._roll(lat_shift, 0, hemi)

        
        
    def add_noise(self):
        pass

    def hemi_corr(self):
        pass
