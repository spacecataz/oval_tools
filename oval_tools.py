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
labels = {'eflux':'Energy Flux', 'avee':'Avg. Energy',
          'ilong':'LBH Long', 'ishort':'LBH Short'}
units = {'eflux':r'$ergs/cm^{2}/s$', 'avee':'$keV$', 'ilong':'Rayleighs',
         'ishort':'Rayleighs'} 

def cc2d(f, t):
    '''
    2D correlation coefficient for two matrices of the same size and shape.
    Based on "Fast Normalized Cross-Correlation" by J. P. Lewis; adapted
    so that no shifting takes place.  Doesn't get more simple than this.
    '''
    corr = ( (f-f.mean()) * (t-t.mean()) ).sum()
    variance = np.sqrt( ((f-f.mean())**2).sum() * ((t-t.mean())**2).sum() )
    
    return corr/variance

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
    >>> from spacepy.plot import set_target
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

    Noise can be added to any variable via the "add_noise" method.  The noise
    is stored separately as self[var+'_n'], so that the original data can 
    still be accessed.  

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

        from numpy import meshgrid, pi, sin, cos
        from scipy.io import readsav

        # Load data from IDL save file.
        temp = readsav(self.filename)

        ## COPY DATA FROM SAVE FILE ##
        # Note that we kill the half-width cells about the midnight boundary.
        # Save lat/lons:
        self['lat'] = temp['lats']
        self['mlt'] = temp['mlts'][1:-1]

        # Calculate variables related to creating permutations:
        self.dMlt = self['mlt'][2] - self['mlt'][1]
        self.dLat = self['lat'][2] - self['lat'][1]
        self.dLon = self.dMlt*15 # This is useful.
        
        # Put fluxes into correct hemispheres:
        self['north']['avee']  = temp['aveenorth' ][:,1:-1]
        self['south']['avee']  = temp['aveesouth' ][:,1:-1]
        self['north']['eflux'] = temp['efluxnorth'][:,1:-1]
        self['south']['eflux'] = temp['efluxsouth'][:,1:-1]
        
        # Calculate some related variables necessary in polar coordinates:
        # Note that phi is defined as the counter-clockwise angle off of the
        # negative SM Y direction.  This is for plotting purposes.
        self.colat = 90.-self['lat']  #radial 
        self.phi   = pi*15/180*self['mlt'] - pi/2.  #Azimuthal+offset.
        self.phi[self.phi<0] += 2*pi

        # Calculate SM XYZ cartesian coords for transformations:
        self.xyz = np.zeros( [3, self.colat.size, self.phi.size] )
        for i, theta in enumerate(pi/180*self.colat):
            for j, phi in enumerate(self.phi):
                self.xyz[:,i,j] = sin(theta)*sin(phi), \
                                  sin(theta)*cos(phi), \
                                  cos(theta)
        
        # Calculate some variables for interpolating:
        self.phi_grid, self.colat_grid = meshgrid(self.phi, self.colat)
                
        # Calculate some variables for mesh plotting:
        self.phi_mesh   = np.linspace(self['mlt'][0] -self.dMlt/2,
                                      self['mlt'][-1]+self.dMlt/2,
                                      self['mlt'].size+1)
        self.colat_mesh = np.linspace(self['lat'][0] -self.dLat/2,
                                      self['lat'][-1]+self.dLat/2,
                                      self['lat'].size+1)
        self.phi_mesh   = pi*15/180*self.phi_mesh - pi/2.
        self.colat_mesh = 90.-self.colat_mesh

        # Calculate some extra variables (e.g., brightness):
        self._calc_bright()
        
        return True

    def _calc_bright(self):
        '''
        The brightness observed in the GAIA LBH passbands (long and short) 
        for given precipitation parameters.
    
        This is an empirical model determined by fitting results obtained by 
        running GLOW with Franck-Condon factors combined with the nominal GAIA 
        LBH passbands, including the effect of O2 absorption.
    
        INPUTS (handled by **self**):
        eflux - [erg/cm^2/s] The energy flux of precipitating electrons.
        e     - [eV]         The average energy of precipitating electrons.
                           A Maxwellian distribution is assumed.
        OUTPUTS (saved as self[hemi][Ilong], etc.):
        Ilong  - [R]   Brightness in the LBH long channel. 
        Ishort - [R]  Brightness in the LBH short channel. 

        The triangular passband is accounted for in all outputs.
        '''

        # Energy dependence:
        # Fitted 6th order polynomial coefficients
        cL = np.array([  3.21561119,   0.31176226,  -9.58963207, -14.42459333,
                         -8.84892016,  27.75758128,  66.92209348])
        cS = np.array([ -25.53341678,   33.15085169,   93.14300319,
                        -66.30232109, -166.32052898,    2.3044538 ,
                        171.21479873])

        for hemi in ['north', 'south']:
            x = np.log10(1000*self[hemi]['avee']+1E-4) - 3.0
            Ilong = np.polyval(cL, x)
            Ishort = np.polyval(cS, x)
            
            # Brightness scales with flux
            self[hemi]['ilong']  = self[hemi]['eflux'] * Ilong
            self[hemi]['ishort'] = self[hemi]['eflux'] * Ishort
        
    def add_dial_plot(self, var, hemi='north', target=None, loc=111,
                      zlim=None, title=None, add_cbar=False, clabel=None,
                      cmap='hot', dolog=False, lat_ticks=15, show_noise=True,
                      *args, **kwargs):
        '''
        Add a dial plot of variable **var** using hemisphere **hemi**.  
        The **target** and **loc** syntax sets the default location of the
        resulting plot.  The resulting matplotlib objects (figure, axes, etc.)
        are returned to the caller for further customization.

        The plotting method is Matplotlib's *pcolormesh*.  The grid points
        in the original file are used as cell centers for each pixel.

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
            Set the color map to be used.  Defaults to 'hot'.
        lat_ticks : int
            Set the cadence for latitude ticks.  Defaults to 15 degrees.
        show_noise : bool
            Control if image noise is included on plot.  Default is True.

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

        # Get variable to plot.  Include noise if requested.
        if var+'_n' not in self[hemi]:
            self[hemi][var+'_n'] = np.zeros( self[hemi][var].shape )
        z = self[hemi][var] + show_noise*self[hemi][var+'_n']
        
        # Get max/min if none given:
        if zlim is None:
            zlim=[0,0]
            zlim[0]=z.min(); zlim[1]=z.max()
            # If log scale, ensure minimum does not cross zero:
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

        # Set contour levels.  Safest to do this "by hand":
        if dolog:
            # Set levels evenly in log space: NOT NEEDED FOR PCOLORMESH
            #levs = np.power(10, np.linspace(np.log10(zlim[0]), 
            #                                np.log10(zlim[1]), nlev))
            z=np.where(z>zlim[0], z, 1.01*zlim[0])
            norm=LogNorm()
            ticks=LogLocator()
            fmt=LogFormatterMathtext()
        else:
            # LEVS NOT NEEDED FOR PCOLORMESH
            #levs = np.linspace(zlim[0], zlim[1], nlev)
            norm=None
            ticks=None
            fmt=None

        ### Create Plot ###
        # Add contour to axes:
        cont = ax.pcolormesh(self.phi_mesh, self.colat_mesh, z,
                             cmap=cmap, *args, norm=norm, **kwargs)
        ax.grid() #pcolor turns off grid.  :P
        #cont = ax.contourf(self.phi, self.colat, z, levs, cmap=cmap, *args,
        #                   norm=norm, **kwargs)

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
    
    def mutate(self, hemi='north', yaw=0., pitch=0.0, roll=0.0, expand=0.,
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

        1.  The pattern will be shifted in latitude (*expand* kwarg)
        2.  The pattern will be rotated in local time (*yaw* kwarg)
        3.  The pattern will be shifted (*pitch* and *roll* kwargs)
        4.  If #3 is performed, the pattern is interpolated back onto
            the original grid (bi-linear interpolation.)

        Multiple calls to this method with one mutation at a time can 
        achieve different orders if desired.

        Parameters
        ==========
        No required params.

        Other Parameters
        ==========
        hemi : string
            Select which hemisphere to change.  Defaults to 'north'.
            Single letters can be used (e.g., "n" for "north").
        expand : float
            Latitude, in degrees, to grow (or shrink if negative) the oval.
        yaw : float
            Angle (degrees) to rotate the oval picture about the magnetic pole.
        roll : float
            Angle (degrees) to rotate the picture about the SM X axis.
        pitch : float
            Angle (degrees) to roate the picture about the SM Y axis.

        '''

        from numpy import pi, cos, sin, matmul, arccos, arctan2
        from scipy.interpolate import LinearNDInterpolator as LinInt
        from scipy.interpolate import griddata
        
        ### Set hemisphere name ###
        if hemi[0].lower()=='n':
            hemi='north'
        else:
            hemi='south'

        # Convert pitch/roll to radians:
        roll  *= pi/180.
        pitch *= pi/180.
            
        ### Expand oval ###
        lat_shift = int( np.round(expand / self.dLat) )
        self._roll(lat_shift, 0, hemi)

        ### Rotate oval ###
        lon_shift = int( np.round(yaw / self.dLon) )
        self._roll(0, lon_shift, hemi)
        
        ### Shift oval ###
        # Roll: rotation about x-axis:
        rot_roll = np.array( [[1,         0,          0],
                              [0, cos(roll), -sin(roll)],
                              [0, sin(roll),  cos(roll)]] )
        # Pitch: rotation about y-axis:
        rot_pitch= np.array( [[ cos(pitch),  0,  sin(pitch)],
                              [          0,  1,           0],
                              [-sin(pitch),  0,  cos(pitch)]] )
        # Combine rotation:
        rot = matmul( rot_roll, rot_pitch )
        
        # Apply rotations to XYZ matrix:
        xyz_rot = np.zeros( [3,self.colat.size, self.phi.size] )
        for i in range(self.colat.size):
            for j in range(self.phi.size):
                    xyz_rot[:,i,j] = matmul(self.xyz[:,i,j], rot)

        new=np.array( [xyz_rot[0,:,:].ravel(),
                       xyz_rot[1,:,:].ravel()]).transpose()
        old=np.array( [self.xyz[0,:,:].ravel(),
                       self.xyz[1,:,:].ravel()]).transpose()
        shape = self[hemi]['eflux'].shape
        for key in self[hemi].keys():
            if key[-2:] == '_n': continue # do not shift noise.
            f =  LinInt(new, self[hemi][key].ravel(), fill_value=0)
            self[hemi][key] = f(old).reshape(shape)

        ### NOTE:  This old portion of the code attempted to do the
        ### interpolation in colat-phi space.  This produced really weird
        ### results.  It didn't work very well.
        
        ## Convert XYZ to new colat/phi:
        #colat = arccos(xyz_rot[2,:,:])
        #phi   = arctan2(xyz_rot[0,:,:], xyz_rot[1,:,:])
        #phi[phi<0] += 2*pi
        #
        ## Re-bin to original grid:
        #old = np.array( [self.colat_grid.ravel()*np.pi/180,
        #                 self.phi_grid.ravel()]).transpose()
        #pts = np.array( [colat.ravel(), phi.ravel()]).transpose()
        #shape = self[hemi]['eflux'].shape
        #for key in self[hemi].keys():
        #    if key[-2:] == '_n': continue # do not shift noise.
        #    f =  LinInt(pts, self[hemi][key].ravel(), fill_value=0)
        #    self[hemi][key] = f(old).reshape(shape)
        #
        ##return colat*180/np.pi, phi

    def add_white_noise(self, var, SNR=10, hemi='north'):
        '''
        Add guassian white noise to variable *var* that creates a 
        signal-to-noise (SNR) ratio of *SNR*.

        Parameters
        ==========
        var : string
            Name of variable onto which noise will be added.

        Other Parameters
        ==========
        SNR : float
            The desired signal-to-noise ratio.  Defaults to 10.
        hemi : string
            Select which hemisphere to change.  Defaults to 'north'.
            Single letters can be used (e.g., "n" for "north").
        '''

        from numpy.random import normal
        
        ### Set hemisphere name ###
        if hemi[0].lower()=='n':
            hemi='north'
        else:
            hemi='south'
        
        # Using SNR = Power_signal / Power_noise, and Power_noise =
        # noise variance for gaussian white noise, calculate the noise
        # standard deviation required to get the desired SNR.
        power   = np.sqrt(np.mean(self[hemi][var]**2))
        std_dev = np.sqrt(power/SNR)

        # Generate noise:
        noise = normal(0, std_dev, self[hemi][var].shape)

        # Save noise related to this variable:
        self[hemi][var+'_n'] = noise

        return True

    def add_bright_noise(self, hemi='n', t=30, n=256):
        '''
        Return realistic 1-sigma noise level for a brightness I. In this 
        case "realistic" is intended to mean "similar to the performance 
        we will be showing in the instrument section of the proposal." 
        Let's double check everything before we assume that statement is 
        true. This uses worst-case QE which is at the long edge of the long 
        channel. 
    
        Intended usage:
        hemi = hemisphere to act on.
        sigI = sigma_brightness(I)
        I_noisy = I + sigI*np.random.randn()
        
        INPUTS:
        hemi   - string   Which hemisphere to act on.
        I      - [R]      Brightness on which noise will be added
        t      - [s]      Integration time (default 30)
        n      - [pixels] Number of pixels horizontally and vertically 
                            on the detector (which is assumed to view a 
                            15x15 deg field of view). (Default 256)
                        
        OUTPUTS:
        Noise is added to both Ilong and Ishort variables.
    
        '''
              
        ### Set hemisphere name ###
        if hemi[0].lower()=='n':
            hemi='north'
        else:
            hemi='south'
        
        opt_eff = 0.0058
        fov = 15./n *  np.pi/180. # rad
        etendue = 3.6 * fov**2 # cm^2 ster
        count_rate = 1e6/(4*np.pi) * etendue * opt_eff # counts/pixel/sec/R
        sens = count_rate * t

        for I in ('ilong', 'ishort'):
            self[hemi][I+'_n'] = np.sqrt(sens*self[hemi][I])/sens
    
    def corr_hemi(self, var, noise=True):
        '''
        Compute the cross-correlation matrix between the northern and southern
        hemisphere for variable *var*.  The noise filter is included in this 
        process.

        Parameters
        ==========
        var : string
            Name of variable onto which noise will be added.
        '''

        from scipy.signal import correlate2d

        # Make sure noise filters are available:
        if var+'_n' not in self['north']:
            self['north'][var+'_n'] = np.zeros( self['north'][var].shape )
        if var+'_n' not in self['south']:
            self['south'][var+'_n'] = np.zeros( self['south'][var].shape )

        # Get variables.  Add noise to signals.
        # We must subtract means to get normalized (pearson) CC.
        x = (self['north'][var] + self['north'][var+'_n'])
        x-=x.mean()
        y = (self['south'][var] + self['south'][var+'_n'])
        y-=y.mean()

        # Normalize final answer by standard deviations:
        return correlate2d(x,y,mode='same',boundary='wrap')/(x.std()*y.std()*x.size)
