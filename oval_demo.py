#!/usr/bin/env python

'''
A demonstration of the oval_tools module.
'''

import numpy as np
from oval_tools import Aurora
import matplotlib.pyplot as plt

# Create figure, add super title, nudge positions.
fig  = plt.figure( figsize=(8.5,11) )
#fig.suptitle('Oval Analysis Demo', size=20)
fig.subplots_adjust(left=.04,right=.99,wspace=.11,hspace=.04,top=.95,bottom=.07)

# Some kwargs for brevity:
kwargs = {'zlim':[0,1000], 'target':fig}

### RAW OVALS + FEATURE
data = Aurora('data/AE_050_100.save')
data.add_feature('data/Aurora_PolarCap/polararc_eflux3_time20.save',
                 expand=[0,-.5])
out1 = data.add_dial_plot('ishort', hemi='n',loc=421, title='North Pole',
                          **kwargs)
out2 = data.add_dial_plot('ishort', hemi='s',loc=422, title='South Pole',
                          **kwargs)

# For the first time we add text, get center of ovals:
x1, x2 = out1[1].get_position().x1, out2[1].get_position().x0
text_x = x1+(x2-x1)/2.
fig.text(text_x, out1[1].get_position().y1, 'Step 1: Load Ovals + PC Arc',
         ha='center', size=16)

### MUTATE
data.mutate('n', expand= 3, yaw= 5, roll= 5)
data.mutate('s', expand=-3, yaw=-5, roll=-5)
out1 = data.add_dial_plot('ishort', hemi='n', loc=423, **kwargs)
out2 = data.add_dial_plot('ishort', hemi='s', loc=424, **kwargs)
fig.text(text_x, out1[1].get_position().y1,
         'Step 2: Transform Oval Position', ha='center', size=16)


### ADD DAYGLOW & NOISE
data.add_dayglow(pitch=9.6)
data.add_bright_noise('n', bkgd=10)
data.add_bright_noise('s', bkgd=10)
out1 = data.add_dial_plot('ishort', hemi='n', loc=425, **kwargs)
out2 = data.add_dial_plot('ishort', hemi='s', loc=426, **kwargs)
fig.text(text_x, out1[1].get_position().y1,
         'Step 3: Add Dayglow & Noise', ha='center', size=16)

### SUBTRACT DAYGLOW
data.remove_dayglow()
out1 = data.add_dial_plot('ishort', hemi='n', loc=427, **kwargs)
out2 = data.add_dial_plot('ishort', hemi='s', loc=428, **kwargs)
fig.text(text_x, out1[1].get_position().y1,
         'Step 4: Remove Dayglow', ha='center', size=16)

### Add a color bar at the bottom.  Use the existing axes' position to guide.
box1 = out1[1].get_position()
box2 = out2[1].get_position()
width = 1.2*(box2.x0 - box1.x1)
x0 = box1.x1 -.1*width

ax = fig.add_axes( [x0, .06, width, .02] )
cb = fig.colorbar(out2[2], cax=ax, orientation='horizontal')
cb.set_label('LBH Short ($R$)', size=12)
plt.show()
