#!/usr/bin/env python

import oval_tools
from importlib import reload
import matplotlib.pyplot as plt

reload(oval_tools)

a = oval_tools.Aurora('./data/AE_100_150.save')
a.add_feature('data/Aurora_PolarCap/polararc_eflux3_time20.save', sens=25.0)
c = a.corr_feature('ilong', debug=True, flip=False)

a.add_dial_plot('ilong', feature=True, add_cbar=True)
plt.show()
