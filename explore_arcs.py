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

