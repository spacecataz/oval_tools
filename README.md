# Oval Tools

This is a basic package for performing investigations into auroral oval
characteristics.

## Installation
The package is terribly simple- just add **oval_tools.py** to your Python path.

## Prerequisites
Prereqs are simple: keep your Python, Numpy, Matplotlib, and Scipy up to date.

|Package | Version |
--------------------
|Python | 3.X  |
|Numpy  | 1.16.X |
|Matplotlib | 3.1.X |
|Scipy | 1.3.X |

Older versions may work, but I'm not about to test that today.

## Included Data
The *data* folder includes data required to run the code and derived analysis
scripts.

The *\*.save* files are IDL-generated files of energy flux
and average energy values across the northern hemisphere obtained via
observations.  Each file are these values averaged over AE bins indicated
in the file name
(i.e., *AE_[lower AE bin boundary]_[upper AE bin boundary].save*).
The values inside the files are binned by magnetic latitude and magnetic
local time.