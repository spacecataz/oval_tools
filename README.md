# Oval Tools

This is a basic package for performing investigations into auroral oval
characteristics.

## Installation
The package is terribly simple- just add **oval_tools.py** to your Python path.

## Prerequisites
Prereqs are simple: keep your Python, Numpy, Matplotlib, and Scipy up to date.

|Package | Version |
|--------|---------|
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

## Usage
The main subroutine is an APL for handling and manipulating auroral data files.
The main object is the *Aurora* object, which is a subclass of *dict*.  Instantiate
with the following syntax:

```python
>>> import oval_tools
>>> data = oval_tools.Aurora('./data/AE_350_400.save')
```

Basic info about the oval object can be obtained using dictionary-like syntax.  Physical
values, such as energy flux and average energy, are nested *per hemisphere*.  For example,
you can look at the latitude/MLT grid and the southern hemisphere energy flux as follows:

```python
# After the above commands are executed...
>>> print(data['mlt'])
>>> print(data['lat'])
>>> print(data['south']['eflux'])
```

Object methods can be used to investigate, alter, and evaluate the information within the
Aurora objects.  Be sure to use tab-complete and docstrings to learn more:

```python
# After the above commands are executed...
>>> data.add_dial_plot('eflux')
>>> data.mutate('north', rotate=90.)
>>> data.add_dial_plot('eflux')
```