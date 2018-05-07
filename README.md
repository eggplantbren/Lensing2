Lensing
=======

Doing gravitational lensing again!

(c) 2014 -- 2018 Brendon J. Brewer.

LICENCE
=======

The contents of this repository are licensed under the GNU General Public
Licence, version 3. See the LICENCE file for details.

ACKNOWLEDGEMENTS
================

This work is supported by a Marsden Fast-Start grant
from the Royal Society of New Zealand.

DEPENDENCIES
============

* Armadillo (http://arma.sourceforge.net/)
* yaml-cpp (https://github.com/jbeder/yaml-cpp)

You can probably get these from your operating system's package manager.
You'll also need the header files, which are sometimes in a -dev or
-devel package.

You'll also need git, a recent version of g++, and gfortran for the
compiler to work.

There are some associated Python scripts as well, which use Python3
(you can try Python 2, but no guarantees), numpy,
matplotlib, and the DNest4 python package.

USAGE
=====

First, clone the repository recursively:
```
git clone --recursive https://github.com/eggplantbren/Lensing2
```
If you don't already have DNest4 on your system, install its Python package
by following the instructions [here](https://github.com/eggplantbren/DNest4)
first.



Then compile all the Fortran and C++:
```
make
```


