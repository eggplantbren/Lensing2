## Lensing

Doing gravitational lensing again!

(c) 2014 -- 2018 Brendon J. Brewer.

## License

The contents of this repository are licensed under the GNU General Public
Licence, version 3. See the LICENCE file for details.

## Acknowledgements

This work is supported by a Marsden Fast-Start grant
from the Royal Society of New Zealand. If you use this software, please cite
the following paper, which details the methodology:

Brewer, Brendon J., David Huijser, and Geraint F. Lewis. "Trans-dimensional Bayesian inference for gravitational lens substructures." Monthly Notices of the Royal Astronomical Society 455, no. 2 (2015): 1819-1829.

## Dependencies

* Armadillo (http://arma.sourceforge.net/)
* yaml-cpp (https://github.com/jbeder/yaml-cpp)

You can probably get these from your operating system's package manager.
You'll also need the header files, which are sometimes in a -dev or
-devel package.

You'll also need git to obtain the source code, and non-ancient versions of g++ and gfortran for it to compile properly.

There are some associated Python scripts as well, which use Python3
(you can try Python 2, but no guarantees), numpy,
matplotlib, and the DNest4 python package.

## Downloading and compiling

First, clone the repository recursively:
```
git clone --recursive https://github.com/eggplantbren/Lensing2
```

Then compile all the Fortran and C++:
```
make
```

If you don't already have the DNest4 Python package on your system,
follow the instructions [here](https://github.com/eggplantbren/DNest4)
(just the bits about Python) to install it.

