## Lensing2

Doing gravitational lensing again!
Why the 2? Because I had an old repository called _Lensing_, and _Lensing2_
was the best I could come up with when I decided to rewrite it from scratch.

(c) 2014 -- 2018 Brendon J. Brewer.

`fastell.f` is courtesy of Rennan Barkana and is included with permission.

## License

The contents of this repository are licensed under the GNU General Public
Licence, version 3. See the LICENCE file for details.

## Acknowledgements

This work is supported by a Marsden Fast-Start grant
from the Royal Society of New Zealand. If you use this software, please cite
the following paper, which details the methodology:

Brewer, Brendon J., David Huijser, and Geraint F. Lewis. "Trans-dimensional Bayesian inference for gravitational lens substructures." Monthly Notices of the Royal Astronomical Society 455, no. 2 (2015): 1819-1829.

You can find the paper for free on the arxiv:
https://arxiv.org/abs/1508.00662

## Dependencies

* Armadillo (http://arma.sourceforge.net/)
* yaml-cpp (https://github.com/jbeder/yaml-cpp)

You can probably get these from your operating system's package manager.
You'll also need the C++ header files, which are sometimes put into a
separate package with the suffix -dev or -devel.

You'll also need git to obtain the source code, and non-ancient versions of
g++ and gfortran for it to compile properly.

There are some associated Python scripts as well, which use Python 3
(you can try Python 2, but no guarantees), numpy,
matplotlib, yaml, and the DNest4 python package. Anaconda's distribution of
Python 3 should work well.

## Downloading and compiling

First, clone the repository recursively:
```
git clone --recursive https://github.com/eggplantbren/Lensing2
```

Then compile all the Fortran and C++:
```
cd Lensing2/src
make
```

If you don't already have the DNest4 Python package on your system,
follow the instructions [here](https://github.com/eggplantbren/DNest4)
(just the bits about Python) to install it.

## Running the example data

To run Lensing2 on the example data using 8 threads (recommended), use

```
./Lensing2 -t 8
```

Lensing2 will run and you will see DNest4 output in the terminal. DNest4 output
will also be written to some text files. Lens modelling is expensive and the
demo image is non-trivial, so give it an hour or so of runtime before
expecting anything interesting from the postprocessing. It's harmless to try
it at any time, though. If you want to try some more aggressive numerical
settings, try this instead:

```
./Lensing2 -t 8 -o OPTIONS_AGGRESSIVE
```

## Postprocessing

Lensing2 will run for a long time. The longer you run it, the more reliable the
output will be (i.e., more posterior samples will have been generated, as
long as everything went to plan).
You can manually terminate it when you like, or you can
do the postprocessing without terminating the main process.

The postprocessing will require
that you've installed the DNest4 Python package
(see [here](https://github.com/eggplantbren/DNest4) for instructions).
Simply invoke

```
python showresults.py
```

This will generate a bunch of output plots
(first the three canonical DNest4 output plots, then lensing stuff).
As you close each plot, more will appear.
Posterior samples will also be saved in a text file
`posterior_sample.txt`. Later I will write a script
to convert the posterior samples to another format for greater convenience,
so you won't have to worry too much about what's in what column. For the
time being, there is at least a header in `posterior_sample.txt` telling you
what everything is.

## The modelling assumptions

Some of the model assumptions have changed a bit since the paper was published.
Ask me for details. Also, some things are still hard-coded (such as the maximum
number of source and lens blobs allowed).

# Running other datasets
To run other datasets, inspect run_files.yaml to see how to set up a run.
You need to provide a YAML file of metadata (see Data/mock_metadata.yaml for
an example), and plain text files with the image, the sigma map, and the PSF.
Let me know if anything is unclear.

Images provided to Lensing2 should be clear of any features such as the lens
galaxy, any non-linear background, etc. You can also use the sigma map to
mask out any non-modellable features by settings those pixels to a very high
standard deviation (> 1E100, and the plotting scripts will treat those pixels
as having been totally masked).

