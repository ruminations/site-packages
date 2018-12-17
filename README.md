# site-packages

- [Synopsis](#synopsis)
- [Roadmap](#roadmap)
- [Installation](#installation)
- [Package Summaries](#package-summaries)
- [License](#license)
- [Contact](#contact)

-----------
## Synopsis

This repository is my personal library of `python` and `cython` code.

----------
## Roadmap

I use the code routinely, and I repair things immediately when they do not do what I want.  This may include arbitrary design changes to the API, although most of it is pretty stable now.  So if you see something you like, it is not likely to break tomorrow.  If a significant user base evolves, I will, of course, have to move to a more formal release cycle.

---------------
## Installation

On `ubuntu` derivatives, `python2.7` add-on packages are installed in:

     /usr/local/lib/python2.7/dist-packages

and the historical directory:

     /usr/local/lib/python2.7/site-packages

is empty.  I use the latter directory to make `python` and `cython` code I authored available for use on my system.  That code is stored as `.py` and `.pyc` (`python`) or `.so` (`cython`) files, as appropriate.  This repository contains the necessary source code and setup files for creating that arrangement, stored in a tree that mirrors the destination of the object code.

If you would like to use this code in a similar fashion, clone the repository to somewhere in your home directory.  For the pure `python` code, issue these commands from your home directory:

     python -m compileall ./site-packages  # recursively compile the clone files
     cd ./site-packages                    # change current directory for searches

To compile the `cython` files, you will, of course, have to install `cython`, if you don't already have it.  Under `ubuntu` and derivatives, install the `cython` package with the commands:

     sudo apt-get update          # freshen local package lists
     sudo apt-get install cython  # will also install any dependencies

For the `cython` code, each `<name>.pyx` file has an associated `setup_<name>.py` file.  First, get a list of the `./<path>/<name>.pyx` files (currently four files):

     find -name "*.pyx"

To generate each `.so` file issue the commands:

     cd ./<path>
     python ./setup_<name>.py build_ext --inplace
     rm ./setup_<name>.py*  # no need to copy these to site-packages
     cd ~/site-packages     # move back to your home directory clone root

This sequence compiles python byte code files and cython shareable objects.

Notice that the modules, `gauss`, `space`, and `space_4D` are doubly implemented: there is a `python` version which establishes correctly functioning algorithms, and then there is a `cython` translation, which provides performance.  You will need to choose which one you prefer and delete the other from your `site-packages` clone directory.  It is not clear what `python` will load if both are present - probably the first one that it finds.

Now you are ready to copy the files to:

     /usr/local/lib/python2.7/site-packages/

Issue the commands:

     # find .py,.pyc, and .pyx files with paths; omit .pyx; pipe them to cpio
     find -name "*.py*" -regex ".*py[c]?$" | \
     # pass-through (-p) to site packages; recreate directory tree (d);
     #   preserve creation/modification times (m); silently overwrite (u);
     #   observing unbuntu ownership conventions (-R)
     sudo cpio -R root:staff -pdmu /usr/local/lib/python2.7/site-packages/
     # do the same for .so files
     find -name "*.so"  | \
     sudo cpio -R root:staff -pdmu /usr/local/lib/python2.7/site-packages/

The first time you install, put that directory in your `PYTHONPATH` by including the following line at the end of your your home directory `~/.bashrc` file:

     export PYTHONPATH=/usr/local/lib/python2.7/site-packages

Then issue the command `. .bashrc` (*notice the space:'`._.`'*) in your home directory to reread the `.bashrc` file without the need to reboot or start a new terminal session.

Finally, delete the clone directory from your home directory.

TODO: Create a make file that automates this.

--------------------
## Package Summaries

All of these modules have some lacunae.  All of them are usable.  The more fundamental ones are quite complete and have been used for various purposes for years.

### `algebra` modules

> #### `calculus.py`

> Facilitates computation of arbitrary `n`-th derivatives of standard math functions and compositions of those functions.  Provides tools for numerical integration and black body calculations.  Contains a brief and correct recursive implementation of the Faà di Bruno formula and some useful statistical tools.

> #### `mathx.pyx`

> A set of convenience functions that are not implemented in the standard `math` module.

### `data` modules

> #### `analysis.py`

> Provides a class for linearly interpolating a dataset and a class for numerically computing derivatives of an analytic base function, e.g. a linear interpolant of a dataset.

> #### `catalog.py`

> Provides a `metaclass` for a read-only singleton dictionary and various data dictionaries.

> #### `cie.py`

> Provides a table of CIE-LMS data and its transformation matrix to CIE-XYZ.  Provides rectification of that data and matrix for use with linear interpolation and trapezoidal integration such that identical CIE-XYZ data is produced.

> #### `sort.py`

> Provides a class implementation of the functionality of the python standard library `heapq.py` module.

### `svg` modules

> #### `dom.py`

> Provides classes for rendering the essential elements of the `svg` DOM in an approximately one to one direct manner.

> #### `plot.py`

> Provides more abstract classes for creating `svg` plots.

### `geometry` modules

> #### `figure.py`

> Provides classes for representing planar figures in three dimensions.

> #### `gauss.py` and `gauss.pyx`

> Re-implements the built-in `complex` number data type functionality in a way that automatically coerces the built-in data type and augments it with geometric methods.

> #### `parametric.py`

> Provides tools for differential geometry calculations on paths in space.

> #### `point.py`

> Provides classes representing 3-dimensional Cartesian geometric points as homogeneous 4-tuples, 3-dimensional affine differentials of geometric points, and quaternions.  Operator overrides adhere strictly to the semantic distinctions facilitating geometrically correct expressions of computations, and informative error messages when the expression is geometrically incorrect.

> #### `sheet.py`

> Provides classes for representing surfaces.  Operations include containment testing and intersection.

> #### `solid.py`

> Provides classes for representing solids.  Operations include containment testing.

> #### `space.pyx`

> Provides linear algebra tools for general multi-dimensional geometric modelling.

> #### `space_4D.pyx`

> Provides linear algebra tools optimized for geometric modelling specific to the 4-dimensional projective space used in computer graphics.  This `cython` implementation executes more quickly than the same semantics implemented with `numpy`.

----------
## License

See:

> [License Link](https://github.com/ruminations/Licenses#design-license)

----------
## Contact
See:

> [Contact Link](https://github.com/ruminations/Contact)

