# Welcome {#mainpage}

Welcome to the documentation of the B.I.O. library.
This web-site contains the API reference (documentation of our interfaces) and more elaborate Tutorials and
How-Tos.

### Some notes on using this documentation

We use [doxygen](https://doxygen.nl) to generate our documentation.
It may not be the most beautiful system, but it works quite well in practice.
If you spot any dead links in the documentation, please open an issue at our bug-tracker (see above) or
directly submit a pull request fixing the problem.

The documentation is versioned together with the library, see https://docs.seqan.de for release-specific
documentation builds.

Since doxygen does not support many modern C++ features, some parts of the documentation may not describe
the interfaces completely. In particular, many *constraints* are only expressed verbally in the documentation of
an interface and not as part of that interface's code. Also *C++ concepts* are currently called "interfaces" throughout
the documentation.
