# Welcome {#mainpage}

Welcome to the documentation of the B.I.O. library.
This web-site contains the API reference (documentation of our interfaces) and more elaborate Tutorials and
How-Tos.


## Overview


### General IO Utilities

|                           |    Description                                                              |
|---------------------------|-----------------------------------------------------------------------------|
| bio::transparent_istream  | Like std::istream / std::ifstream, but automatically decompresses           |
| bio::transparent_ostream  | Like std::ostream / std::ofstream, but automatically compresses             |


The transparent streams can be used in place of the standard library streams. They transparently add/remove
compressions such as GZip, BZip2 and BGZip.


### Readers and Writers


| Reader                    | Writer                |    Description                                                |
|---------------------------|-----------------------|---------------------------------------------------------------|
| bio::plain_io::reader     | bio::plain_io::writer | Plaintext files, CSV, TSV; simple VCF or SAM                  |
| TODO                      | TODO                  | SAM, BAM and CRAM                                             |
| bio::seq_io::reader       | TODO                  | FastA and FastQ files; also "sequence access" to SAM/BAM/CRAM |
| bio::var_io::reader       | bio::var_io::writer   | VCF and BCF                                                   |


B.I.O. offers several modules for typical IO-"use cases". These "use cases" are implemented as abstractions over actual
file formats or layouts, so instead of "reading a FastA file", you create a bio::seq_io::reader that supports
reading various different formats and detects the format automatically.

This has the advantage that your application doesn't have to handle any special cases and seemlessly supports the
given formats. It also means that you can easily convert between them.

The actual implementation for format reading/writing is in \link format \endlink. Most users never need to access
the `*_input_handler` and `*_output_handler` classes.


## Some notes on using this documentation

We use [doxygen](https://doxygen.nl) to generate our documentation.
It may not be the most beautiful system, but it works quite well in practice.
If you spot any dead links in the documentation, please open an issue at our bug-tracker or
directly submit a pull request fixing the problem.

The documentation is versioned together with the library, see https://docs.seqan.de for release-specific
documentation builds.

Since doxygen does not support many modern C++ features, some parts of the documentation may not describe
the interfaces completely. In particular, many *constraints* are only expressed verbally in the documentation of
an interface and not as part of that interface's code. Also *C++ concepts* are currently called "interfaces" in many
parts of the documentation.
