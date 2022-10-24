# Welcome {#mainpage}

Welcome to the documentation of the BioC++ I/O library.
This web-site contains the API reference (documentation of our interfaces) and some small Tutorials and HowTos.

The I/O library requires the BioC++ core library, and it is recommended to have a look at its documentation first.


## Overview

This section contains a very short overview of the most important parts of the library.


### General IO Utilities

|                               |    Description                                                              |
|-------------------------------|-----------------------------------------------------------------------------|
| bio::io::transparent_istream  | Like std::istream / std::ifstream, but automatically decompresses           |
| bio::io::transparent_ostream  | Like std::ostream / std::ofstream, but automatically compresses             |


The transparent streams can be used in place of the standard library streams. They transparently add/remove
compressions such as GZip, BZip2 and BGZip.


### Record-based I/O


| Reader                    | Writer                    |    Description                                |
|---------------------------|---------------------------|-----------------------------------------------|
| bio::io::txt::reader      | bio::io::txt::writer      | Plaintext files, CSV, TSV; simple VCF or SAM  |
| bio::io::seq::reader      | TODO                      | FastA and FastQ files                         |
| bio::io::var::reader      | bio::io::var::writer      | VCF and BCF                                   |


The I/O library offers several modules for typical IO-"use cases". These "use cases" are implemented as abstractions over actual
file formats or layouts, so instead of "reading a FastA file", you create a bio::io::seq::reader that supports
reading various different formats and detects the format automatically.

This has the advantage that your application doesn't have to handle any special cases and seemlessly supports the
given formats. It also means that you can easily convert between them.

The actual implementation for format reading/writing is in \link format \endlink. Most users never need to access
the `*_input_handler` and `*_output_handler` classes.

Please note that SAM/BAM/CRAM support will be added at a later point.
