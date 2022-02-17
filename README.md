# B.I.O. – the Biological Input/Output library

B.I.O. is a C++ library for reading and writing files in the field of Bioinformatics and in particular sequence
analysis. It provides easy-to-use interfaces for the following formats:

  * Plain I/O: plain-text, CSV, TSV, …
  * Map I/O: SAM, BAM, …
  * Seq I/O: FastA, FastQ, …
  * Var I/O: VCF, BCF, …

The primary goal of this library is to offer higher level abstractions than the C libraries typically used in this
domain (e.g. htslib) while at the same time offering an excellent performance.
It hopes to offer a modern, well-integrated design that covers most typical I/O use-cases Bioinformaticians encounter.

The library relies strongly on *Modern C++* and plays well with other Modern C++ libraries.

Please see the [online documentation](TODO) for more details.

## Current state

The library is currently under heavy development. There is no release, yet, and all interfaces are subject to change.

## Dependencies

|                   | requirement                               | version  | comment                                     |
|-------------------|-------------------------------------------|----------|---------------------------------------------|
|**compiler**       | [GCC](https://gcc.gnu.org)                | ≥ 10     | no other compiler is currently supported!   |
|**required libs**  | [SeqAn3](https://github.com/seqan/seqan3) | ≥ 3      |                                             |
|**optional libs**  | [zlib](https://github.com/madler/zlib)    | ≥ 1.2    | required for `*.gz` and `.bam` file support |
|                   | [bzip2](https://www.sourceware.org/bzip2) | ≥ 1.0    | required for `*.bz2` file support           |

## Usage

* Using the library entails no build-steps, it is header-only and can be used as-is.
* A single-header version is available (TODO).
* CMake files are provided for easy integration into applications (and automatic detection/inclusion of dependencies).
