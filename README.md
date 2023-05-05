<img align="left" src="test/documentation/biocpp_logo.svg">

# The BioC++ Input/Output Library

This is the I/O library of the [BioC++ project](https://github.com/biocpp/biocpp). It provides the `io` module
which offers easy-to-use interfaces for the following formats:

  * **txt** formats: plain-text files (line-wise reading) and delimited files (e.g. CSV and TSV).
  * **seq** formats: FastA and FastQ.
  * **var** formats: VCF and BCF.

The primary goal of this library is to offer higher level abstractions than the C libraries typically used in this
domain (e.g. htslib) while at the same time offering an excellent performance.
It hopes to offer a modern, well-integrated design that covers most typical I/O use-cases Bioinformaticians encounter.

The library provides stand-alone implementations of the formats (and does not require htslib). Support for reading and
writing SAM/BAM/CRAM will be available in separate library until full implementations are available here.

Please see the [online documentation](https://biocpp.github.io) for more details.

**Attention:** this library is currently a work-in-progress, and interfaces are not yet stable.


## Example

Simple reading of a FastA-file that is transparently decompressed:

```cpp
bio::io::seq::reader reader{"example.fasta.gz"};

for (auto & rec : reader)
{
  fmt::print("ID:  {}\n", rec.id);
  fmt::print("Seq: {}\n", rec.seq);
}
```

Reading a variant file and writing a new one that only contains variants that "PASS":

```cpp
bio::io::var::reader reader{"example.vcf.gz"};
bio::io::var::writer writer{"example.bcf"};

for (auto & rec : reader)
  if (rec.filter.empty() || (rec.filter.size() == 1 && rec.filter[0] == "PASS"))
    writer.push_back(rec);
```
The format is transparently converted from compressed VCF to BCF if files have the respective extensions / magic
headers.

## Easy to use

  * Header-only → just drop it in your source code or include it as a git submodule!
  * The BioC++ core library is the only hard requirement.
  * No build-system and no configure steps required.
  * Optional CMake support available.
  * Integrates well with the C++20 standard library.

## Dependencies

|                   | requirement                                          | version  | comment                                     |
|-------------------|------------------------------------------------------|----------|---------------------------------------------|
|**compiler**       | [GCC](https://gcc.gnu.org)                           | ≥ 11     | no other compiler is currently supported!   |
|**required libs**  | [BioC++ core](https://github.com/biocpp/biocpp-core) | = 0.7    |                                             |
|**optional libs**  | [zlib](https://github.com/madler/zlib)               | ≥ 1.2    | required for `*.gz` and `*.bcf` support     |
|                   | [bzip2](https://www.sourceware.org/bzip2)            | ≥ 1.0    | required for `*.bz2` file support           |

## Quick-Setup

  * The recommended way to install BioC++ libraries is to follow the instructions in the meta-repository: https://github.com/biocpp/biocpp
  * CMake is recommended, because it finds dependencies automatically and sets required flags/paths.
  * Quick instructions without CMake:
    1. Clone biocpp-core and biocpp-io into e.g. `~/devel`.
    2. Assuming you have ZLib and Bzip2 in system paths, do:

```sh
g++ -O3 -DNDEBUG -Wall -Wextra -std=c++20           \
    -I ~/devel/biocpp-core/include                  \
    -I ~/devel/biocpp-io/include                    \
    -DBIOCPP_IO_HAS_ZLIB=1 -DBIOCPP_IO_HAS_BZIP2=1  \
    your_file.cpp
```
  * If you want to use {fmt}, add `-I ~/devel/fmt/include -D FMT_HEADER_ONLY=1`.

