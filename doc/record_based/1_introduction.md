# Introduction {#record_based_intro}

Most files in bioinformatics are comprised of *records*, i.e. multiple, indepedent entries that each consist of one or
more *fields*.
For example, a FastA file contains one or more sequence records that each contain an ID field and sequence field.

[TOC]

```
>myseq1
ACGT

>myseq2
GAGGA

>myseq3
ACTA
```

<center>
↓↓↓↓↓↓↓
</center>


| ID field   | sequence field |
|:----------:|:--------------:|
| "myseq1"   | "ACGT"         |
| "myseq2"   | "GAGGA"        |
| "myseq3"   | "ACTA"         |

Each line in this table is conceptionally "a record", and each file is modeled as a series of these records.
The process of "reading a file", is transforming the on-disk representation displayed above into the "abstraction" shown below.
The process of "writing a file" is the reverse.

Details on how records are defined is available here: \ref record_faq

## Readers

So called *readers* are responsible for detecting the format and decoding a file into a series of records:

\snippet test/snippet/seq_io/seq_io_reader.cpp simple_usage_file

The reader is an *input range* which is C++ terminology for "something that you can iterate over (once)".
The last bit is important, it implies that once you reach the end, the reader will be "empty". To iterate over it again, you need to recreate it.

<!-- Details on how readers are defined is available here: \ref reader_writer_faq -->

## Writers

TODO
