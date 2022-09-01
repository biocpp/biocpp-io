// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#define BIOCPP_IO_NO_SYNC_CHECK 1

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>
#include <bio/test/tmp_filename.hpp>

#include <bio/io/var_io/writer.hpp>

// this test would throw an exception if the macro above was not defined
TEST(biocpp_io_issue_53, standalone)
{
    using namespace bio::alphabet::literals;

    bio::io::var_io::header hdr{};
    hdr.file_format   = "VCFv4.3";
    hdr.column_labels = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "x"};

    bio::io::var_io::writer writer{std::cout, bio::io::vcf{}};
    writer.set_header(hdr);
    bio::io::var_io::default_record<> record{};
    record.chrom()     = "chr1";
    record.pos()       = 11111;
    record.id()        = "test";
    record.ref()       = "ATC"_dna5;
    record.alt()       = {"AGC", "A"};
    record.qual()      = 1.F;
    record.filter()    = {"PASS"};
    record.genotypes() = {};
    record.info()      = {};
    writer.push_back(record);
}
