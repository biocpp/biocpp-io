// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>

#include <bio/io/var_io/header.hpp>

#include "../format/vcf_data.hpp"

// using namespace bio::io::literals;
using namespace std::literals;

TEST(var_io_header, spec_from_text)
{
    using svpair = std::pair<std::string const, std::string>;
    bio::io::var_io::header hdr{example_from_spec_header};

    EXPECT_EQ(hdr.file_format, "VCFv4.3");

    // filters
    ASSERT_EQ(hdr.filters.size(), 3ull);

    // filter 0
    EXPECT_EQ(hdr.filters[0].id, "PASS");
    EXPECT_EQ(hdr.filters[0].description, "\"All filters passed\"");
    EXPECT_EQ(hdr.filters[0].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.filters[0].idx, 0);

    // filter 1
    EXPECT_EQ(hdr.filters[1].id, "q10");
    EXPECT_EQ(hdr.filters[1].description, "\"Quality below 10\"");
    EXPECT_EQ(hdr.filters[1].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.filters[1].idx, 7);

    // filter 2
    EXPECT_EQ(hdr.filters[2].id, "s50");
    EXPECT_EQ(hdr.filters[2].description, "\"Less than 50% of samples have data\"");
    EXPECT_EQ(hdr.filters[2].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.filters[2].idx, 8);

    // infos
    ASSERT_EQ(hdr.infos.size(), 6ull);

    // info 0
    EXPECT_EQ(hdr.infos[0].id, "NS");
    EXPECT_EQ(hdr.infos[0].number, 1);
    EXPECT_EQ(hdr.infos[0].type_id, bio::io::var_io::value_type_id::int32);
    EXPECT_EQ(hdr.infos[0].description, "\"Number of Samples With Data\"");
    EXPECT_EQ(hdr.infos[0].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos[0].idx, 1);

    // info 1
    EXPECT_EQ(hdr.infos[1].id, "DP");
    EXPECT_EQ(hdr.infos[1].number, 1);
    EXPECT_EQ(hdr.infos[1].type_id, bio::io::var_io::value_type_id::int32);
    EXPECT_EQ(hdr.infos[1].description, "\"Total Depth\"");
    EXPECT_EQ(hdr.infos[1].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos[1].idx, 2);

    // info 2
    EXPECT_EQ(hdr.infos[2].id, "AF");
    EXPECT_EQ(hdr.infos[2].number, bio::io::var_io::header_number::A);
    EXPECT_EQ(hdr.infos[2].type_id, bio::io::var_io::value_type_id::vector_of_float32);
    EXPECT_EQ(hdr.infos[2].description, "\"Allele Frequency\"");
    EXPECT_EQ(hdr.infos[2].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos[2].idx, 3);

    // info 3
    EXPECT_EQ(hdr.infos[3].id, "AA");
    EXPECT_EQ(hdr.infos[3].number, 1);
    EXPECT_EQ(hdr.infos[3].type_id, bio::io::var_io::value_type_id::string);
    EXPECT_EQ(hdr.infos[3].description, "\"Ancestral Allele\"");
    EXPECT_EQ(hdr.infos[3].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos[3].idx, 4);

    // info 4
    EXPECT_EQ(hdr.infos[4].id, "DB");
    EXPECT_EQ(hdr.infos[4].number, 0);
    EXPECT_EQ(hdr.infos[4].type_id, bio::io::var_io::value_type_id::flag);
    EXPECT_EQ(hdr.infos[4].description, "\"dbSNP membership, build 129\"");
    EXPECT_EQ(hdr.infos[4].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos[4].idx, 5);

    // info 5
    EXPECT_EQ(hdr.infos[5].id, "H2");
    EXPECT_EQ(hdr.infos[5].number, 0);
    EXPECT_EQ(hdr.infos[5].type_id, bio::io::var_io::value_type_id::flag);
    EXPECT_EQ(hdr.infos[5].description, "\"HapMap2 membership\"");
    EXPECT_EQ(hdr.infos[5].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos[5].idx, 6);

    // formats
    ASSERT_EQ(hdr.formats.size(), 4ull);

    // format 0
    EXPECT_EQ(hdr.formats[0].id, "GT");
    EXPECT_EQ(hdr.formats[0].number, 1);
    EXPECT_EQ(hdr.formats[0].type_id, bio::io::var_io::value_type_id::string);
    EXPECT_EQ(hdr.formats[0].description, "\"Genotype\"");
    EXPECT_EQ(hdr.formats[0].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.formats[0].idx, 9);

    // format 1
    EXPECT_EQ(hdr.formats[1].id, "GQ");
    EXPECT_EQ(hdr.formats[1].number, 1);
    EXPECT_EQ(hdr.formats[1].type_id, bio::io::var_io::value_type_id::int32);
    EXPECT_EQ(hdr.formats[1].description, "\"Genotype Quality\"");
    EXPECT_EQ(hdr.formats[1].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.formats[1].idx, 10);

    // format 2
    EXPECT_EQ(hdr.formats[2].id, "DP");
    EXPECT_EQ(hdr.formats[2].number, 1);
    EXPECT_EQ(hdr.formats[2].type_id, bio::io::var_io::value_type_id::int32);
    EXPECT_EQ(hdr.formats[2].description, "\"Read Depth\"");
    EXPECT_EQ(hdr.formats[2].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.formats[2].idx, 2);

    // format 3
    EXPECT_EQ(hdr.formats[3].id, "HQ");
    EXPECT_EQ(hdr.formats[3].number, 2);
    EXPECT_EQ(hdr.formats[3].type_id, bio::io::var_io::value_type_id::vector_of_int32);
    EXPECT_EQ(hdr.formats[3].description, "\"Haplotype Quality\"");
    EXPECT_EQ(hdr.formats[3].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.formats[3].idx, 11);

    // contigs
    ASSERT_EQ(hdr.contigs.size(), 1ull);
    EXPECT_EQ(hdr.contigs[0].id, "20");
    EXPECT_EQ(hdr.contigs[0].length, 62435964);
    ASSERT_EQ(hdr.contigs[0].other_fields.size(), 4ull);
    auto it = hdr.contigs[0].other_fields.begin();
    EXPECT_TRUE(*it == (svpair{"assembly", "B36"}));
    EXPECT_TRUE(*++it == (svpair{"md5", "f126cdf8a6e0c7f379d618ff66beb2da"}));
    EXPECT_TRUE(*++it == (svpair{"species", "\"Homo sapiens\""}));
    EXPECT_TRUE(*++it == (svpair{"taxonomy", "x"}));

    // other lines in header
    std::vector<std::string_view> other_lines_cmp{"fileDate=20090805",
                                                  "source=myImputationProgramV3.1",
                                                  "reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta",
                                                  "phasing=partial"};
    ASSERT_EQ(hdr.other_lines.size(), other_lines_cmp.size());
    for (size_t i = 0; i < other_lines_cmp.size(); ++i)
        EXPECT_EQ(hdr.other_lines[i], other_lines_cmp[i]);

    // columns labels
    std::vector<std::string_view> column_labels_cmp{"CHROM",
                                                    "POS",
                                                    "ID",
                                                    "REF",
                                                    "ALT",
                                                    "QUAL",
                                                    "FILTER",
                                                    "INFO",
                                                    "FORMAT",
                                                    "NA00001",
                                                    "NA00002",
                                                    "NA00003"};
    ASSERT_EQ(hdr.column_labels.size(), column_labels_cmp.size());
    for (size_t i = 0; i < column_labels_cmp.size(); ++i)
        EXPECT_EQ(hdr.column_labels[i], column_labels_cmp[i]);

    // TODO check hash-tables

    EXPECT_EQ(hdr.to_plaintext(), example_from_spec_header_regenerated);
    EXPECT_EQ(hdr.to_plaintext_without_idx(), example_from_spec_header_regenerated_no_IDX);
}

// TODO add test with mixed-order header entries, check that IDX is correct
