// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>

#include <bio/io/var/header.hpp>

#include "../format/vcf_data.hpp"

// using namespace bio::io::literals;
using namespace std::literals;

TEST(var_header, spec_from_text)
{
    using svpair = std::pair<std::string const, std::string>;
    bio::io::var::header hdr{example_from_spec_header};

    EXPECT_EQ(hdr.file_format, "VCFv4.3");

    // filters
    ASSERT_EQ(hdr.filters.size(), 3ull);

    // filter 0
    EXPECT_EQ(hdr.filters.find("PASS") - hdr.filters.begin(), 0);
    EXPECT_EQ(hdr.filters["PASS"].description, "\"All filters passed\"");
    EXPECT_EQ(hdr.filters["PASS"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.filters["PASS"].idx, 0);

    // filter 1
    EXPECT_EQ(hdr.filters.find("q10") - hdr.filters.begin(), 1);
    EXPECT_EQ(hdr.filters["q10"].description, "\"Quality below 10\"");
    EXPECT_EQ(hdr.filters["q10"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.filters["q10"].idx, 7);

    // filter 2
    EXPECT_EQ(hdr.filters.find("s50") - hdr.filters.begin(), 2);
    EXPECT_EQ(hdr.filters["s50"].description, "\"Less than 50% of samples have data\"");
    EXPECT_EQ(hdr.filters["s50"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.filters["s50"].idx, 8);

    // infos
    ASSERT_EQ(hdr.infos.size(), 6ull);

    // info 0
    EXPECT_EQ(hdr.infos.find("NS") - hdr.infos.begin(), 0);
    EXPECT_EQ(hdr.infos["NS"].number, 1);
    EXPECT_EQ(hdr.infos["NS"].type_id, bio::io::var::value_type_id::int32);
    EXPECT_EQ(hdr.infos["NS"].description, "\"Number of Samples With Data\"");
    EXPECT_EQ(hdr.infos["NS"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos["NS"].idx, 1);

    // info 1
    EXPECT_EQ(hdr.infos.find("DP") - hdr.infos.begin(), 1);
    EXPECT_EQ(hdr.infos["DP"].number, 1);
    EXPECT_EQ(hdr.infos["DP"].type_id, bio::io::var::value_type_id::int32);
    EXPECT_EQ(hdr.infos["DP"].description, "\"Total Depth\"");
    EXPECT_EQ(hdr.infos["DP"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos["DP"].idx, 2);

    // info 2
    EXPECT_EQ(hdr.infos.find("AF") - hdr.infos.begin(), 2);
    EXPECT_EQ(hdr.infos["AF"].number, bio::io::var::header_number::A);
    EXPECT_EQ(hdr.infos["AF"].type_id, bio::io::var::value_type_id::vector_of_float32);
    EXPECT_EQ(hdr.infos["AF"].description, "\"Allele Frequency\"");
    EXPECT_EQ(hdr.infos["AF"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos["AF"].idx, 3);

    // info 3
    EXPECT_EQ(hdr.infos.find("AA") - hdr.infos.begin(), 3);
    EXPECT_EQ(hdr.infos["AA"].number, 1);
    EXPECT_EQ(hdr.infos["AA"].type_id, bio::io::var::value_type_id::string);
    EXPECT_EQ(hdr.infos["AA"].description, "\"Ancestral Allele\"");
    EXPECT_EQ(hdr.infos["AA"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos["AA"].idx, 4);

    // info 4
    EXPECT_EQ(hdr.infos.find("DB") - hdr.infos.begin(), 4);
    EXPECT_EQ(hdr.infos["DB"].number, 0);
    EXPECT_EQ(hdr.infos["DB"].type_id, bio::io::var::value_type_id::flag);
    EXPECT_EQ(hdr.infos["DB"].description, "\"dbSNP membership, build 129\"");
    EXPECT_EQ(hdr.infos["DB"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos["DB"].idx, 5);

    // info 5
    EXPECT_EQ(hdr.infos.find("H2") - hdr.infos.begin(), 5);
    EXPECT_EQ(hdr.infos["H2"].number, 0);
    EXPECT_EQ(hdr.infos["H2"].type_id, bio::io::var::value_type_id::flag);
    EXPECT_EQ(hdr.infos["H2"].description, "\"HapMap2 membership\"");
    EXPECT_EQ(hdr.infos["H2"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.infos["H2"].idx, 6);

    // formats
    ASSERT_EQ(hdr.formats.size(), 4ull);

    // format 0
    EXPECT_EQ(hdr.formats.find("GT") - hdr.formats.begin(), 0);
    EXPECT_EQ(hdr.formats["GT"].number, 1);
    EXPECT_EQ(hdr.formats["GT"].type_id, bio::io::var::value_type_id::string);
    EXPECT_EQ(hdr.formats["GT"].description, "\"Genotype\"");
    EXPECT_EQ(hdr.formats["GT"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.formats["GT"].idx, 9);

    // format 1
    EXPECT_EQ(hdr.formats.find("GQ") - hdr.formats.begin(), 1);
    EXPECT_EQ(hdr.formats["GQ"].number, 1);
    EXPECT_EQ(hdr.formats["GQ"].type_id, bio::io::var::value_type_id::int32);
    EXPECT_EQ(hdr.formats["GQ"].description, "\"Genotype Quality\"");
    EXPECT_EQ(hdr.formats["GQ"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.formats["GQ"].idx, 10);

    // format 2
    EXPECT_EQ(hdr.formats.find("DP") - hdr.formats.begin(), 2);
    EXPECT_EQ(hdr.formats["DP"].number, 1);
    EXPECT_EQ(hdr.formats["DP"].type_id, bio::io::var::value_type_id::int32);
    EXPECT_EQ(hdr.formats["DP"].description, "\"Read Depth\"");
    EXPECT_EQ(hdr.formats["DP"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.formats["DP"].idx, 2);

    // format 3
    EXPECT_EQ(hdr.formats.find("HQ") - hdr.formats.begin(), 3);
    EXPECT_EQ(hdr.formats["HQ"].number, 2);
    EXPECT_EQ(hdr.formats["HQ"].type_id, bio::io::var::value_type_id::vector_of_int32);
    EXPECT_EQ(hdr.formats["HQ"].description, "\"Haplotype Quality\"");
    EXPECT_EQ(hdr.formats["HQ"].other_fields.size(), 0ull);
    EXPECT_EQ(hdr.formats["HQ"].idx, 11);

    // contigs
    ASSERT_EQ(hdr.contigs.size(), 1ull);
    EXPECT_EQ(hdr.contigs.find("20") - hdr.contigs.begin(), 0);
    EXPECT_EQ(hdr.contigs["20"].length, 62435964);
    ASSERT_EQ(hdr.contigs["20"].other_fields.size(), 4ull);
    auto it = hdr.contigs["20"].other_fields.begin();
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
