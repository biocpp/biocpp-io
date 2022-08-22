// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>

#include <bio/io/format/vcf_input_handler.hpp>
#include <bio/io/var_io/reader.hpp>

#include "vcf_data.hpp"

enum class style
{
    def,
    vcf,
    bcf
};

template <style s, bio::io::ownership own>
void field_types()
{
    std::istringstream istr{std::string{example_from_spec}};

    bio::io::format_input_handler<bio::io::vcf> handler{istr, bio::io::var_io::reader_options{}};

    bio::io::var_io::record_private_data priv{&handler.get_header()};

    using fields_t = std::conditional_t<s == style::def,
                                        decltype(bio::io::var_io::field_types<own>),
                                        decltype(bio::io::var_io::field_types_bcf_style<own>)>;
    using record_t = bio::io::record<decltype(bio::io::var_io::default_field_ids), fields_t>;

    std::vector<record_t> recs;

    if constexpr (s == style::def)
        recs = example_records_default_style<own>();
    else
        recs = example_records_bcf_style<own>();

    for (auto & rec : recs)
        get<bio::io::field::_private>(rec) = priv;

    record_t rec;

    handler.parse_next_record_into(rec);
    EXPECT_EQ(rec, recs[0]);

    handler.parse_next_record_into(rec);
    EXPECT_EQ(rec, recs[1]);

    handler.parse_next_record_into(rec);
    EXPECT_EQ(rec, recs[2]);

    handler.parse_next_record_into(rec);
    EXPECT_EQ(rec, recs[3]);

    handler.parse_next_record_into(rec);
    EXPECT_EQ(rec, recs[4]);
}

TEST(vcf, field_types_default_style_shallow)
{
    field_types<style::def, bio::io::ownership::shallow>();
}

TEST(vcf, field_types_default_style_deep)
{
    field_types<style::def, bio::io::ownership::deep>();
}

TEST(vcf, field_types_bcf_style_shallow)
{
    field_types<style::bcf, bio::io::ownership::shallow>();
}

TEST(vcf, field_types_bcf_style_deep)
{
    field_types<style::bcf, bio::io::ownership::deep>();
}

TEST(vcf, incomplete_header)
{
    /* This test checks the ability of VCF input to "learn" missing contigs, infos, formats... while reading the file.
     * The header is updated while reading, so this test checks the respective header code, too, i.e.
     * BOTH the records are verified and the parts of the header that are supposed to have changed */

    std::istringstream istr{incomplete_header_before + example_from_spec_records};

    using record_t =
      bio::io::record<decltype(bio::io::var_io::default_field_ids), decltype(bio::io::var_io::field_types<>)>;

    bio::io::format_input_handler<bio::io::vcf> handler{istr, bio::io::var_io::reader_options{.print_warnings = false}};

    bio::io::var_io::record_private_data priv{&handler.get_header()};

    bio::io::var_io::header const & hdr = handler.get_header();

    auto recs = example_records_default_style<bio::io::ownership::shallow>();

    for (auto & rec : recs)
        get<bio::io::field::_private>(rec) = priv;

    EXPECT_EQ(hdr.to_plaintext(), incomplete_header_before);

    bio::io::var_io::header::filter_t filter_compare;
    bio::io::var_io::header::info_t   info_compare;
    bio::io::var_io::header::format_t format_compare;
    record_t                          rec;

    /* FIRST RECORD */
    ASSERT_EQ(hdr.contigs.size(), 0ull);
    ASSERT_EQ(hdr.infos.size(), 1ull);
    ASSERT_EQ(hdr.formats.size(), 1ull);

    handler.parse_next_record_into(rec); // add contigs, infos and formats to header
    EXPECT_EQ(rec, recs[0]);

    ASSERT_EQ(hdr.contigs.size(), 1ull);
    ASSERT_EQ(hdr.infos.size(), 5ull);
    ASSERT_EQ(hdr.formats.size(), 4ull);

    EXPECT_EQ(hdr.contigs[0].id, "20");
    EXPECT_EQ(hdr.contigs[0].idx, 0);

    info_compare     = bio::io::var_io::reserved_infos.at("DP");
    info_compare.idx = 3;
    EXPECT_TRUE(hdr.infos[1] == info_compare);

    info_compare     = bio::io::var_io::reserved_infos.at("AF");
    info_compare.idx = 4;
    EXPECT_TRUE(hdr.infos[2] == info_compare);

    info_compare     = bio::io::var_io::reserved_infos.at("DB");
    info_compare.idx = 5;
    EXPECT_TRUE(hdr.infos[3] == info_compare);

    info_compare     = bio::io::var_io::reserved_infos.at("H2");
    info_compare.idx = 6;
    EXPECT_TRUE(hdr.infos[4] == info_compare);

    format_compare     = bio::io::var_io::reserved_formats.at("GQ");
    format_compare.idx = 7;
    EXPECT_TRUE(hdr.formats[1] == format_compare);

    format_compare     = bio::io::var_io::reserved_formats.at("DP");
    format_compare.idx = 3;
    EXPECT_TRUE(hdr.formats[2] == format_compare);

    format_compare     = bio::io::var_io::reserved_formats.at("HQ");
    format_compare.idx = 8;
    EXPECT_TRUE(hdr.formats[3] == format_compare);

    /* SECOND RECORD */
    ASSERT_EQ(hdr.filters.size(), 1ull);

    handler.parse_next_record_into(rec); // add filter to header
    EXPECT_EQ(rec, recs[1]);

    ASSERT_EQ(hdr.filters.size(), 2ull);
    EXPECT_EQ(hdr.filters[1].id, "q10");
    EXPECT_EQ(hdr.filters[1].description, "\"Automatically added by SeqAn3.\"");
    EXPECT_EQ(hdr.filters[1].idx, 9);

    /* THIRD RECORD */
    ASSERT_EQ(hdr.infos.size(), 5ull);

    handler.parse_next_record_into(rec); // one new info added here
    EXPECT_EQ(rec, recs[2]);

    ASSERT_EQ(hdr.infos.size(), 6ull);

    info_compare     = bio::io::var_io::reserved_infos.at("AA");
    info_compare.idx = 10;
    EXPECT_TRUE(hdr.infos[5] == info_compare);

    /* fourth and fifth don't add anything */

    EXPECT_EQ(hdr.to_plaintext(), incomplete_header_after);
}
