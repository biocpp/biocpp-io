// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

// #include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <bio/format/bcf_input_handler.hpp>
#include <bio/var_io/reader.hpp>

#include "bcf_data.hpp"
#include "vcf_data.hpp"

//-----------------------------------------------------------------------------
// low level stuff
//-----------------------------------------------------------------------------

TEST(bcf, iterator)
{
    std::istringstream       istream{static_cast<std::string>(example_from_spec_bcf)};
    bio::transparent_istream str{istream};

    bio::detail::bcf_input_iterator it{str};

    EXPECT_EQ(it.header.text, example_from_spec_bcf_header);

    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 91);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 76);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 98);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 71);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 80);

    ++it;
    EXPECT_TRUE(it == std::default_sentinel);
}

TEST(bcf, iterator_underflow)
{
    seqan3::test::tmp_filename filename{"bcf_iterator_overflow.unbcf"};

    {
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        filecreator.write(example_from_spec_bcf_unbgzf.data(), example_from_spec_bcf_unbgzf.size());
    }

    // buffers size of 50 results in the header and all records not fitting into a buffer
    // this tests the iterators behaviour to correctly cache parts of old buffer
    bio::transparent_istream str{filename.get_path(), {.buffer1_size = 50}};

    bio::detail::bcf_input_iterator it{str};

    EXPECT_EQ(it.header.text, example_from_spec_bcf_header);

    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 91);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 76);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 98);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 71);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 80);

    ++it;
    EXPECT_TRUE(it == std::default_sentinel);
}

//-----------------------------------------------------------------------------
// actual parsing
//-----------------------------------------------------------------------------

enum class style
{
    def,
    vcf,
    bcf
};

template <style s, bio::ownership own>
void field_types()
{
    std::istringstream       istr{std::string{example_from_spec_bcf}};
    bio::transparent_istream str{istr};

    bio::format_input_handler<bio::bcf> handler{str, bio::var_io::reader_options{}};

    bio::var_io::record_private_data priv{&handler.get_header()};

    using fields_t = std::conditional_t<s == style::def,
                                        decltype(bio::var_io::field_types<own>),
                                        std::conditional_t<s == style::vcf,
                                                           decltype(bio::var_io::field_types_vcf_style<own>),
                                                           decltype(bio::var_io::field_types_bcf_style<own>)>>;
    using record_t = bio::record<decltype(bio::var_io::default_field_ids), fields_t>;

    using int_t       = int8_t;
    constexpr auto mv = bio::var_io::missing_value<int_t>;

    std::vector<record_t> recs;

    if constexpr (s == style::def)
        recs = example_records_default_style<own, int_t>();
    else if constexpr (s == style::vcf)
        recs = example_records_vcf_style<own, int_t>();
    else
        recs = example_records_bcf_style<own, int_t>();

    // this workaround is pending clarification in https://github.com/samtools/hts-specs/issues/593
    if constexpr (s == style::vcf)
    {
        bio::detail::get_second(recs[1].genotypes()).back().push_back(std::vector<int_t>{mv});
        bio::detail::get_second(recs[2].genotypes()).back().push_back(std::vector<int_t>{mv});
        bio::detail::get_second(recs[3].genotypes()).back().push_back(std::vector<int_t>{mv});
    }
    else
    {
        std::get<std::vector<std::vector<int_t>>>(bio::detail::get_second(recs[1].genotypes().back())).push_back({mv});
        std::get<std::vector<std::vector<int_t>>>(bio::detail::get_second(recs[2].genotypes().back())).push_back({mv});
        std::get<std::vector<std::vector<int_t>>>(bio::detail::get_second(recs[3].genotypes().back())).push_back({mv});
    }

    for (auto & rec : recs)
        get<bio::field::_private>(rec) = priv;

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

TEST(bcf, field_types_default_style_shallow)
{
    field_types<style::def, bio::ownership::shallow>();
}

TEST(bcf, field_types_default_style_deep)
{
    field_types<style::def, bio::ownership::deep>();
}

TEST(bcf, field_types_vcf_style_shallow)
{
    field_types<style::vcf, bio::ownership::shallow>();
}

TEST(bcf, field_types_vcf_style_deep)
{
    field_types<style::vcf, bio::ownership::deep>();
}

TEST(bcf, field_types_bcf_style_shallow)
{
    field_types<style::bcf, bio::ownership::shallow>();
}

TEST(bcf, field_types_bcf_style_deep)
{
    field_types<style::bcf, bio::ownership::deep>();
}
