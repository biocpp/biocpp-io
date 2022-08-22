// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>
#include <bio/test/tmp_filename.hpp>

#include <bio/io/format/bcf_input_handler.hpp>
#include <bio/io/var_io/reader.hpp>

#include "bcf_data.hpp"
#include "vcf_data.hpp"

//-----------------------------------------------------------------------------
// low level stuff
//-----------------------------------------------------------------------------

TEST(bcf, iterator)
{
    std::istringstream           istream{static_cast<std::string>(example_from_spec_bcf)};
    bio::io::transparent_istream str{istream};

    bio::io::detail::bcf_input_iterator it{str};

    EXPECT_EQ(it.header.text, example_from_spec_bcf_header);

    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 91ull);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 76ull);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 98ull);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 71ull);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 80ull);

    ++it;
    EXPECT_TRUE(it == std::default_sentinel);
}

TEST(bcf, iterator_underflow)
{
    bio::test::tmp_filename filename{"bcf_iterator_overflow.unbcf"};

    {
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        filecreator.write(example_from_spec_bcf_unbgzf.data(), example_from_spec_bcf_unbgzf.size());
    }

    // buffers size of 50 results in the header and all records not fitting into a buffer
    // this tests the iterators behaviour to correctly cache parts of old buffer
    bio::io::transparent_istream str{filename.get_path(), {.buffer1_size = 50}};

    bio::io::detail::bcf_input_iterator it{str};

    EXPECT_EQ(it.header.text, example_from_spec_bcf_header);

    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 91ull);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 76ull);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 98ull);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 71ull);

    ++it;
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_EQ(it->first.size(), 80ull);

    ++it;
    EXPECT_TRUE(it == std::default_sentinel);
}

//-----------------------------------------------------------------------------
// actual parsing
//-----------------------------------------------------------------------------

enum class style
{
    def,
    bcf
};

template <style s, bio::io::ownership own>
void field_types()
{
    std::istringstream           istr{std::string{example_from_spec_bcf}};
    bio::io::transparent_istream str{istr};

    bio::io::format_input_handler<bio::io::bcf> handler{str, bio::io::var_io::reader_options{}};

    bio::io::var_io::record_private_data priv{&handler.get_header()};

    using fields_t = std::conditional_t<s == style::def,
                                        decltype(bio::io::var_io::field_types<own>),
                                        decltype(bio::io::var_io::field_types_bcf_style<own>)>;
    using record_t = bio::io::record<decltype(bio::io::var_io::default_field_ids), fields_t>;

    using int_t       = int8_t;
    using vec_t       = bio::ranges::concatenated_sequences<std::vector<int_t>>;
    constexpr auto mv = bio::io::var_io::missing_value<int_t>;

    std::vector<record_t> recs;

    if constexpr (s == style::def)
        recs = example_records_default_style<own, int_t>();
    else
        recs = example_records_bcf_style<own, int_t>();

    // this workaround is pending clarification in https://github.com/samtools/hts-specs/issues/593
    std::get<vec_t>(bio::io::detail::get_second(recs[1].genotypes().back())).push_back(std::vector{mv});
    std::get<vec_t>(bio::io::detail::get_second(recs[2].genotypes().back())).push_back(std::vector{mv});
    std::get<vec_t>(bio::io::detail::get_second(recs[3].genotypes().back())).push_back(std::vector{mv});

    for (auto & rec : recs)
        get<bio::io::field::_private>(rec) = priv;

    record_t rec;

    handler.parse_next_record_into(rec);
    get<bio::io::field::_private>(rec).raw_record  = nullptr;
    get<bio::io::field::_private>(rec).record_core = nullptr;
    EXPECT_EQ(rec, recs[0]);

    handler.parse_next_record_into(rec);
    get<bio::io::field::_private>(rec).raw_record  = nullptr;
    get<bio::io::field::_private>(rec).record_core = nullptr;
    EXPECT_EQ(rec, recs[1]);

    handler.parse_next_record_into(rec);
    get<bio::io::field::_private>(rec).raw_record  = nullptr;
    get<bio::io::field::_private>(rec).record_core = nullptr;
    EXPECT_EQ(rec, recs[2]);

    handler.parse_next_record_into(rec);
    get<bio::io::field::_private>(rec).raw_record  = nullptr;
    get<bio::io::field::_private>(rec).record_core = nullptr;
    EXPECT_EQ(rec, recs[3]);

    handler.parse_next_record_into(rec);
    get<bio::io::field::_private>(rec).raw_record  = nullptr;
    get<bio::io::field::_private>(rec).record_core = nullptr;
    EXPECT_EQ(rec, recs[4]);
}

TEST(bcf, field_types_default_style_shallow)
{
    field_types<style::def, bio::io::ownership::shallow>();
}

TEST(bcf, field_types_default_style_deep)
{
    field_types<style::def, bio::io::ownership::deep>();
}

TEST(bcf, field_types_bcf_style_shallow)
{
    field_types<style::bcf, bio::io::ownership::shallow>();
}

TEST(bcf, field_types_bcf_style_deep)
{
    field_types<style::bcf, bio::io::ownership::deep>();
}
