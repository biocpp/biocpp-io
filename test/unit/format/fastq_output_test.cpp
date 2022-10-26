// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <sstream>

#include <gtest/gtest.h>

#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/quality/phred42.hpp>
#include <bio/meta/tag/vtag.hpp>
#include <bio/test/expect_range_eq.hpp>
#include <bio/test/expect_same_type.hpp>

#include <bio/io/format/fastq_output_handler.hpp>

#include "fastq_data.hpp"

using namespace bio::alphabet::literals;
using namespace bio::meta::literals;
using std::literals::string_view_literals::operator""sv;

template <bool deep, typename salph_t, typename qalph_t>
std::string do_test(writer_options opt)
{
    std::ostringstream ostr{};

    {
        bio::io::format_output_handler<bio::io::fastq> handler{ostr, opt};

        auto recs = example_records<deep, salph_t, qalph_t>();

        for (auto & rec : recs)
            handler.write_record(rec);
    }

    return ostr.str();
}

template <typename T>
class fastq_write : public ::testing::Test
{};

using salphs = bio::meta::type_list<char, bio::alphabet::dna5>;
using qalphs = bio::meta::type_list<char, bio::alphabet::phred42>;

using numbers = ::testing::Types<bio::meta::vtag_t<0>, bio::meta::vtag_t<1>>;

TYPED_TEST_SUITE(fastq_write, numbers, );

TYPED_TEST(fastq_write, do_it)
{
    using salph_t       = bio::meta::list_traits::at<TypeParam::first_value, salphs>;
    using qalph_t       = bio::meta::list_traits::at<TypeParam::first_value, qalphs>;
    constexpr bool deep = TypeParam::first_value % 2;

    std::string s = do_test<deep, salph_t, qalph_t>(writer_options{});

    EXPECT_EQ(s, fastq_default_output);
}

TYPED_TEST(fastq_write, double_id)
{
    using salph_t       = bio::meta::list_traits::at<TypeParam::first_value, salphs>;
    using qalph_t       = bio::meta::list_traits::at<TypeParam::first_value, qalphs>;
    constexpr bool deep = TypeParam::first_value % 2;

    std::string s = do_test<deep, salph_t, qalph_t>(writer_options{.double_id = true});

    EXPECT_EQ(s, fastq_default_output_double_id);
}

// ----------------------------------------------------------------------------
// failure
// ----------------------------------------------------------------------------

TEST(fastq_write_manual, size_mismatch)
{
    std::ostringstream                             ostr{};
    bio::io::format_output_handler<bio::io::fastq> handler{ostr};
    auto recs = example_records<1, bio::alphabet::dna5, bio::alphabet::phred42>();
    recs[0].qual.resize(2);
    EXPECT_THROW(handler.write_record(recs[0]), bio::io::format_error);
}
