// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include <bio/io/detail/views_eager_split.hpp>

using seqan3::operator""_dna4;

TEST(view_eager_split, basic)
{
    {
        std::string s  = "FOO|BAR|BAX|BAT";
        auto        v  = s | bio::detail::eager_split('|');
        auto        it = v.begin();

        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "FOO");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAR");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAX");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAT");

        ++it;
        ASSERT_TRUE(it == v.end());
    }

    {
        std::string s  = "|FOO||BAR|BAX|BAT||";
        auto        v  = s | bio::detail::eager_split('|');
        auto        it = v.begin();

        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "FOO");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAR");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAX");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAT");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "");

        ++it;
        ASSERT_TRUE(it == v.end());
    }
}

TEST(view_eager_split, quotes)
{
    std::string s = "FOO,BAR\",BAX,BAT\",BAZ\",BA\"";

    {
        /* ignores quotation marks by default */
        auto v  = s | bio::detail::eager_split(',');
        auto it = v.begin();

        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "FOO");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAR\"");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAX");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAT\"");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAZ\"");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BA\"");

        ++it;
        ASSERT_TRUE(it == v.end());
    }

    {
        /* skips the delimeter inside quotation marks */
        auto v  = s | bio::detail::eager_split(',', true);
        auto it = v.begin();

        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "FOO");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAR\",BAX,BAT\"");

        ++it;
        ASSERT_FALSE(it == v.end());
        EXPECT_EQ(*it, "BAZ\",BA\"");

        ++it;
        ASSERT_TRUE(it == v.end());
    }
}

TEST(view_eager_split, concepts)
{
    std::string s = "FOO|BAR|BAX|BAT";
    auto        v = s | bio::detail::eager_split('|');

    EXPECT_TRUE(std::ranges::input_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v)>); // could be true but currently isn't
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v)>);

    EXPECT_TRUE(std::ranges::view<decltype(v)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v), std::string_view>));
}
