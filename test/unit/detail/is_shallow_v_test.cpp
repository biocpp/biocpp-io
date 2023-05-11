// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2022, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/meta/tuple.hpp>
#include <bio/ranges/container/dictionary.hpp>
#include <bio/test/expect_range_eq.hpp>
#include <bio/test/tmp_directory.hpp>
#include <bio/test/tmp_filename.hpp>

#include <bio/io/detail/utility.hpp>

TEST(is_shallow_v, basic)
{
    /* regular types and references */
    EXPECT_FALSE(bio::io::detail::is_shallow_v<int>);
    EXPECT_TRUE(bio::io::detail::is_shallow_v<int &>);

    /* tuples */
    EXPECT_FALSE((bio::io::detail::is_shallow_v<std::tuple<int, float>>));
    EXPECT_TRUE((bio::io::detail::is_shallow_v<std::tuple<int &, float>>));
    EXPECT_FALSE((bio::io::detail::is_shallow_v<std::pair<int, float>>));
    EXPECT_TRUE((bio::io::detail::is_shallow_v<std::pair<int &, float>>));
    EXPECT_FALSE((bio::io::detail::is_shallow_v<bio::meta::tuple<int, float>>));
    EXPECT_TRUE((bio::io::detail::is_shallow_v<bio::meta::tuple<int &, float>>));

    /* ranges */
    EXPECT_FALSE((bio::io::detail::is_shallow_v<std::string>));
    EXPECT_TRUE((bio::io::detail::is_shallow_v<std::string_view>));
    EXPECT_FALSE((bio::io::detail::is_shallow_v<std::vector<int>>));
    EXPECT_TRUE((bio::io::detail::is_shallow_v<std::vector<std::string_view>>));
    EXPECT_FALSE((bio::io::detail::is_shallow_v<bio::ranges::dictionary<std::string, std::string>>));
    EXPECT_TRUE((bio::io::detail::is_shallow_v<bio::ranges::dictionary<std::string_view, std::string>>));
}
