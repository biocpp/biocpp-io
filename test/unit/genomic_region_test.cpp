// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <sstream>

#include <gtest/gtest.h>

#include <bio/genomic_region.hpp>

struct genomic_region : public ::testing::Test
{
    bio::genomic_region<> reg0{"chr20", 100, 200};
    bio::genomic_region<> reg1{"chr20", 100, 200};
    bio::genomic_region<> reg2{"chr20", 100, 150};
    bio::genomic_region<> reg3{"chr20", 150, 200};
    bio::genomic_region<> reg4{"chr20", 50, 100};
    bio::genomic_region<> reg5{"chr20", 200, 250};
    bio::genomic_region<> reg6{"chr21", 100, 200};
};

TEST_F(genomic_region, construction_deduction)
{
    bio::genomic_region reg0{"chr20", 100, 200};
    EXPECT_TRUE((std::same_as<decltype(reg0), bio::genomic_region<>>));

    bio::genomic_region<bio::ownership::deep> reg1{"chr20", 100, 200};
    EXPECT_TRUE((std::same_as<decltype(reg1), bio::genomic_region<>>));
    EXPECT_TRUE((std::same_as<decltype(reg1)::string_t, std::string>));

    bio::genomic_region<bio::ownership::shallow> reg2{"chr20", 100, 200};
    EXPECT_TRUE((std::same_as<decltype(reg2)::string_t, std::string_view>));
}

// just test some basics
TEST_F(genomic_region, default_compare)
{
    EXPECT_EQ(reg0, reg1);

    EXPECT_LT(reg0, reg3);
    EXPECT_LT(reg0, reg5);
    EXPECT_LT(reg0, reg6);
    EXPECT_LT(reg1, reg3);

    EXPECT_GT(reg0, reg4);
    EXPECT_GT(reg0, reg4);
    EXPECT_GT(reg6, reg1);
}

TEST_F(genomic_region, relative_to_region)
{
    EXPECT_TRUE(reg0.relative_to(reg1) == std::weak_ordering::equivalent);
    EXPECT_TRUE(reg0.relative_to(reg2) == std::weak_ordering::equivalent);
    EXPECT_TRUE(reg0.relative_to(reg3) == std::weak_ordering::equivalent);

    EXPECT_TRUE(reg0.relative_to(reg4) == std::weak_ordering::greater);
    EXPECT_TRUE(reg3.relative_to(reg4) == std::weak_ordering::greater);
    EXPECT_TRUE(reg5.relative_to(reg4) == std::weak_ordering::greater);

    EXPECT_TRUE(reg0.relative_to(reg5) == std::weak_ordering::less);
    EXPECT_TRUE(reg0.relative_to(reg6) == std::weak_ordering::less);
    EXPECT_TRUE(reg4.relative_to(reg5) == std::weak_ordering::less);
}

TEST_F(genomic_region, relative_to_point)
{
    EXPECT_TRUE(reg0.relative_to("chr20", 100) == std::weak_ordering::equivalent);
    EXPECT_TRUE(reg0.relative_to("chr20", 150) == std::weak_ordering::equivalent);
    EXPECT_TRUE(reg0.relative_to("chr20", 200) != std::weak_ordering::equivalent);

    EXPECT_TRUE(reg0.relative_to("chr20", 99) == std::weak_ordering::greater);
    EXPECT_TRUE(reg3.relative_to("chr19", 666) == std::weak_ordering::greater);
    EXPECT_TRUE(reg5.relative_to("", 200) == std::weak_ordering::greater);

    EXPECT_TRUE(reg0.relative_to("chr20", 200) == std::weak_ordering::less);
    EXPECT_TRUE(reg0.relative_to("chr20", 201) == std::weak_ordering::less);
    EXPECT_TRUE(reg4.relative_to("chr21", 1) == std::weak_ordering::less);
}

TEST_F(genomic_region, distance)
{
    EXPECT_TRUE(reg0.distance(reg1) == -100);
    EXPECT_TRUE(reg0.distance(reg2) == -50);
    EXPECT_TRUE(reg0.distance(reg3) == -50);

    EXPECT_TRUE(reg0.distance(reg4) == 0);
    EXPECT_TRUE(reg3.distance(reg4) == 50);
    EXPECT_TRUE(reg5.distance(reg4) == 100);

    EXPECT_TRUE(reg0.distance(reg5) == 0);
    EXPECT_TRUE(reg0.distance(reg6) == std::numeric_limits<int64_t>::max());
    EXPECT_TRUE(reg4.distance(reg5) == 100);
}
