// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>
#include <bio/test/expect_same_type.hpp>
#include <bio/test/tmp_directory.hpp>
#include <bio/test/tmp_filename.hpp>

#include <bio/io/var/record.hpp>

#include "../format/vcf_data.hpp"

template <typename T>
class var_record : public ::testing::Test
{};

using var_record_types = ::testing::Types<bio::io::var::record_deep,
                                          bio::io::var::record_shallow,
                                          bio::io::var::record_idx_deep,
                                          bio::io::var::record_idx_shallow>;

TYPED_TEST_SUITE(var_record, var_record_types, );

TYPED_TEST(var_record, reader_requirements)
{
    EXPECT_TRUE(bio::io::var::detail::record_read_concept_checker(std::type_identity<TypeParam>{}));
}

TYPED_TEST(var_record, writer_requirements)
{
    EXPECT_TRUE(bio::io::var::detail::record_write_concept_checker(std::type_identity<TypeParam>{}));
}

TEST(var_record, concepts)
{
    EXPECT_TRUE(std::copy_constructible<bio::io::var::record_deep>);
    EXPECT_TRUE(!std::copy_constructible<bio::io::var::record_shallow>);
}

TEST(var_record_dictionary, heterogeneous_access)
{
    using namespace bio::meta::literals;

    auto records = example_records_default_style<bio::io::ownership::deep>();

    auto & record1 = std::get<0>(records);

    auto & af = record1.info["AF"];
    EXPECT_SAME_TYPE(decltype(af), bio::io::var::info_variant_deep &);

    auto & af_get = get<"AF">(af);
    EXPECT_SAME_TYPE(decltype(af_get), std::vector<float> &);

    auto & af_het = record1.info["AF"_vtag];
    EXPECT_SAME_TYPE(decltype(af_het), std::vector<float> &);
}
