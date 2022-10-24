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

#include <bio/io/var_io/record.hpp>

template <typename T>
class var_io_record : public ::testing::Test
{};

using var_io_record_types = ::testing::Types<bio::io::var_io::record_default,
                                             bio::io::var_io::record_default_shallow,
                                             bio::io::var_io::record_idx,
                                             bio::io::var_io::record_idx_shallow>;

TYPED_TEST_SUITE(var_io_record, var_io_record_types, );

TYPED_TEST(var_io_record, reader_requirements)
{
    EXPECT_TRUE(bio::io::var_io::detail::record_read_concept_checker(std::type_identity<TypeParam>{}));
}
