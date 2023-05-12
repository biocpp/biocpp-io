// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2023, Hannes Hauswedell
// Copyright (c) 2020-2022, deCODE Genetics
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/biocpp/biocpp-io/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <concepts>

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>
#include <bio/test/expect_same_type.hpp>
#include <bio/test/tmp_filename.hpp>

#include <bio/io/seq/record.hpp>

TEST(seq_record, concepts)
{
    EXPECT_TRUE(std::copy_constructible<bio::io::seq::record_dna_deep>);
    EXPECT_TRUE(!std::copy_constructible<bio::io::seq::record_dna_shallow>);
}
