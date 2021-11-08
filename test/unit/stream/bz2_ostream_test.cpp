// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/stream/compression.hpp>
#include <bio/stream/detail/bz2_ostream.hpp>
#include <bio/stream/detail/make_stream.hpp>

#include "data.hpp"
#include "ostream_test_template.hpp"

TEST(bz2_ostream, regular)
{
    regular<bio::compression_format::bz2>();
}

TEST(bz2_ostream, type_erased)
{
    type_erased<bio::compression_format::bz2>();
}
