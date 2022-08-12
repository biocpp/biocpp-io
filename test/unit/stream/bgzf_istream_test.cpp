// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/io/stream/compression.hpp>
#include <bio/io/stream/detail/bgzf_istream.hpp>
#include <bio/io/stream/detail/make_stream.hpp>

#include "data.hpp"
#include "istream_test_template.hpp"

TEST(bgzf_istream, regular)
{
    regular<bio::compression_format::bgzf>();
}

TEST(bgzf_istream, type_erased)
{
    type_erased<bio::compression_format::bgzf>();
}
