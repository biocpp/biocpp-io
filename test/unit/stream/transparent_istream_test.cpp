// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/io/stream/compression.hpp>
#include <bio/io/stream/detail/make_stream.hpp>
#include <bio/io/stream/transparent_istream.hpp>

#include "data.hpp"
#include "istream_test_template.hpp"

TEST(transparent_istream, regular_none)
{
    regular<bio::io::compression_format::none, bio::io::transparent_istream>();
}

TEST(transparent_istream, type_erased_none)
{
    type_erased<bio::io::compression_format::none, bio::io::transparent_istream>();
}

#if BIOCPP_IO_HAS_ZLIB
TEST(transparent_istream, regular_bgzf)
{
    regular<bio::io::compression_format::bgzf, bio::io::transparent_istream>();
}

TEST(transparent_istream, type_erased_bgzf)
{
    type_erased<bio::io::compression_format::bgzf, bio::io::transparent_istream>();
}

TEST(transparent_istream, regular_gz)
{
    regular<bio::io::compression_format::gz, bio::io::transparent_istream>();
}

TEST(transparent_istream, type_erased_gz)
{
    type_erased<bio::io::compression_format::gz, bio::io::transparent_istream>();
}
#endif

#if BIOCPP_IO_HAS_BZIP2
TEST(transparent_istream, regular_bz2)
{
    regular<bio::io::compression_format::bz2, bio::io::transparent_istream>();
}

TEST(transparent_istream, type_erased_bz2)
{
    type_erased<bio::io::compression_format::bz2, bio::io::transparent_istream>();
}
#endif
