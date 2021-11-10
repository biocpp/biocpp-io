// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/stream/compression.hpp>
#include <bio/stream/detail/make_stream.hpp>
#include <bio/stream/transparent_ostream.hpp>

#include "data.hpp"
#include "ostream_test_template.hpp"

TEST(transparent_ostream, regular_none)
{
    regular<bio::compression_format::none, bio::transparent_ostream>();
}

TEST(transparent_ostream, type_erased_none)
{
    type_erased<bio::compression_format::none, bio::transparent_ostream>();
}

#if BIO_HAS_ZLIB
TEST(transparent_ostream, regular_bgzf)
{
    regular<bio::compression_format::bgzf, bio::transparent_ostream>();
}

TEST(transparent_ostream, type_erased_bgzf)
{
    type_erased<bio::compression_format::bgzf, bio::transparent_ostream>();
}

TEST(transparent_ostream, regular_gz)
{
    regular<bio::compression_format::gz, bio::transparent_ostream>();
}

TEST(transparent_ostream, type_erased_gz)
{
    type_erased<bio::compression_format::gz, bio::transparent_ostream>();
}
#endif

#if BIO_HAS_BZIP2
TEST(transparent_ostream, regular_bz2)
{
    regular<bio::compression_format::bz2, bio::transparent_ostream>();
}

TEST(transparent_ostream, type_erased_bz2)
{
    type_erased<bio::compression_format::bz2, bio::transparent_ostream>();
}
#endif
