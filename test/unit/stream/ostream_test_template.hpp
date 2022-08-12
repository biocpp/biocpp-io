// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <string>

#include <bio/io/misc.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <bio/io/stream/compression.hpp>
#include <bio/io/stream/detail/make_stream.hpp>
#include <bio/io/stream/transparent_ostream.hpp>

#include "data.hpp"

template <bio::io::compression_format f, typename stream_t = typename bio::io::detail::compression_stream<f>::ostream>
void regular()
{
    seqan3::test::tmp_filename filename{"ostream_test"};

    {
        std::ofstream of{filename.get_path()};
        if constexpr (std::same_as<stream_t, bio::io::transparent_ostream>)
        {
            stream_t ogzf{of, {.compression = f}};
            ogzf << uncompressed << std::flush;
        }
        else
        {
            stream_t ogzf{of};
            ogzf << uncompressed << std::flush;
        }
    }

    std::ifstream fi{filename.get_path(), std::ios::binary};
    std::string   buffer{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};

    if constexpr (f == bio::io::compression_format::bgzf)
        buffer[9] = '\x00'; // zero-out the OS byte.

    EXPECT_EQ(buffer, compressed<f>);
}

template <bio::io::compression_format f, typename stream_t = typename bio::io::detail::compression_stream<f>::ostream>
void type_erased()
{
    seqan3::test::tmp_filename filename{"ostream_test"};

    {
        std::ofstream of{filename.get_path()};

        if constexpr (std::same_as<stream_t, bio::io::transparent_ostream>)
        {
            std::unique_ptr<std::ostream> ogzf{
              new stream_t{of, {.compression = f}}
            };
            *ogzf << uncompressed << std::flush;
        }
        else
        {
            std::unique_ptr<std::ostream> ogzf{new stream_t{of}};
            *ogzf << uncompressed << std::flush;
        }
    }

    std::ifstream fi{filename.get_path(), std::ios::binary};
    std::string   buffer{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};

    if constexpr (f == bio::io::compression_format::bgzf)
        buffer[9] = '\x00'; // zero-out the OS byte.

    EXPECT_EQ(buffer, compressed<f>);
}
