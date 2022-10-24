// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2022, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>
#include <bio/test/tmp_directory.hpp>
#include <bio/test/tmp_filename.hpp>

#include <bio/io/detail/index_tabix.hpp>
#include <bio/io/stream/detail/fast_streambuf_iterator.hpp>

#ifndef BIOCPP_IO_DATA_DIR
#    error "BIOCPP_IO_DATA_DIR not defined. This is required."
#endif

TEST(index_tabix, read_write)
{
    std::filesystem::path input = BIOCPP_IO_DATA_DIR;
    input /= "../format/1000G_chr10_sample.vcf.gz.tbi";

    bio::test::tmp_filename output{"out.tbi"};

    using stream_rng_t =
      std::ranges::subrange<bio::io::detail::fast_istreambuf_iterator<char>, std::default_sentinel_t>;

    {
        bio::io::detail::tabix_index idx{};
        idx.read(input);
        idx.write(output.get_path());
    }

    // verify by comparing the decompressed contents
    {
        bio::io::transparent_istream input_f{input};
        stream_rng_t                 input_s{bio::io::detail::fast_istreambuf_iterator<char>{input_f}, {}};

        bio::io::transparent_istream output_f{output.get_path()};
        stream_rng_t                 output_s{bio::io::detail::fast_istreambuf_iterator<char>{output_f}, {}};

        EXPECT_TRUE((std::ranges::equal(input_s, output_s)));
    }
}

// TODO explicit tests for reg2chunks? Is implicitly tests by indexed var tests
