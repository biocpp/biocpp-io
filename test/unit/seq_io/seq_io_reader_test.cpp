// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <sstream>

#include <gtest/gtest.h>

#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/test/expect_range_eq.hpp>
#include <bio/test/expect_same_type.hpp>
#include <bio/test/tmp_filename.hpp>

#include <bio/io/seq_io/reader.hpp>
#include <bio/io/stream/detail/fast_streambuf_iterator.hpp>

#include "data.hpp"

TEST(seq_io_reader, concepts)
{
    using t = bio::io::seq_io::reader<>;
    EXPECT_TRUE((std::ranges::input_range<t>));

    using ct = bio::io::seq_io::reader<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::input_range<ct>));
}

void seq_io_reader_filename_constructor(bool ext_check, auto &&... args)
{
    /* just the filename */
    {
        bio::test::tmp_filename filename{"seq_io_reader_constructor.fasta"};
        std::ofstream           filecreator{filename.get_path(), std::ios::out | std::ios::binary};

        EXPECT_NO_THROW((bio::io::seq_io::reader{filename.get_path(), std::forward<decltype(args)>(args)...}));
    }

    // correct format check is done by tests of that format

    /* non-existent file */
    {
        EXPECT_THROW((bio::io::seq_io::reader{"/dev/nonexistant/foobarOOO", std::forward<decltype(args)>(args)...}),
                     bio::io::file_open_error);
    }

    /* wrong extension */
    if (ext_check)
    {
        bio::test::tmp_filename filename{"seq_io_reader_constructor.xyz"};
        std::ofstream           filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW((bio::io::seq_io::reader{filename.get_path(), std::forward<decltype(args)>(args)...}),
                     bio::io::unhandled_extension_error);
    }
}

TEST(seq_io_reader, constructor1_just_filename)
{
    seq_io_reader_filename_constructor(true);
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq_io::reader{""}), bio::io::seq_io::reader<>>));
}

TEST(seq_io_reader, constructor1_with_opts)
{
    bio::io::seq_io::reader_options opt{.field_types = bio::io::seq_io::field_types_protein};
    seq_io_reader_filename_constructor(true, std::move(opt));

    using control_t = bio::io::seq_io::reader<std::remove_cvref_t<decltype(bio::io::seq_io::default_field_ids)>,
                                              std::remove_cvref_t<decltype(bio::io::seq_io::field_types_protein)>,
                                              bio::meta::type_list<bio::io::fasta, bio::io::fastq>>;
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq_io::reader{"", opt}), control_t>));
}

TEST(seq_io_reader, constructor2_just_filename_direct_format)
{
    seq_io_reader_filename_constructor(false, bio::io::fasta{});
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq_io::reader{"", bio::io::fasta{}}), bio::io::seq_io::reader<>>));
}

TEST(seq_io_reader, constructor2_with_opts_direct_format)
{
    bio::io::seq_io::reader_options opt{.field_types = bio::io::seq_io::field_types_dna};
    seq_io_reader_filename_constructor(false, bio::io::fasta{}, std::move(opt));

    using control_t = bio::io::seq_io::reader<std::remove_cvref_t<decltype(bio::io::seq_io::default_field_ids)>,
                                              std::remove_cvref_t<decltype(bio::io::seq_io::field_types_dna)>,
                                              bio::meta::type_list<bio::io::fasta, bio::io::fastq>>;
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq_io::reader{"", bio::io::fasta{}, opt}), control_t>));
}

TEST(seq_io_reader, constructor2_just_filename_format_variant)
{
    bio::io::seq_io::reader<>::format_type var{};

    seq_io_reader_filename_constructor(false, var);
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq_io::reader{"", var}), bio::io::seq_io::reader<>>));
}

TEST(seq_io_reader, constructor2_with_opts_format_variant)
{
    bio::io::seq_io::reader<>::format_type var{};
    bio::io::seq_io::reader_options        opt{.field_types = bio::io::seq_io::field_types_dna};
    seq_io_reader_filename_constructor(false, var, std::move(opt));

    using control_t = bio::io::seq_io::reader<std::remove_cvref_t<decltype(bio::io::seq_io::default_field_ids)>,
                                              std::remove_cvref_t<decltype(bio::io::seq_io::field_types_dna)>,
                                              bio::meta::type_list<bio::io::fasta, bio::io::fastq>>;
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq_io::reader{"", var, std::move(opt)}), control_t>));
}

TEST(seq_io_reader, constructor3)
{
    std::istringstream str;

    EXPECT_NO_THROW((bio::io::seq_io::reader{str, bio::io::fasta{}}));
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq_io::reader{str, bio::io::fasta{}}), bio::io::seq_io::reader<>>));
}

TEST(seq_io_reader, constructor3_with_opts)
{
    std::istringstream              str;
    bio::io::seq_io::reader_options opt{.field_types = bio::io::seq_io::field_types_dna};
    EXPECT_NO_THROW((bio::io::seq_io::reader{str, bio::io::fasta{}, opt}));

    using control_t = bio::io::seq_io::reader<std::remove_cvref_t<decltype(bio::io::seq_io::default_field_ids)>,
                                              std::remove_cvref_t<decltype(bio::io::seq_io::field_types_dna)>,
                                              bio::meta::type_list<bio::io::fasta, bio::io::fastq>>;
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq_io::reader{str, bio::io::fasta{}, opt}), control_t>));
}

TEST(seq_io_reader, constructor4)
{
    std::istringstream str;

    EXPECT_NO_THROW((bio::io::seq_io::reader{std::move(str), bio::io::fasta{}}));
    EXPECT_TRUE(
      (std::same_as<decltype(bio::io::seq_io::reader{std::move(str), bio::io::fasta{}}), bio::io::seq_io::reader<>>));
}

TEST(seq_io_reader, constructor4_with_opts)
{
    std::istringstream              str;
    bio::io::seq_io::reader_options opt{.field_types = bio::io::seq_io::field_types_dna};
    EXPECT_NO_THROW((bio::io::seq_io::reader{std::move(str), bio::io::fasta{}, opt}));

    using control_t = bio::io::seq_io::reader<std::remove_cvref_t<decltype(bio::io::seq_io::default_field_ids)>,
                                              std::remove_cvref_t<decltype(bio::io::seq_io::field_types_dna)>,
                                              bio::meta::type_list<bio::io::fasta, bio::io::fastq>>;
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq_io::reader{std::move(str), bio::io::fasta{}, opt}), control_t>));
}

TEST(seq_io_reader, iteration)
{
    {
        std::istringstream      str{static_cast<std::string>(input)};
        bio::io::seq_io::reader reader{str, bio::io::fasta{}};

        EXPECT_EQ(std::ranges::distance(reader), 5);
    }

    {
        std::istringstream      str{static_cast<std::string>(input)};
        bio::io::seq_io::reader reader{str, bio::io::fasta{}};

        size_t count = 0;
        for (auto & rec : reader)
        {
            ++count;
            EXPECT_TRUE(rec.id().starts_with("ID"));
            // only very basic check here, rest in format test
        }
        EXPECT_EQ(count, 5ull);
    }
}

TEST(seq_io_reader, empty_file)
{
    {
        bio::test::tmp_filename filename{"seq_io_reader_constructor.fasta"};
        std::ofstream           filecreator{filename.get_path(), std::ios::out | std::ios::binary};

        bio::io::seq_io::reader reader{filename.get_path()};

        EXPECT_THROW(reader.begin(), bio::io::file_open_error);
    }
}

TEST(seq_io_reader, empty_stream)
{
    {
        std::istringstream      str{""};
        bio::io::seq_io::reader reader{str, bio::io::fasta{}};

        EXPECT_THROW(reader.begin(), bio::io::file_open_error);
    }
}

TEST(seq_io_reader, custom_field_types)
{
    bio::io::seq_io::reader_options opt{.field_types = bio::io::seq_io::field_types<bio::io::ownership::deep>};

    std::istringstream      str{static_cast<std::string>(input)};
    bio::io::seq_io::reader reader{str, bio::io::fasta{}, opt};

    EXPECT_TRUE((std::same_as<decltype(reader.front().seq()), std::vector<bio::alphabet::dna5> &>));
    EXPECT_TRUE((std::same_as<decltype(reader.front().id()), std::string &>));
}

TEST(seq_io_reader, custom_field_ids_structured_bindings)
{
    bio::io::seq_io::reader_options opt{.field_ids   = bio::meta::vtag<bio::io::field::seq, bio::io::field::id>,
                                        .field_types = bio::meta::ttag<std::string, std::string>};

    std::istringstream      str{static_cast<std::string>(input)};
    bio::io::seq_io::reader reader{str, bio::io::fasta{}, opt};

    for (auto & [seq, id] : reader)
        EXPECT_TRUE(id.starts_with("ID"));
}

TEST(seq_io_reader, decompression_filename)
{
    bio::test::tmp_filename filename{"seq_io_reader.fasta.gz"};

    {
        std::ofstream                             filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        bio::io::detail::fast_ostreambuf_iterator it{filecreator};
        it.write_range(input_bgzipped);
    }

    bio::io::seq_io::reader reader{filename.get_path()};

    size_t count = 0;
    for (auto & rec : reader)
    {
        ++count;
        EXPECT_TRUE(rec.id().starts_with("ID"));
        // only very basic check here, rest in format test
    }
    EXPECT_EQ(count, 5ull);
}

TEST(seq_io_reader, decompression_stream)
{
    std::istringstream str{static_cast<std::string>(input_bgzipped)};

    bio::io::seq_io::reader reader{str, bio::io::fasta{}};

    size_t count = 0;
    for (auto & rec : reader)
    {
        ++count;
        EXPECT_TRUE(rec.id().starts_with("ID"));
        // only very basic check here, rest in format test
    }
    EXPECT_EQ(count, 5ull);
}

// The following neads to cause a static assertion
// TEST(seq_io_reader, option_fail)
// {
//     bio::io::seq_io::reader_options opt{.field_types = bio::meta::ttag<int, int, int>};
// }
