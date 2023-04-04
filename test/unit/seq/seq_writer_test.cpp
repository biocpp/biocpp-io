// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>
#include <bio/test/tmp_filename.hpp>

#include <bio/io/format/fasta.hpp>
#include <bio/io/seq/writer.hpp>
#include <bio/io/stream/transparent_istream.hpp>

#include "../format/seq_output_detail.hpp"

//TODO these should go somewhere central instead
inline std::string_view fasta_default_output = R"(>ID1
ACGTTTTTTTTTTTTTTT

>ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTT

>ID3 lala
ACGTTTA
)";

TEST(seq_writer, concepts)
{
    using rec_t = bio::io::seq::record_dna_deep;

    using t = bio::io::seq::writer<>;
    EXPECT_TRUE((std::ranges::output_range<t, rec_t>));

    using ct = bio::io::seq::writer<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::output_range<ct, rec_t>));
}

void seq_writer_filename_constructor(bool ext_check, auto &&... args)
{
    using t =
      decltype(bio::io::seq::writer{std::declval<std::filesystem::path &>(), std::forward<decltype(args)>(args)...});
    [[maybe_unused]] t * ptr = nullptr;

    /* just the filename */
    {
        bio::test::tmp_filename filename{"seq_writer_constructor.fasta"};

        // constructor
        EXPECT_NO_THROW((ptr = new bio::io::seq::writer{filename.get_path(), std::forward<decltype(args)>(args)...}));

        // destructor
        EXPECT_NO_THROW(delete ptr);
        ptr = nullptr;
    }

    /* wrong extension */
    if (ext_check)
    {
        bio::test::tmp_filename filename{"seq_writer_constructor.xyz"};
        EXPECT_THROW((ptr = new bio::io::seq::writer{filename.get_path(), std::forward<decltype(args)>(args)...}),
                     bio::io::unhandled_extension_error);

        // destructor, nothrow because already thrown during construction
        EXPECT_NO_THROW(delete ptr);
    }
}

TEST(seq_writer, constructor1_just_filename)
{
    seq_writer_filename_constructor(true);
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq::writer{""}), bio::io::seq::writer<>>));
}

TEST(seq_writer, constructor1_with_opts)
{
    bio::io::seq::writer_options opt{.formats = bio::meta::ttag<bio::io::fasta>};
    using control_t = bio::io::seq::writer<bio::meta::type_list<bio::io::fasta>>;

    seq_writer_filename_constructor(true, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq::writer{"", opt}), control_t>));
}

TEST(seq_writer, constructor2_just_filename_direct_format)
{
    seq_writer_filename_constructor(false, bio::io::fasta{});
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq::writer{"", bio::io::fasta{}}), bio::io::seq::writer<>>));
}

TEST(seq_writer, constructor2_with_opts_direct_format)
{
    bio::io::seq::writer_options opt{.formats = bio::meta::ttag<bio::io::fasta>};
    using control_t = bio::io::seq::writer<bio::meta::type_list<bio::io::fasta>>;

    seq_writer_filename_constructor(false, bio::io::fasta{}, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq::writer{"", bio::io::fasta{}, opt}), control_t>));
}

TEST(seq_writer, constructor2_just_filename_format_variant)
{
    std::variant<bio::io::fasta, bio::io::fastq> var{};

    seq_writer_filename_constructor(false, var);
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq::writer{"", var}), bio::io::seq::writer<>>));
}

TEST(seq_writer, constructor2_with_opts_format_variant)
{
    std::variant<bio::io::fasta> var{};
    bio::io::seq::writer_options opt{.formats = bio::meta::ttag<bio::io::fasta>};
    using control_t = bio::io::seq::writer<bio::meta::type_list<bio::io::fasta>>;

    seq_writer_filename_constructor(false, var, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq::writer{"", var, std::move(opt)}), control_t>));
}

void seq_writer_stream_constructor(auto &&... args)
{
    using t                  = decltype(bio::io::seq::writer{std::forward<decltype(args)>(args)...});
    [[maybe_unused]] t * ptr = nullptr;

    {
        std::ostringstream str;

        // constructor
        EXPECT_NO_THROW((ptr = new bio::io::seq::writer{std::forward<decltype(args)>(args)...}));

        // destructor
        EXPECT_NO_THROW(delete ptr);
        ptr = nullptr;
    }
}

TEST(seq_writer, constructor3)
{
    std::ostringstream str;
    seq_writer_stream_constructor(str, bio::io::fasta{});

    EXPECT_TRUE((std::same_as<decltype(bio::io::seq::writer{str, bio::io::fasta{}}), bio::io::seq::writer<>>));
}

TEST(seq_writer, constructor3_with_opts)
{
    std::ostringstream           str;
    bio::io::seq::writer_options opt{.formats = bio::meta::ttag<bio::io::fasta>};
    seq_writer_stream_constructor(str, bio::io::fasta{}, opt);

    using control_t = bio::io::seq::writer<bio::meta::type_list<bio::io::fasta>>;
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq::writer{str, bio::io::fasta{}, opt}), control_t>));
}

TEST(seq_writer, constructor4)
{
    std::ostringstream str;
    seq_writer_stream_constructor(std::move(str), bio::io::fasta{});

    EXPECT_TRUE(
      (std::same_as<decltype(bio::io::seq::writer{std::move(str), bio::io::fasta{}}), bio::io::seq::writer<>>));
}

TEST(seq_writer, constructor4_with_opts)
{
    std::ostringstream           str;
    bio::io::seq::writer_options opt{.formats = bio::meta::ttag<bio::io::fasta>};
    seq_writer_stream_constructor(std::move(str), bio::io::fasta{}, opt);

    using control_t = bio::io::seq::writer<bio::meta::type_list<bio::io::fasta>>;
    EXPECT_TRUE((std::same_as<decltype(bio::io::seq::writer{std::move(str), bio::io::fasta{}, opt}), control_t>));
}

template <size_t i>
void write_record_test_impl()
{
    std::ostringstream   stream{};
    bio::io::seq::writer writer{stream, bio::io::fasta{}};

    auto records = example_records<false, char, char>();

    if constexpr (i == 0)
    {
        writer.push_back(records[0]);
        writer.push_back(records[1]);
        writer.push_back(records[2]);
    }
    else if constexpr (i == 1)
    {
        auto it = writer.begin();
        it      = records[0];
        it      = records[1];
        it      = records[2];
    }
    else if constexpr (i == 2)
    {
        auto it = writer.begin();
        *it     = records[0];
        *it     = records[1];
        *it     = records[2];
    }
    else if constexpr (i == 3)
    {
        auto fn = [&writer](auto & r) { writer.emplace_back(r.id, r.seq, r.qual); };
        fn(records[0]);
        fn(records[1]);
        fn(records[2]);
    }

    EXPECT_EQ(stream.str(), fasta_default_output);
}

TEST(seq_writer, push_back_record)
{
    write_record_test_impl<0>();
}

TEST(seq_writer, assign_to_iterator)
{
    write_record_test_impl<1>();
}

TEST(seq_writer, assign_to_deref_iterator)
{
    write_record_test_impl<2>();
}

TEST(seq_writer, emplace_back)
{
    write_record_test_impl<3>();
}

TEST(seq_writer, compression)
{
    std::ostringstream stream{};

    {
        bio::io::seq::writer writer{stream,
                                    bio::io::fasta{},
                                    bio::io::seq::writer_options{.stream_options = bio::io::transparent_ostream_options{
                                                                   .compression = bio::io::compression_format::bgzf}}};

        auto records = example_records<false, char, char>();

        writer.push_back(records[0]);
        writer.push_back(records[1]);
        writer.push_back(records[2]);
    }

    std::string str = stream.str();
    EXPECT_TRUE(str.starts_with("\x1f\x8b\x08")); // Gzip header

    std::istringstream           control_stream{str};
    bio::io::transparent_istream decompressor{control_stream};
    std::string                  buffer(std::istreambuf_iterator<char>{decompressor}, std::istreambuf_iterator<char>{});
    EXPECT_RANGE_EQ(buffer, fasta_default_output);
}
