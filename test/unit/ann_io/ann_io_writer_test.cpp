// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <bio/stream/transparent_istream.hpp>
#include <bio/ann_io/writer.hpp>

#include "../format/ann_data.hpp"

using custom_field_ids_t = bio::vtag_t<bio::field::chrom, bio::field::chromStart, bio::field::chromEnd>;

TEST(ann_io_writer, concepts)
{
    using rec_t = bio::record<decltype(bio::ann_io::default_field_ids), decltype(bio::ann_io::field_types<>)>;

    using t = bio::ann_io::writer<>;
    EXPECT_TRUE((std::ranges::output_range<t, rec_t>));

    using ct = bio::ann_io::writer<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::output_range<ct, rec_t>));
}

void ann_io_writer_filename_constructor(bool ext_check, auto &&... args)
{
    using t =
      decltype(bio::ann_io::writer{std::declval<std::filesystem::path &>(), std::forward<decltype(args)>(args)...});
    [[maybe_unused]] t * ptr = nullptr;

    /* just the filename */
    {
        seqan3::test::tmp_filename filename{"ann_io_writer_constructor.bed"};

        // constructor
        EXPECT_NO_THROW((ptr = new bio::ann_io::writer{filename.get_path(), std::forward<decltype(args)>(args)...}));

        // destructor
        EXPECT_NO_THROW(delete ptr);
        ptr = nullptr;
    }

    /* wrong extension */
    if (ext_check)
    {
        seqan3::test::tmp_filename filename{"ann_io_writer_constructor.xyz"};
        EXPECT_THROW((ptr = new bio::ann_io::writer{filename.get_path(), std::forward<decltype(args)>(args)...}),
                     bio::unhandled_extension_error);

        // destructor, nothrow because already thrown during construction
        EXPECT_NO_THROW(delete ptr);
    }
}

TEST(ann_io_writer, constructor1_just_filename)
{
    ann_io_writer_filename_constructor(true);
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{""}), bio::ann_io::writer<>>));
}

TEST(ann_io_writer, constructor1_with_opts)
{
    bio::ann_io::writer_options opt{.formats = bio::ttag<bio::bed>};
    using control_t = bio::ann_io::writer<seqan3::type_list<bio::bed>>;

    ann_io_writer_filename_constructor(true, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{"", opt}), control_t>));
}

TEST(ann_io_writer, constructor2_just_filename_direct_format)
{
    ann_io_writer_filename_constructor(false, bio::bed{});
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{"", bio::bed{}}), bio::ann_io::writer<>>));
}

TEST(ann_io_writer, constructor2_with_opts_direct_format)
{
    bio::ann_io::writer_options opt{.formats = bio::ttag<bio::bed>};
    using control_t = bio::ann_io::writer<seqan3::type_list<bio::bed>>;

    ann_io_writer_filename_constructor(false, bio::bed{}, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{"", bio::bed{}, opt}), control_t>));
}

TEST(ann_io_writer, constructor2_just_filename_format_variant)
{
    std::variant<bio::bed> var{};

    ann_io_writer_filename_constructor(false, var);
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{"", var}), bio::ann_io::writer<>>));
}

TEST(ann_io_writer, constructor2_with_opts_format_variant)
{
    std::variant<bio::bed>      var{};
    bio::ann_io::writer_options opt{.formats = bio::ttag<bio::bed>};
    using control_t = bio::ann_io::writer<seqan3::type_list<bio::bed>>;

    ann_io_writer_filename_constructor(false, var, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{"", var, std::move(opt)}), control_t>));
}

void ann_io_writer_stream_constructor(auto &&... args)
{
    using t                  = decltype(bio::ann_io::writer{std::forward<decltype(args)>(args)...});
    [[maybe_unused]] t * ptr = nullptr;

    {
        std::ostringstream str;

        // constructor
        EXPECT_NO_THROW((ptr = new bio::ann_io::writer{std::forward<decltype(args)>(args)...}));

        // destructor
        EXPECT_NO_THROW(delete ptr);
        ptr = nullptr;
    }
}

TEST(ann_io_writer, constructor3)
{
    std::ostringstream str;
    ann_io_writer_stream_constructor(str, bio::bed{});

    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{str, bio::bed{}}), bio::ann_io::writer<>>));
}

TEST(ann_io_writer, constructor3_with_opts)
{
    std::ostringstream          str;
    bio::ann_io::writer_options opt{.formats = bio::ttag<bio::bed>};
    ann_io_writer_stream_constructor(str, bio::bed{}, opt);

    using control_t = bio::ann_io::writer<seqan3::type_list<bio::bed>>;
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{str, bio::bed{}, opt}), control_t>));
}

TEST(ann_io_writer, constructor4)
{
    std::ostringstream str;
    ann_io_writer_stream_constructor(std::move(str), bio::bed{});

    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{std::move(str), bio::bed{}}), bio::ann_io::writer<>>));
}

TEST(ann_io_writer, constructor4_with_opts)
{
    std::ostringstream          str;
    bio::ann_io::writer_options opt{.formats = bio::ttag<bio::bed>};
    ann_io_writer_stream_constructor(std::move(str), bio::bed{}, opt);

    using control_t = bio::ann_io::writer<seqan3::type_list<bio::bed>>;
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::writer{std::move(str), bio::bed{}, opt}), control_t>));
}

template <size_t i>
void write_record_test_impl()
{
    bio::ann_io::header hdr{minimal_example_header_regenerated};

    std::ostringstream  stream{};
    bio::ann_io::writer writer{stream, bio::bed{}};

    auto records = example_records_default_style<bio::ownership::shallow>();

    writer.set_header(hdr);

    if constexpr (i == 0)
    {
        writer.push_back(records[0]);
        writer.push_back(records[1]);
        writer.push_back(records[2]);
        writer.push_back(records[3]);
        writer.push_back(records[4]);
        writer.push_back(records[5]);
        writer.push_back(records[6]);
        writer.push_back(records[7]);
        writer.push_back(records[8]);
    }
    else if constexpr (i == 1)
    {
        auto it = writer.begin();
        it      = records[0];
        it      = records[1];
        it      = records[2];
        it      = records[3];
        it      = records[4];
        it      = records[5];
        it      = records[6];
        it      = records[7];
        it      = records[8];
    }
    else if constexpr (i == 2)
    {
        auto it = writer.begin();
        *it     = records[0];
        *it     = records[1];
        *it     = records[2];
        *it     = records[3];
        *it     = records[4];
        *it     = records[5];
        *it     = records[6];
        *it     = records[7];
        *it     = records[8];
    }
    else if constexpr (i == 4)
    {
        auto fn = [&writer](auto &... args) { writer.emplace_back(args...); };
        std::apply(fn, records[0]);
        std::apply(fn, records[1]);
        std::apply(fn, records[2]);
        std::apply(fn, records[3]);
        std::apply(fn, records[4]);
        std::apply(fn, records[5]);
        std::apply(fn, records[6]);
        std::apply(fn, records[7]);
        std::apply(fn, records[8]);
    }

    EXPECT_EQ(stream.str(), minimal_example_header_regenerated + minimal_example);
}

TEST(ann_io_writer, push_back_record)
{
    write_record_test_impl<0>();
}

TEST(ann_io_writer, assign_to_iterator)
{
    write_record_test_impl<1>();
}

TEST(ann_io_writer, assign_to_deref_iterator)
{
    write_record_test_impl<2>();
}

TEST(ann_io_writer, emplace_back)
{
    write_record_test_impl<4>();
}

// TEST(ann_io_writer, minimal_fields)
// {
//     std::ostringstream  stream{};
//     bio::ann_io::writer writer{stream, bio::bed{}};
//
//     writer.set_header(bio::ann_io::header{minimal_example_header_regenerated});
//
//     writer.emplace_back(custom_field_ids_t{}, "20", 14370, "G");
//     writer.emplace_back(custom_field_ids_t{}, "20", 17330, "T");
//     writer.emplace_back(custom_field_ids_t{}, "20", 1110696, "A");
//     writer.emplace_back(custom_field_ids_t{}, "20", 1230237, "T");
//     writer.emplace_back(custom_field_ids_t{}, "20", 1234567, "GTC");
//
//     std::string compare = example_from_spec_header_regenerated_no_IDX;
//     compare += minimal_field_rows;
//     EXPECT_EQ(stream.str(), compare);
// }

TEST(ann_io_writer, no_header1) // record contains header_ptr but this is == nullptr
{
    std::ostringstream stream{};
    auto *             writer = new bio::ann_io::writer{stream, bio::bed{}};

    auto records = example_records_default_style<bio::ownership::shallow>();

    EXPECT_NO_THROW(writer->push_back(records[0]));

    // destructor
    EXPECT_NO_THROW(delete writer);
    writer = nullptr;
}

TEST(var_io_writer, no_header2) // record does not contain header_ptr
{
    std::ostringstream stream{};
    auto *             writer = new bio::ann_io::writer{stream, bio::bed{}};

    EXPECT_NO_THROW(writer->emplace_back(custom_field_ids_t{}, 20, 1, 5));

    // destructor
    EXPECT_NO_THROW(delete writer);
    writer = nullptr;
}

// TEST(var_io_writer, compression)
// {
//     std::ostringstream stream{};
//
//     {
//         bio::var_io::writer writer{stream,
//                                    bio::vcf{},
//                                    bio::var_io::writer_options{.stream_options = bio::transparent_ostream_options{
//                                                                  .compression = bio::compression_format::bgzf}}};
//
//         writer.set_header(bio::var_io::header{example_from_spec_header});
//
//         auto records = example_records_bcf_style<bio::ownership::shallow>();
//
//         writer.push_back(records[0]);
//         writer.push_back(records[1]);
//         writer.push_back(records[2]);
//         writer.push_back(records[3]);
//         writer.push_back(records[4]);
//     }
//
//     std::string str = stream.str();
//     EXPECT_TRUE(str.starts_with("\x1f\x8b\x08")); // Gzip header
//
//     std::istringstream       control_stream{str};
//     bio::transparent_istream decompressor{control_stream};
//     std::string              buffer(std::istreambuf_iterator<char>{decompressor}, std::istreambuf_iterator<char>{});
//     EXPECT_RANGE_EQ(buffer, example_from_spec_header_regenerated_no_IDX + example_from_spec_records);
// }
