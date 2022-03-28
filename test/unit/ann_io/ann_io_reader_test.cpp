// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <bio/ann_io/reader.hpp>

#include "../format/ann_data.hpp"

TEST(ann_io_reader, concepts)
{
    using t = bio::ann_io::reader<>;
    EXPECT_TRUE((std::ranges::input_range<t>));

    using ct = bio::ann_io::reader<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::input_range<ct>));
}

void ann_io_reader_filename_constructor(bool ext_check, auto &&... args)
{
    /* just the filename */
    {
        seqan3::test::tmp_filename filename{"ann_io_reader_constructor.bed"};
        std::ofstream              filecreator{filename.get_path(), std::ios::out | std::ios::binary};

        EXPECT_NO_THROW((bio::ann_io::reader{filename.get_path(), std::forward<decltype(args)>(args)...}));
    }

    // correct format check is done by tests of that format

    /* non-existent file */
    {
        EXPECT_THROW((bio::ann_io::reader{"/dev/nonexistant/foobarOOO", std::forward<decltype(args)>(args)...}),
                     bio::file_open_error);
    }

    /* wrong extension */
    if (ext_check)
    {
        seqan3::test::tmp_filename filename{"ann_io_reader_constructor.xyz"};
        std::ofstream              filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW((bio::ann_io::reader{filename.get_path(), std::forward<decltype(args)>(args)...}),
                     bio::unhandled_extension_error);
    }
}

TEST(ann_io_reader, constructor1_just_filename)
{
    ann_io_reader_filename_constructor(true);
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{""}), bio::ann_io::reader<>>));
}

TEST(ann_io_reader, constructor1_with_opts)
{
    bio::ann_io::reader_options opt{.field_types = bio::ann_io::field_types<>};
    using control_t = bio::ann_io::reader<decltype(bio::ann_io::default_field_ids),
                                          decltype(bio::ann_io::field_types<bio::ownership::shallow>),
                                          seqan3::type_list<bio::bed>>;

    ann_io_reader_filename_constructor(true, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{"", opt}), control_t>));
}

TEST(ann_io_reader, constructor2_just_filename_direct_format)
{
    ann_io_reader_filename_constructor(false, bio::bed{});
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{"", bio::bed{}}), bio::ann_io::reader<>>));
}

TEST(ann_io_reader, constructor2_with_opts_direct_format)
{
    bio::ann_io::reader_options opt{.field_types = bio::ann_io::field_types<>};
    using control_t = bio::ann_io::reader<decltype(bio::ann_io::default_field_ids),
                                          decltype(bio::ann_io::field_types<bio::ownership::shallow>),
                                          seqan3::type_list<bio::bed>>;

    ann_io_reader_filename_constructor(false, bio::bed{}, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{"", bio::bed{}, opt}), control_t>));
}

TEST(ann_io_reader, constructor2_just_filename_format_variant)
{
    std::variant<bio::bed> var{};

    ann_io_reader_filename_constructor(false, var);
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{"", var}), bio::ann_io::reader<>>));
}

TEST(ann_io_reader, constructor2_with_opts_format_variant)
{
    std::variant<bio::bed> var{};
    bio::ann_io::reader_options      opt{.field_types = bio::ann_io::field_types<>};
    using control_t = bio::ann_io::reader<decltype(bio::ann_io::default_field_ids),
                                          decltype(bio::ann_io::field_types<bio::ownership::shallow>),
                                          seqan3::type_list<bio::bed>>;

    ann_io_reader_filename_constructor(false, var, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{"", var, std::move(opt)}), control_t>));
}

TEST(ann_io_reader, constructor3)
{
    std::istringstream str;

    EXPECT_NO_THROW((bio::ann_io::reader{str, bio::bed{}}));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{str, bio::bed{}}), bio::ann_io::reader<>>));
}

TEST(ann_io_reader, constructor3_with_opts)
{
    std::istringstream          str;
    bio::ann_io::reader_options opt{.field_types = bio::ann_io::field_types<>};
    using control_t = bio::ann_io::reader<decltype(bio::ann_io::default_field_ids),
                                          decltype(bio::ann_io::field_types<bio::ownership::shallow>),
                                          seqan3::type_list<bio::bed>>;

    EXPECT_NO_THROW((bio::ann_io::reader{str, bio::bed{}, opt}));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{str, bio::bed{}, opt}), control_t>));
}

TEST(ann_io_reader, constructor4)
{
    std::istringstream str;

    EXPECT_NO_THROW((bio::ann_io::reader{std::move(str), bio::bed{}}));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{std::move(str), bio::bed{}}), bio::ann_io::reader<>>));
}

TEST(ann_io_reader, constructor4_with_opts)
{
    std::istringstream          str;
    bio::ann_io::reader_options opt{.field_types = bio::ann_io::field_types<>};
    using control_t = bio::ann_io::reader<decltype(bio::ann_io::default_field_ids),
                                          decltype(bio::ann_io::field_types<bio::ownership::shallow>),
                                          seqan3::type_list<bio::bed>>;

    EXPECT_NO_THROW((bio::ann_io::reader{std::move(str), bio::bed{}, opt}));
    EXPECT_TRUE((std::same_as<decltype(bio::ann_io::reader{std::move(str), bio::bed{}, opt}), control_t>));
}

TEST(ann_io_reader, iteration)
{
    {
        std::istringstream  str{static_cast<std::string>(minimal_example)};
        bio::ann_io::reader reader{str, bio::bed{}};

        EXPECT_EQ(std::ranges::distance(reader), 9);
    }

    {
        std::istringstream  str{static_cast<std::string>(minimal_example)};
        bio::ann_io::reader reader{str, bio::bed{}};

        size_t count = 0;
        for (auto & rec : reader)
        {
            ++count;
            EXPECT_EQ(rec.chrom(), "chr7");
            // only very basic check here, rest in format test
        }
        EXPECT_EQ(count, 9);
    }
}

TEST(ann_io_reader, empty_file)
{
    {
        seqan3::test::tmp_filename filename{"ann_io_reader_constructor.bed"};
        std::ofstream              filecreator{filename.get_path(), std::ios::out | std::ios::binary};

        bio::ann_io::reader reader{filename.get_path()};

        EXPECT_THROW(reader.begin(), bio::file_open_error);
    }
}

TEST(ann_io_reader, empty_stream)
{
    {
        std::istringstream  str{""};
        bio::ann_io::reader reader{str, bio::bed{}};

        EXPECT_THROW(reader.begin(), bio::file_open_error);
    }
}

TEST(ann_io_reader, get_header)
{
    // get header before calling begin()
    {
        std::istringstream  str{static_cast<std::string>(minimal_example_with_header)};
        bio::ann_io::reader reader{str, bio::bed{}};

        bio::ann_io::header const & hdr = reader.header();

        EXPECT_EQ(hdr.to_plaintext(), minimal_example_header_regenerated);
    }

    // get header after calling begin()
    {
        std::istringstream  str{static_cast<std::string>(minimal_example_with_header)};
        bio::ann_io::reader reader{str, bio::bed{}};

        auto it = reader.begin();
        EXPECT_EQ(it->chrom(), "chr7");

        bio::ann_io::header const & hdr = reader.header();

        EXPECT_EQ(hdr.to_plaintext(), minimal_example_header_regenerated);
    }
}

TEST(ann_io_reader, custom_field_types)
{
    bio::ann_io::reader_options opt{.field_types = bio::ann_io::field_types<bio::ownership::deep>};

    std::istringstream  str{static_cast<std::string>(minimal_example)};
    bio::ann_io::reader reader{str, bio::bed{}, opt};

    EXPECT_TRUE((std::same_as<std::ranges::range_value_t<decltype(reader)>,
                              bio::record<decltype(bio::ann_io::default_field_ids),
                                          decltype(bio::ann_io::field_types<bio::ownership::deep>)>>));
}

TEST(ann_io_reader, custom_field_ids_structured_bindings)
{
    bio::ann_io::reader_options opt{.field_ids   = bio::vtag<bio::field::chrom, bio::field::chromStart, bio::field::chromEnd>,
                                    .field_types = bio::ttag<std::string, uint32_t, uint32_t>};

    std::istringstream  str{static_cast<std::string>(minimal_example)};
    bio::ann_io::reader reader{str, bio::bed{}, opt};

    for (auto & [chrom, start, end] : reader)
        EXPECT_EQ(chrom, "chr7");
}
