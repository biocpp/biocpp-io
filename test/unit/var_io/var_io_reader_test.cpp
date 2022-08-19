// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/test/expect_range_eq.hpp>
#include <bio/test/tmp_directory.hpp>
#include <bio/test/tmp_filename.hpp>

#include <bio/io/var_io/reader.hpp>

#include "../format/vcf_data.hpp"

TEST(var_io_reader, concepts)
{
    using t = bio::io::var_io::reader<>;
    EXPECT_TRUE((std::ranges::input_range<t>));

    using ct = bio::io::var_io::reader<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::input_range<ct>));
}

void var_io_reader_filename_constructor(bool ext_check, auto &&... args)
{
    /* just the filename */
    {
        bio::test::tmp_filename filename{"var_io_reader_constructor.vcf"};
        std::ofstream           filecreator{filename.get_path(), std::ios::out | std::ios::binary};

        EXPECT_NO_THROW((bio::io::var_io::reader{filename.get_path(), std::forward<decltype(args)>(args)...}));
    }

    // correct format check is done by tests of that format

    /* non-existent file */
    {
        EXPECT_THROW((bio::io::var_io::reader{"/dev/nonexistant/foobarOOO", std::forward<decltype(args)>(args)...}),
                     bio::io::file_open_error);
    }

    /* wrong extension */
    if (ext_check)
    {
        bio::test::tmp_filename filename{"var_io_reader_constructor.xyz"};
        std::ofstream           filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW((bio::io::var_io::reader{filename.get_path(), std::forward<decltype(args)>(args)...}),
                     bio::io::unhandled_extension_error);
    }
}

TEST(var_io_reader, constructor1_just_filename)
{
    var_io_reader_filename_constructor(true);
    EXPECT_TRUE((std::same_as<decltype(bio::io::var_io::reader{""}), bio::io::var_io::reader<>>));
}

TEST(var_io_reader, constructor1_with_opts)
{
    bio::io::var_io::reader_options opt{.field_types = bio::io::var_io::field_types_bcf_style<>};
    using control_t =
      bio::io::var_io::reader<decltype(bio::io::var_io::default_field_ids),
                              decltype(bio::io::var_io::field_types_bcf_style<bio::io::ownership::shallow>),
                              bio::meta::type_list<bio::io::vcf, bio::io::bcf>>;

    var_io_reader_filename_constructor(true, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::io::var_io::reader{"", opt}), control_t>));
}

TEST(var_io_reader, constructor2_just_filename_direct_format)
{
    var_io_reader_filename_constructor(false, bio::io::vcf{});
    EXPECT_TRUE((std::same_as<decltype(bio::io::var_io::reader{"", bio::io::vcf{}}), bio::io::var_io::reader<>>));
}

TEST(var_io_reader, constructor2_with_opts_direct_format)
{
    bio::io::var_io::reader_options opt{.field_types = bio::io::var_io::field_types_bcf_style<>};
    using control_t =
      bio::io::var_io::reader<decltype(bio::io::var_io::default_field_ids),
                              decltype(bio::io::var_io::field_types_bcf_style<bio::io::ownership::shallow>),
                              bio::meta::type_list<bio::io::vcf, bio::io::bcf>>;

    var_io_reader_filename_constructor(false, bio::io::vcf{}, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::io::var_io::reader{"", bio::io::vcf{}, opt}), control_t>));
}

TEST(var_io_reader, constructor2_just_filename_format_variant)
{
    std::variant<bio::io::vcf, bio::io::bcf> var{};

    var_io_reader_filename_constructor(false, var);
    EXPECT_TRUE((std::same_as<decltype(bio::io::var_io::reader{"", var}), bio::io::var_io::reader<>>));
}

TEST(var_io_reader, constructor2_with_opts_format_variant)
{
    std::variant<bio::io::vcf, bio::io::bcf> var{};
    bio::io::var_io::reader_options          opt{.field_types = bio::io::var_io::field_types_bcf_style<>};
    using control_t =
      bio::io::var_io::reader<decltype(bio::io::var_io::default_field_ids),
                              decltype(bio::io::var_io::field_types_bcf_style<bio::io::ownership::shallow>),
                              bio::meta::type_list<bio::io::vcf, bio::io::bcf>>;

    var_io_reader_filename_constructor(false, var, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::io::var_io::reader{"", var, std::move(opt)}), control_t>));
}

TEST(var_io_reader, constructor3)
{
    std::istringstream str;

    EXPECT_NO_THROW((bio::io::var_io::reader{str, bio::io::vcf{}}));
    EXPECT_TRUE((std::same_as<decltype(bio::io::var_io::reader{str, bio::io::vcf{}}), bio::io::var_io::reader<>>));
}

TEST(var_io_reader, constructor3_with_opts)
{
    std::istringstream              str;
    bio::io::var_io::reader_options opt{.field_types = bio::io::var_io::field_types_bcf_style<>};
    using control_t =
      bio::io::var_io::reader<decltype(bio::io::var_io::default_field_ids),
                              decltype(bio::io::var_io::field_types_bcf_style<bio::io::ownership::shallow>),
                              bio::meta::type_list<bio::io::vcf, bio::io::bcf>>;

    EXPECT_NO_THROW((bio::io::var_io::reader{str, bio::io::vcf{}, opt}));
    EXPECT_TRUE((std::same_as<decltype(bio::io::var_io::reader{str, bio::io::vcf{}, opt}), control_t>));
}

TEST(var_io_reader, constructor4)
{
    std::istringstream str;

    EXPECT_NO_THROW((bio::io::var_io::reader{std::move(str), bio::io::vcf{}}));
    EXPECT_TRUE(
      (std::same_as<decltype(bio::io::var_io::reader{std::move(str), bio::io::vcf{}}), bio::io::var_io::reader<>>));
}

TEST(var_io_reader, constructor4_with_opts)
{
    std::istringstream              str;
    bio::io::var_io::reader_options opt{.field_types = bio::io::var_io::field_types_bcf_style<>};
    using control_t =
      bio::io::var_io::reader<decltype(bio::io::var_io::default_field_ids),
                              decltype(bio::io::var_io::field_types_bcf_style<bio::io::ownership::shallow>),
                              bio::meta::type_list<bio::io::vcf, bio::io::bcf>>;

    EXPECT_NO_THROW((bio::io::var_io::reader{std::move(str), bio::io::vcf{}, opt}));
    EXPECT_TRUE((std::same_as<decltype(bio::io::var_io::reader{std::move(str), bio::io::vcf{}, opt}), control_t>));
}

TEST(var_io_reader, iteration)
{
    {
        std::istringstream      str{static_cast<std::string>(example_from_spec)};
        bio::io::var_io::reader reader{str, bio::io::vcf{}};

        EXPECT_EQ(std::ranges::distance(reader), 5);
    }

    {
        std::istringstream      str{static_cast<std::string>(example_from_spec)};
        bio::io::var_io::reader reader{str, bio::io::vcf{}};

        size_t count = 0;
        for (auto & rec : reader)
        {
            ++count;
            EXPECT_EQ(rec.chrom(), "20");
            // only very basic check here, rest in format test
        }
        EXPECT_EQ(count, 5ull);
    }
}

TEST(var_io_reader, empty_file)
{
    {
        bio::test::tmp_filename filename{"var_io_reader_constructor.vcf"};
        std::ofstream           filecreator{filename.get_path(), std::ios::out | std::ios::binary};

        bio::io::var_io::reader reader{filename.get_path()};

        EXPECT_THROW(reader.begin(), bio::io::file_open_error);
    }
}

TEST(var_io_reader, empty_stream)
{
    {
        std::istringstream      str{""};
        bio::io::var_io::reader reader{str, bio::io::vcf{}};

        EXPECT_THROW(reader.begin(), bio::io::file_open_error);
    }
}

TEST(var_io_reader, get_header)
{
    // get header before calling begin()
    {
        std::istringstream      str{static_cast<std::string>(example_from_spec)};
        bio::io::var_io::reader reader{str, bio::io::vcf{}};

        bio::io::var_io::header const & hdr = reader.header();

        EXPECT_EQ(hdr.to_plaintext(), example_from_spec_header_regenerated);
    }

    // get header after calling begin()
    {
        std::istringstream      str{static_cast<std::string>(example_from_spec)};
        bio::io::var_io::reader reader{str, bio::io::vcf{}};

        auto it = reader.begin();
        EXPECT_EQ(it->chrom(), "20");

        bio::io::var_io::header const & hdr = reader.header();

        EXPECT_EQ(hdr.to_plaintext(), example_from_spec_header_regenerated);
    }
}

TEST(var_io_reader, custom_field_types)
{
    bio::io::var_io::reader_options opt{.field_types =
                                          bio::io::var_io::field_types_bcf_style<bio::io::ownership::deep>};

    std::istringstream      str{static_cast<std::string>(example_from_spec)};
    bio::io::var_io::reader reader{str, bio::io::vcf{}, opt};

    EXPECT_TRUE(
      (std::same_as<std::ranges::range_value_t<decltype(reader)>,
                    bio::io::record<decltype(bio::io::var_io::default_field_ids),
                                    decltype(bio::io::var_io::field_types_bcf_style<bio::io::ownership::deep>)>>));
}

TEST(var_io_reader, custom_field_ids_structured_bindings)
{
    bio::io::var_io::reader_options opt{
      .field_ids   = bio::meta::vtag<bio::io::field::chrom, bio::io::field::pos, bio::io::field::ref>,
      .field_types = bio::meta::ttag<std::string, uint32_t, std::string>};

    std::istringstream      str{static_cast<std::string>(example_from_spec)};
    bio::io::var_io::reader reader{str, bio::io::vcf{}, opt};

    for (auto & [chrom, pos, ref] : reader)
        EXPECT_EQ(chrom, "20");
}

TEST(var_io_reader, decompression_filename)
{
    bio::test::tmp_filename filename{"var_io_reader.vcf.gz"};

    {
        std::ofstream                             filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        bio::io::detail::fast_ostreambuf_iterator it{filecreator};
        it.write_range(example_from_spec_bgzipped);
    }

    bio::io::var_io::reader reader{filename.get_path()};

    size_t count = 0;
    for (auto & rec : reader)
    {
        ++count;
        EXPECT_EQ(rec.chrom(), "20");
        // only very basic check here, rest in format test
    }
    EXPECT_EQ(count, 5ull);
}

TEST(var_io_reader, decompression_stream)
{
    std::istringstream str{static_cast<std::string>(example_from_spec_bgzipped)};

    bio::io::var_io::reader reader{str, bio::io::vcf{}};

    size_t count = 0;
    for (auto & rec : reader)
    {
        ++count;
        EXPECT_EQ(rec.chrom(), "20");
        // only very basic check here, rest in format test
    }
    EXPECT_EQ(count, 5ull);
}

TEST(var_io_reader, region_filter)
{
    bio::io::genomic_region<>       region{.chrom = "20", .beg = 17000, .end = 1230300};
    bio::io::var_io::reader_options options{.region = region};

    bio::test::tmp_directory dir{};

    {
        std::ofstream os{dir.path() / "example.vcf.gz", std::ios::binary};
        os << example_from_spec_bgzipped;
    }

    {
        std::ofstream os{dir.path() / "example.vcf.gz.tbi", std::ios::binary};
        os << example_from_spec_bgzipped_tbi;
    }

    {
        bio::io::var_io::reader reader{dir.path() / "example.vcf.gz", options};

        size_t count = 0;
        for (auto & rec : reader)
        {
            ++count;
            EXPECT_EQ(rec.chrom(), "20");
            EXPECT_GE(rec.pos(), region.beg);
            EXPECT_LT(rec.pos(), region.end);
        }
        EXPECT_EQ(count, 3ull);
    }

    std::filesystem::remove(dir.path() / "example.vcf.gz");
    std::filesystem::remove(dir.path() / "example.vcf.gz.tbi");
}

// TODO region_filter_filename

TEST(var_io_reader, region_filter_linear)
{
    bio::io::genomic_region<>       region{.chrom = "20", .beg = 17000, .end = 1230300};
    bio::io::var_io::reader_options options{.region = region, .region_index_optional = true};

    {
        std::istringstream      str{static_cast<std::string>(example_from_spec)};
        bio::io::var_io::reader reader{str, bio::io::vcf{}, options};

        EXPECT_EQ(std::ranges::distance(reader), 3);
    }

    {
        std::istringstream      str{static_cast<std::string>(example_from_spec)};
        bio::io::var_io::reader reader{str, bio::io::vcf{}, options};

        size_t count = 0;
        for (auto & rec : reader)
        {
            ++count;
            EXPECT_EQ(rec.chrom(), "20");
            EXPECT_GE(rec.pos(), region.beg);
            EXPECT_LT(rec.pos(), region.end);
        }
        EXPECT_EQ(count, 3ull);
    }
}

TEST(var_io_reader, reopen)
{
    bio::test::tmp_directory dir{};

    {
        std::ofstream os{dir.path() / "example.vcf.gz", std::ios::binary};
        os << example_from_spec_bgzipped;
    }

    {
        bio::io::var_io::reader reader{dir.path() / "example.vcf.gz"};

        size_t count = 0;
        for (auto & rec : reader)
        {
            ++count;
            EXPECT_EQ(rec.chrom(), "20");
        }
        EXPECT_EQ(count, 5ull);

        reader.reopen();

        count = 0;
        for (auto & rec : reader)
        {
            ++count;
            EXPECT_EQ(rec.chrom(), "20");
        }
        EXPECT_EQ(count, 5ull);
    }

    std::filesystem::remove(dir.path() / "example.vcf.gz");
}

// TODO tests for reopen(region) on 1000G_chr10_sample.vcf.gz
