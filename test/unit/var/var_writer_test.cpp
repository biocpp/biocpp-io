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
#include <bio/io/stream/transparent_istream.hpp>
#include <bio/io/var/writer.hpp>

#include "../format/vcf_data.hpp"

TEST(var_writer, concepts)
{
    using rec_t = bio::io::var::record_idx_deep;

    using t = bio::io::var::writer<>;
    EXPECT_TRUE((std::ranges::output_range<t, rec_t>));

    using ct = bio::io::var::writer<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::output_range<ct, rec_t>));
}

void var_writer_filename_constructor(bool ext_check, auto &&... args)
{
    using t =
      decltype(bio::io::var::writer{std::declval<std::filesystem::path &>(), std::forward<decltype(args)>(args)...});
    [[maybe_unused]] t * ptr = nullptr;

    /* just the filename */
    {
        bio::test::tmp_filename filename{"var_writer_constructor.vcf"};

        // constructor
        EXPECT_NO_THROW((ptr = new bio::io::var::writer{filename.get_path(), std::forward<decltype(args)>(args)...}));

        // destructor
        EXPECT_THROW(delete ptr, bio::io::missing_header_error);
        ptr = nullptr;
    }

    /* wrong extension */
    if (ext_check)
    {
        bio::test::tmp_filename filename{"var_writer_constructor.xyz"};
        EXPECT_THROW((ptr = new bio::io::var::writer{filename.get_path(), std::forward<decltype(args)>(args)...}),
                     bio::io::unhandled_extension_error);

        // destructor, nothrow because already thrown during construction
        EXPECT_NO_THROW(delete ptr);
    }
}

TEST(var_writer, constructor1_just_filename)
{
    var_writer_filename_constructor(true);
    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{""}), bio::io::var::writer<>>));
}

TEST(var_writer, constructor1_with_opts)
{
    bio::io::var::writer_options opt{.formats = bio::meta::ttag<bio::io::vcf>};
    using control_t = bio::io::var::writer<bio::meta::type_list<bio::io::vcf>>;

    var_writer_filename_constructor(true, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{"", opt}), control_t>));
}

TEST(var_writer, constructor2_just_filename_direct_format)
{
    var_writer_filename_constructor(false, bio::io::vcf{});
    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{"", bio::io::vcf{}}), bio::io::var::writer<>>));
}

TEST(var_writer, constructor2_with_opts_direct_format)
{
    bio::io::var::writer_options opt{.formats = bio::meta::ttag<bio::io::vcf>};
    using control_t = bio::io::var::writer<bio::meta::type_list<bio::io::vcf>>;

    var_writer_filename_constructor(false, bio::io::vcf{}, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{"", bio::io::vcf{}, opt}), control_t>));
}

TEST(var_writer, constructor2_just_filename_format_variant)
{
    std::variant<bio::io::bcf, bio::io::vcf> var{};

    var_writer_filename_constructor(false, var);
    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{"", var}), bio::io::var::writer<>>));
}

TEST(var_writer, constructor2_with_opts_format_variant)
{
    std::variant<bio::io::vcf>   var{};
    bio::io::var::writer_options opt{.formats = bio::meta::ttag<bio::io::vcf>};
    using control_t = bio::io::var::writer<bio::meta::type_list<bio::io::vcf>>;

    var_writer_filename_constructor(false, var, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{"", var, std::move(opt)}), control_t>));
}

void var_writer_stream_constructor(auto &&... args)
{
    using t                  = decltype(bio::io::var::writer{std::forward<decltype(args)>(args)...});
    [[maybe_unused]] t * ptr = nullptr;

    {
        std::ostringstream str;

        // constructor
        EXPECT_NO_THROW((ptr = new bio::io::var::writer{std::forward<decltype(args)>(args)...}));

        // destructor
        EXPECT_THROW(delete ptr, bio::io::missing_header_error);
        ptr = nullptr;
    }
}

TEST(var_writer, constructor3)
{
    std::ostringstream str;
    var_writer_stream_constructor(str, bio::io::vcf{});

    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{str, bio::io::vcf{}}), bio::io::var::writer<>>));
}

TEST(var_writer, constructor3_with_opts)
{
    std::ostringstream           str;
    bio::io::var::writer_options opt{.formats = bio::meta::ttag<bio::io::vcf>};
    var_writer_stream_constructor(str, bio::io::vcf{}, opt);

    using control_t = bio::io::var::writer<bio::meta::type_list<bio::io::vcf>>;
    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{str, bio::io::vcf{}, opt}), control_t>));
}

TEST(var_writer, constructor4)
{
    std::ostringstream str;
    var_writer_stream_constructor(std::move(str), bio::io::vcf{});

    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{std::move(str), bio::io::vcf{}}), bio::io::var::writer<>>));
}

TEST(var_writer, constructor4_with_opts)
{
    std::ostringstream           str;
    bio::io::var::writer_options opt{.formats = bio::meta::ttag<bio::io::vcf>};
    var_writer_stream_constructor(std::move(str), bio::io::vcf{}, opt);

    using control_t = bio::io::var::writer<bio::meta::type_list<bio::io::vcf>>;
    EXPECT_TRUE((std::same_as<decltype(bio::io::var::writer{std::move(str), bio::io::vcf{}, opt}), control_t>));
}

template <size_t i>
void write_record_test_impl()
{
    bio::io::var::header hdr{example_from_spec_header};

    bio::io::var::record_private_data priv{};

    std::ostringstream   stream{};
    bio::io::var::writer writer{stream, bio::io::vcf{}};

    auto records = example_records_bcf_style<bio::io::ownership::shallow>();

    if constexpr (i == 3)
    {
        priv.header_ptr = &hdr;
        std::apply([&](auto &... rec) { ((rec._private = priv), ...); }, records);
    }
    else
    {
        writer.set_header(hdr);
    }

    if constexpr (i == 0 || i == 3)
    {
        std::apply([&](auto &... rec) { (writer.push_back(rec), ...); }, records);
    }
    else if constexpr (i == 1)
    {
        auto it = writer.begin();
        std::apply([&](auto &... rec) { ((it = rec), ...); }, records);
    }
    else if constexpr (i == 2)
    {
        auto it = writer.begin();
        std::apply([&](auto &... rec) { ((*it = rec), ...); }, records);
    }
    else if constexpr (i == 4)
    {
        std::apply(
          [&](auto &... r)
          { (writer.emplace_back(r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, r.info, r.genotypes), ...); },
          records);
    }

    EXPECT_EQ(stream.str(), example_from_spec_header_regenerated_no_IDX + example_from_spec_records);
}

TEST(var_writer, push_back_record)
{
    write_record_test_impl<0>();
}

TEST(var_writer, assign_to_iterator)
{
    write_record_test_impl<1>();
}

TEST(var_writer, assign_to_deref_iterator)
{
    write_record_test_impl<2>();
}

TEST(var_writer, implicit_header)
{
    write_record_test_impl<3>();
}

TEST(var_writer, emplace_back)
{
    write_record_test_impl<4>();
}

TEST(var_writer, minimal_fields)
{
    std::ostringstream   stream{};
    bio::io::var::writer writer{stream, bio::io::vcf{}};

    writer.set_header(bio::io::var::header{example_from_spec_header});

    bio::io::var::record r{.chrom     = {},
                           .pos       = {},
                           .id        = std::ignore,
                           .ref       = std::string{},
                           .alt       = std::ignore,
                           .qual      = std::ignore,
                           .filter    = std::ignore,
                           .info      = std::ignore,
                           .genotypes = std::ignore};
    r.chrom = "20";
    r.pos   = 14370;
    r.ref   = "G";
    writer.push_back(r);
    r.chrom = "20";
    r.pos   = 17330;
    r.ref   = "T";
    writer.push_back(r);
    r.chrom = "20";
    r.pos   = 1110696;
    r.ref   = "A";
    writer.push_back(r);
    r.chrom = "20";
    r.pos   = 1230237;
    r.ref   = "T";
    writer.push_back(r);
    r.chrom = "20";
    r.pos   = 1234567;
    r.ref   = "GTC";
    writer.push_back(r);

    std::string compare = example_from_spec_header_regenerated_no_IDX;
    compare += minimal_field_rows;
    EXPECT_EQ(stream.str(), compare);
}

TEST(var_writer, no_header1) // record contains header_ptr but this is == nullptr
{
    std::ostringstream stream{};
    auto *             writer = new bio::io::var::writer{stream, bio::io::vcf{}};

    auto records = example_records_bcf_style<bio::io::ownership::shallow>();

    EXPECT_THROW(writer->push_back(std::get<0>(records)), bio::io::missing_header_error);

    // destructor
    EXPECT_THROW(delete writer, bio::io::missing_header_error);
    writer = nullptr;
}

TEST(var_writer, no_header2) // record does not contain header_ptr
{
    std::ostringstream stream{};
    auto *             writer  = new bio::io::var::writer{stream, bio::io::vcf{}};
    auto               records = example_records_bcf_style<bio::io::ownership::shallow>();
    auto &             r       = std::get<0>(records);

    EXPECT_THROW((writer->emplace_back(r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, r.info, r.genotypes)),
                 bio::io::missing_header_error);

    // destructor
    EXPECT_THROW(delete writer, bio::io::missing_header_error);
    writer = nullptr;
}

TEST(var_writer, compression)
{
    std::ostringstream stream{};

    {
        bio::io::var::writer writer{stream,
                                    bio::io::vcf{},
                                    bio::io::var::writer_options{.stream_options = bio::io::transparent_ostream_options{
                                                                   .compression = bio::io::compression_format::bgzf}}};

        writer.set_header(bio::io::var::header{example_from_spec_header});

        auto records = example_records_bcf_style<bio::io::ownership::shallow>();

        std::apply([&](auto &... rec) { (writer.push_back(rec), ...); }, records);
    }

    std::string str = stream.str();
    EXPECT_TRUE(str.starts_with("\x1f\x8b\x08")); // Gzip header

    std::istringstream           control_stream{str};
    bio::io::transparent_istream decompressor{control_stream};
    std::string                  buffer(std::istreambuf_iterator<char>{decompressor}, std::istreambuf_iterator<char>{});
    EXPECT_RANGE_EQ(buffer, example_from_spec_header_regenerated_no_IDX + example_from_spec_records);
}

TEST(var_writer, biocpp_io_issue_53)
{
    using namespace bio::alphabet::literals;

    {
        using writer_t = decltype(bio::io::var::writer{std::cout, bio::io::vcf{}});
        writer_t * ptr = nullptr;

        EXPECT_THROW((ptr = new writer_t{std::cout, bio::io::vcf{}}), bio::io::sync_with_stdio_detected);
        delete ptr;
        ptr = nullptr;
    }

    std::ios::sync_with_stdio(false);

    {
        bio::io::var::header hdr{};
        hdr.file_format   = "VCFv4.3";
        hdr.column_labels = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "x"};

        bio::io::var::writer writer{std::cout, bio::io::vcf{}};
        writer.set_header(hdr);

        bio::io::var::record record{};
        record.chrom     = "chr1";
        record.pos       = 11111;
        record.id        = "test";
        record.ref       = "ATC"_dna5;
        record.alt       = {"AGC", "A"};
        record.qual      = 1.F;
        record.filter    = {"PASS"};
        record.genotypes = {};
        record.info      = {};
        writer.push_back(record);
    }
}

TEST(var_writer, get_header)
{
    {
        bio::io::var::writer writer{"var_writer_constructor.vcf"};
        bio::io::var::header hdr;
        writer.set_header(std::move(hdr));
        EXPECT_NO_THROW(writer.header());
    }

    {
        bio::io::var::writer writer{"var_writer_constructor.vcf"};
        bio::io::var::header hdr;
        writer.set_header(std::move(hdr));
        {
            bio::io::var::writer writer2{std::move(writer)};
        }
        EXPECT_THROW(writer.header(), bio::io::bio_error); // in moved-from state
    }
}
