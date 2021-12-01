// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/detail/debug_stream_byte.hpp>
#include <seqan3/io/stream/istreambuf_view.hpp>
#include <seqan3/io/stream/transparent_istream.hpp>
#include <seqan3/io/variant_io/writer.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include "../format/vcf_data.hpp"

using custom_field_ids_t = seqan3::tag_t<seqan3::field::chrom, seqan3::field::pos, seqan3::field::ref>;

TEST(var_io_writer, concepts)
{
    using rec_t = seqan3::record<seqan3::list_traits::repeat<seqan3::var_io::default_field_ids.size, std::string>,
                                 std::remove_cvref_t<decltype(seqan3::var_io::default_field_ids)>>;

    using t = seqan3::var_io::writer<>;
    EXPECT_TRUE((std::ranges::output_range<t, rec_t>));

    using ct = seqan3::var_io::writer<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::output_range<ct, rec_t>));
}

void var_io_writer_filename_constructor(bool ext_check, auto &&... args)
{
    /* just the filename */
    {
        seqan3::test::tmp_filename filename{"var_io_writer_constructor.vcf"};
        EXPECT_NO_THROW((seqan3::var_io::writer{filename.get_path(), std::forward<decltype(args)>(args)...}));
    }

    /* wrong extension */
    if (ext_check)
    {
        seqan3::test::tmp_filename filename{"var_io_writer_constructor.xyz"};
        EXPECT_THROW((seqan3::var_io::writer{filename.get_path(), std::forward<decltype(args)>(args)...}),
                     seqan3::unhandled_extension_error);
    }
}

TEST(var_io_writer, constructor1_just_filename)
{
    var_io_writer_filename_constructor(true);
    EXPECT_TRUE((std::same_as<decltype(seqan3::var_io::writer{""}), seqan3::var_io::writer<>>));
}

TEST(var_io_writer, constructor1_with_opts)
{
    seqan3::var_io::writer_options opt{.field_ids = custom_field_ids_t{}};
    using control_t = seqan3::var_io::writer<custom_field_ids_t, seqan3::type_list<seqan3::format_vcf>>;

    var_io_writer_filename_constructor(true, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(seqan3::var_io::writer{"", opt}), control_t>));
}

TEST(var_io_writer, constructor2_just_filename_direct_format)
{
    var_io_writer_filename_constructor(false, seqan3::format_vcf{});
    EXPECT_TRUE((std::same_as<decltype(seqan3::var_io::writer{"", seqan3::format_vcf{}}), seqan3::var_io::writer<>>));
}

TEST(var_io_writer, constructor2_with_opts_direct_format)
{
    seqan3::var_io::writer_options opt{.field_ids = custom_field_ids_t{}};
    using control_t = seqan3::var_io::writer<custom_field_ids_t, seqan3::type_list<seqan3::format_vcf>>;

    var_io_writer_filename_constructor(false, seqan3::format_vcf{}, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(seqan3::var_io::writer{"", seqan3::format_vcf{}, opt}), control_t>));
}

TEST(var_io_writer, constructor2_just_filename_format_variant)
{
    std::variant<seqan3::format_vcf> var{};

    var_io_writer_filename_constructor(false, var);
    EXPECT_TRUE((std::same_as<decltype(seqan3::var_io::writer{"", var}), seqan3::var_io::writer<>>));
}

TEST(var_io_writer, constructor2_with_opts_format_variant)
{
    std::variant<seqan3::format_vcf> var{};
    seqan3::var_io::writer_options   opt{.field_ids = custom_field_ids_t{}};
    using control_t = seqan3::var_io::writer<custom_field_ids_t, seqan3::type_list<seqan3::format_vcf>>;

    var_io_writer_filename_constructor(false, var, std::move(opt));
    EXPECT_TRUE((std::same_as<decltype(seqan3::var_io::writer{"", var, std::move(opt)}), control_t>));
}

TEST(var_io_writer, constructor3)
{
    std::ostringstream str;

    EXPECT_NO_THROW((seqan3::var_io::writer{str, seqan3::format_vcf{}}));
    EXPECT_TRUE((std::same_as<decltype(seqan3::var_io::writer{str, seqan3::format_vcf{}}), seqan3::var_io::writer<>>));
}

TEST(var_io_writer, constructor3_with_opts)
{
    std::ostringstream             str;
    seqan3::var_io::writer_options opt{.field_ids = custom_field_ids_t{}};
    using control_t = seqan3::var_io::writer<custom_field_ids_t, seqan3::type_list<seqan3::format_vcf>>;

    EXPECT_NO_THROW((seqan3::var_io::writer{str, seqan3::format_vcf{}, opt}));
    EXPECT_TRUE((std::same_as<decltype(seqan3::var_io::writer{str, seqan3::format_vcf{}, opt}), control_t>));
}

TEST(var_io_writer, constructor4)
{
    std::ostringstream str;

    EXPECT_NO_THROW((seqan3::var_io::writer{std::move(str), seqan3::format_vcf{}}));
    EXPECT_TRUE(
      (std::same_as<decltype(seqan3::var_io::writer{std::move(str), seqan3::format_vcf{}}), seqan3::var_io::writer<>>));
}

TEST(var_io_writer, constructor4_with_opts)
{
    std::ostringstream             str;
    seqan3::var_io::writer_options opt{.field_ids = custom_field_ids_t{}};
    using control_t = seqan3::var_io::writer<custom_field_ids_t, seqan3::type_list<seqan3::format_vcf>>;

    EXPECT_NO_THROW((seqan3::var_io::writer{std::move(str), seqan3::format_vcf{}, opt}));
    EXPECT_TRUE((std::same_as<decltype(seqan3::var_io::writer{std::move(str), seqan3::format_vcf{}, opt}), control_t>));
}

template <size_t i>
void write_record_test_impl()
{
    using rec_t = seqan3::record<seqan3::type_list<std::string,
                                                   std::string,
                                                   std::string,
                                                   std::string,
                                                   std::string,
                                                   std::string,
                                                   std::string,
                                                   std::string,
                                                   std::vector<std::string>,
                                                   seqan3::var_io::record_private_data>,
                                 std::remove_cvref_t<decltype(seqan3::var_io::default_field_ids)>>;

    seqan3::var_io::header hdr{example_from_spec_header};

    seqan3::var_io::record_private_data priv{};

    std::ostringstream     stream{};
    seqan3::var_io::writer writer{stream, seqan3::format_vcf{}};

    if constexpr (i == 3)
        priv.header_ptr = &hdr;
    else
        writer.set_header(hdr);

    rec_t rec1{
      "20",
      "14370",
      "rs6054257",
      "G",
      "A",
      "29",
      "PASS",
      "NS=3;DP=14;AF=0.5;DB;H2",
      {"GT:GQ:DP:HQ", "0|0:48:1:51,51", "1|0:48:8:51,51", "1/1:43:5:.,."},
      priv
    };
    rec_t rec2{
      "20",
      "17330",
      ".",
      "T",
      "A",
      "3",
      "q10",
      "NS=3;DP=11;AF=0.017",
      {"GT:GQ:DP:HQ", "0|0:49:3:58,50", "0|1:3:5:65,3", "0/0:41:3"},
      priv
    };
    rec_t rec3{
      "20",
      "1110696",
      "rs6040355",
      "A",
      "G,T",
      "67",
      "PASS",
      "NS=2;DP=10;AF=0.333,0.667;AA=T;DB",
      {"GT:GQ:DP:HQ", "1|2:21:6:23,27", "2|1:2:0:18,2", "2/2:35:4"},
      priv
    };
    rec_t rec4{
      "20",
      "1230237",
      ".",
      "T",
      ".",
      "47",
      "PASS",
      "NS=3;DP=13;AA=T",
      {"GT:GQ:DP:HQ", "0|0:54:7:56,60", "0|0:48:4:51,51", "0/0:61:2"},
      priv
    };
    rec_t rec5{
      "20",
      "1234567",
      "microsat1",
      "GTC",
      "G,GTCT",
      "50",
      "PASS",
      "NS=3;DP=9;AA=G",
      {"GT:GQ:DP", "0/1:35:4", "0/2:17:2", "1/1:40:3"},
      priv
    };

    if constexpr (i == 0 || i == 3)
    {
        writer.push_back(rec1);
        writer.push_back(rec2);
        writer.push_back(rec3);
        writer.push_back(rec4);
        writer.push_back(rec5);
    }
    else if constexpr (i == 1)
    {
        auto it = writer.begin();
        it      = rec1;
        it      = rec2;
        it      = rec3;
        it      = rec4;
        it      = rec5;
    }
    else if constexpr (i == 2)
    {
        auto it = writer.begin();
        *it     = rec1;
        *it     = rec2;
        *it     = rec3;
        *it     = rec4;
        *it     = rec5;
    }

    EXPECT_EQ(stream.str(), example_from_spec_header_regenerated_no_IDX + example_from_spec_records);
}

TEST(var_io_writer, push_back_record)
{
    write_record_test_impl<0>();
}

TEST(var_io_writer, assign_to_iterator)
{
    write_record_test_impl<1>();
}

TEST(var_io_writer, assign_to_deref_iterator)
{
    write_record_test_impl<2>();
}

TEST(var_io_writer, implicit_header)
{
    write_record_test_impl<3>();
}

TEST(var_io_writer, emplace_back)
{
    std::ostringstream     stream{};
    seqan3::var_io::writer writer{stream, seqan3::format_vcf{}};

    writer.set_header(seqan3::var_io::header{example_from_spec_header});

    using gs = seqan3::var_io::genotypes_as_strings<>;

    writer.emplace_back("20",
                        "14370",
                        "rs6054257",
                        "G",
                        "A",
                        "29",
                        "PASS",
                        "NS=3;DP=14;AF=0.5;DB;H2",
                        gs{"GT:GQ:DP:HQ", "0|0:48:1:51,51", "1|0:48:8:51,51", "1/1:43:5:.,."},
                        seqan3::var_io::record_private_data{});
    writer.emplace_back("20",
                        "17330",
                        ".",
                        "T",
                        "A",
                        "3",
                        "q10",
                        "NS=3;DP=11;AF=0.017",
                        gs{"GT:GQ:DP:HQ", "0|0:49:3:58,50", "0|1:3:5:65,3", "0/0:41:3"},
                        seqan3::var_io::record_private_data{});
    writer.emplace_back("20",
                        "1110696",
                        "rs6040355",
                        "A",
                        "G,T",
                        "67",
                        "PASS",
                        "NS=2;DP=10;AF=0.333,0.667;AA=T;DB",
                        gs{"GT:GQ:DP:HQ", "1|2:21:6:23,27", "2|1:2:0:18,2", "2/2:35:4"},
                        seqan3::var_io::record_private_data{});
    writer.emplace_back("20",
                        "1230237",
                        ".",
                        "T",
                        ".",
                        "47",
                        "PASS",
                        "NS=3;DP=13;AA=T",
                        gs{"GT:GQ:DP:HQ", "0|0:54:7:56,60", "0|0:48:4:51,51", "0/0:61:2"},
                        seqan3::var_io::record_private_data{});
    writer.emplace_back("20",
                        "1234567",
                        "microsat1",
                        "GTC",
                        "G,GTCT",
                        "50",
                        "PASS",
                        "NS=3;DP=9;AA=G",
                        gs{"GT:GQ:DP", "0/1:35:4", "0/2:17:2", "1/1:40:3"},
                        seqan3::var_io::record_private_data{});

    EXPECT_EQ(stream.str(), example_from_spec_header_regenerated_no_IDX + example_from_spec_records);
}

TEST(var_io_writer, minimal_fields)
{
    std::ostringstream     stream{};
    seqan3::var_io::writer writer{stream,
                                  seqan3::format_vcf{},
                                  seqan3::var_io::writer_options{.field_ids = custom_field_ids_t{}}};

    writer.set_header(seqan3::var_io::header{example_from_spec_header});

    writer.emplace_back("20", "14370", "G");
    writer.emplace_back("20", "17330", "T");
    writer.emplace_back("20", "1110696", "A");
    writer.emplace_back("20", "1230237", "T");
    writer.emplace_back("20", "1234567", "GTC");

    std::string compare = example_from_spec_header_regenerated_no_IDX;
    compare += minimal_field_rows;
    EXPECT_EQ(stream.str(), compare);
}

TEST(var_io_writer, no_header1)
{
    std::ostringstream     stream{};
    seqan3::var_io::writer writer{stream, seqan3::format_vcf{}};

    EXPECT_THROW((writer.emplace_back("20",
                                      "14370",
                                      ".",
                                      "G",
                                      ".",
                                      ".",
                                      ".",
                                      ".",
                                      ".\t.\t.\t.",
                                      seqan3::var_io::record_private_data{})),
                 std::runtime_error);
}

TEST(var_io_writer, no_header2) // this triggers a different path were the private data is not available
{
    std::ostringstream     stream{};
    auto                   fields_without_private = seqan3::tag<seqan3::field::chrom,
                                              seqan3::field::pos,
                                              seqan3::field::id,
                                              seqan3::field::ref,
                                              seqan3::field::alt,
                                              seqan3::field::qual,
                                              seqan3::field::filter,
                                              seqan3::field::info,
                                              seqan3::field::genotypes>;
    seqan3::var_io::writer writer{stream,
                                  seqan3::format_vcf{},
                                  seqan3::var_io::writer_options{.field_ids = fields_without_private}};

    EXPECT_THROW((writer.emplace_back("20", "14370", ".", "G", ".", ".", ".", ".", ".\t.\t.\t.")), std::runtime_error);
}

TEST(var_io_writer, compression)
{
    std::ostringstream stream{};

    {
        seqan3::var_io::writer writer{
          stream,
          seqan3::format_vcf{},
          seqan3::var_io::writer_options{
            .stream_options = seqan3::transparent_ostream_options{.compression = seqan3::compression_format::bgzf}}};

        writer.set_header(seqan3::var_io::header{example_from_spec_header});

        using gs = seqan3::var_io::genotypes_as_strings<>;

        writer.emplace_back("20",
                            "14370",
                            "rs6054257",
                            "G",
                            "A",
                            "29",
                            "PASS",
                            "NS=3;DP=14;AF=0.5;DB;H2",
                            gs{"GT:GQ:DP:HQ", "0|0:48:1:51,51", "1|0:48:8:51,51", "1/1:43:5:.,."},
                            seqan3::var_io::record_private_data{});
        writer.emplace_back("20",
                            "17330",
                            ".",
                            "T",
                            "A",
                            "3",
                            "q10",
                            "NS=3;DP=11;AF=0.017",
                            gs{"GT:GQ:DP:HQ", "0|0:49:3:58,50", "0|1:3:5:65,3", "0/0:41:3"},
                            seqan3::var_io::record_private_data{});
        writer.emplace_back("20",
                            "1110696",
                            "rs6040355",
                            "A",
                            "G,T",
                            "67",
                            "PASS",
                            "NS=2;DP=10;AF=0.333,0.667;AA=T;DB",
                            gs{"GT:GQ:DP:HQ", "1|2:21:6:23,27", "2|1:2:0:18,2", "2/2:35:4"},
                            seqan3::var_io::record_private_data{});
        writer.emplace_back("20",
                            "1230237",
                            ".",
                            "T",
                            ".",
                            "47",
                            "PASS",
                            "NS=3;DP=13;AA=T",
                            gs{"GT:GQ:DP:HQ", "0|0:54:7:56,60", "0|0:48:4:51,51", "0/0:61:2"},
                            seqan3::var_io::record_private_data{});
        writer.emplace_back("20",
                            "1234567",
                            "microsat1",
                            "GTC",
                            "G,GTCT",
                            "50",
                            "PASS",
                            "NS=3;DP=9;AA=G",
                            gs{"GT:GQ:DP", "0/1:35:4", "0/2:17:2", "1/1:40:3"},
                            seqan3::var_io::record_private_data{});
    }

    std::string str = stream.str();
    EXPECT_TRUE(str.starts_with("\x1f\x8b\x08")); // Gzip header

    std::istringstream          control_stream{str};
    seqan3::transparent_istream decompressor{control_stream};
    EXPECT_RANGE_EQ(seqan3::views::istreambuf(decompressor),
                    example_from_spec_header_regenerated_no_IDX + example_from_spec_records);
}
