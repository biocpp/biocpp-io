// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream/byte.hpp>
#include <seqan3/core/debug_stream/optional.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/core/debug_stream/variant.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/input_format_concept.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sam_file/output_format_concept.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include <bio/map_io/reader_options.hpp>

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;
using seqan3::operator""_phred42;
using seqan3::operator""_tag;

// global variables for reuse
bio::map_io::reader_options default_options{};
constexpr std::string_view  warning_prefix{"[B.I.O sam format warning in line "};

struct sam_file_data : public ::testing::Test
{
    sam_file_data()
    {
        ref_sequences = std::vector<seqan3::dna5_vector>{ref_seq};
        ref_ids       = std::vector<std::string>{ref_id};
        header        = seqan3::sam_file_header{ref_ids};
        header.ref_id_info.emplace_back(ref_seq.size(), "");
        header.ref_dict[header.ref_ids()[0]] = 0; // set up header which is otherwise done on file level
    }

    std::vector<seqan3::dna5_vector> seqs{"ACGT"_dna5, "AGGCTGNAG"_dna5, "GGAGTATA"_dna5};

    std::vector<std::string> ids{"read1", "read2", "read3"};

    std::vector<std::vector<seqan3::phred42>> quals{
      {"!##$"_phred42},
      {"!##$&'()*"_phred42},
      {"!!*+,-./"_phred42},
    };

    std::vector<int32_t> offsets{1, 0, 1};

    seqan3::dna5_vector ref_seq = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;

    std::string ref_id = "ref";

    std::vector<int32_t> positions // 1-based in b.i.o. 0-based in seqan3
      {1, 2, 3};

    // clang-format off
    std::vector<std::vector<seqan3::cigar>> cigars
    {
        // cigar read1
        {{1, 'S'_cigar_operation},
         {1, 'M'_cigar_operation},
         {1, 'D'_cigar_operation},
         {1, 'M'_cigar_operation},
         {1, 'I'_cigar_operation}},
        // cigar read2
        {{1, 'H'_cigar_operation},
         {7, 'M'_cigar_operation},
         {1, 'D'_cigar_operation},
         {1, 'M'_cigar_operation},
         {1, 'S'_cigar_operation},
         {2, 'H'_cigar_operation}},
        // cigar read3
        {{1, 'S'_cigar_operation},
         {1, 'M'_cigar_operation},
         {1, 'P'_cigar_operation},
         {1, 'M'_cigar_operation},
         {1, 'I'_cigar_operation},
         {1, 'M'_cigar_operation},
         {1, 'I'_cigar_operation},
         {1, 'D'_cigar_operation},
         {1, 'M'_cigar_operation},
         {1, 'S'_cigar_operation}}
    };
    // clang-format on

    std::vector<bio::map_io::sam_flag> flags{bio::map_io::sam_flag{41u},
                                             bio::map_io::sam_flag{42u},
                                             bio::map_io::sam_flag{43u}};

    std::vector<uint8_t> mapqs{61u, 62u, 63u};

    std::vector<std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>> mates // position 1-based
      {
        {0, 10, 300},
        {0, 10, 300},
        {0, 10, 300}
    };

    std::vector<bio::map_io::sam_tag_dictionary> tag_dicts{bio::map_io::sam_tag_dictionary{},
                                                           bio::map_io::sam_tag_dictionary{},
                                                           bio::map_io::sam_tag_dictionary{}};

    std::vector<seqan3::dna5_vector>                  ref_sequences{};
    std::vector<std::string>                          ref_ids{};
    seqan3::sam_file_header<std::vector<std::string>> header{};
};

template <typename format_t>
struct sam_file_read : public sam_file_data
{};

TYPED_TEST_SUITE_P(sam_file_read);

// ----------------------------------------------------------------------------
// sam_file_read
// ----------------------------------------------------------------------------

TYPED_TEST_P(sam_file_read, full_data_set)
{
    // prepare tag dictionary
    this->tag_dicts[0]["NM"_tag] = -7;
    this->tag_dicts[0]["AS"_tag] = 2;
    this->tag_dicts[0]["CC"_tag] = 300;
    this->tag_dicts[0]["cc"_tag] = -300;
    this->tag_dicts[0]["aa"_tag] = 'c';
    this->tag_dicts[0]["ff"_tag] = 3.1f;
    this->tag_dicts[0]["zz"_tag] = "str";
    this->tag_dicts[1]["bc"_tag] = std::vector<int8_t>{-3};
    this->tag_dicts[1]["bC"_tag] = std::vector<uint8_t>{3u, 200u};
    this->tag_dicts[1]["bs"_tag] = std::vector<int16_t>{-3, 200, -300};
    this->tag_dicts[1]["bS"_tag] = std::vector<uint16_t>{300u, 40u, 500u};
    this->tag_dicts[1]["bi"_tag] = std::vector<int32_t>{-3, 200, -66000};
    this->tag_dicts[1]["bI"_tag] = std::vector<uint32_t>{294967296u};
    this->tag_dicts[1]["bf"_tag] = std::vector<float>{3.5f, 0.1f, 43.8f};
    this->tag_dicts[1]["bH"_tag] = std::vector<std::byte>{std::byte{0x1A}, std::byte{0xE3}, std::byte{0x01}};

    // read data
    typename TestFixture::stream_type istream{this->verbose_reads_input};
    using record_t = bio::record<decltype(default_options.field_ids), decltype(default_options.field_types)>;
    bio::format_input_handler<bio::sam> input_handler{istream, default_options};
    record_t                            rec;

    for (unsigned i = 0; i < 3; ++i)
    {
        input_handler.parse_next_record_into(rec);

        EXPECT_EQ(rec.id(), this->ids[i]) << "failed at i = " << i << std::endl;
        EXPECT_EQ(rec.flag(), this->flags[i]) << "failed at i = " << i << std::endl;
        EXPECT_EQ(rec.rname(), this->ref_id) << "failed at i = " << i << std::endl;
        EXPECT_EQ(rec.pos(), this->positions[i]) << "failed at i = " << i << std::endl;
        EXPECT_EQ(rec.mapq(), this->mapqs[i]) << "failed at i = " << i << std::endl;
        EXPECT_RANGE_EQ(rec.cigar(), this->cigars[i]);
        EXPECT_EQ(rec.rnext(), this->ref_id) << "failed at i = " << i << std::endl;
        EXPECT_EQ(rec.pnext(), std::get<1>(this->mates[i])) << "failed at i = " << i << std::endl;
        EXPECT_EQ(rec.tlen(), std::get<2>(this->mates[i])) << "failed at i = " << i << std::endl;
        EXPECT_RANGE_EQ(rec.seq(), this->seqs[i]);
        EXPECT_RANGE_EQ(rec.qual(), this->quals[i]);
        EXPECT_EQ(rec.tags(), this->tag_dicts[i]) << "failed at i = " << i << std::endl;
    }
}

TYPED_TEST_P(sam_file_read, all_missing_data)
{
    typename TestFixture::stream_type istream{this->empty_input};

    using record_t = bio::record<decltype(default_options.field_ids), decltype(default_options.field_types)>;
    bio::format_input_handler<bio::sam> input_handler{istream, default_options};
    record_t                            rec;

    input_handler.parse_next_record_into(rec);

    EXPECT_TRUE(rec.id().empty());
    EXPECT_TRUE(rec.rname().empty());
    EXPECT_TRUE(rec.rnext().empty());
    EXPECT_TRUE(rec.cigar().empty());
    EXPECT_TRUE(rec.seq().empty());
    EXPECT_TRUE(rec.qual().empty());
    EXPECT_TRUE(rec.tags().empty());

    EXPECT_EQ(rec.flag(), bio::map_io::sam_flag{0u});
    EXPECT_EQ(rec.pos(), 0);
    EXPECT_EQ(rec.pnext(), 0);
    EXPECT_EQ(rec.mapq(), 0u);
    EXPECT_EQ(rec.tlen(), 0);
}

TYPED_TEST_P(sam_file_read, select_fields)
{
    typename TestFixture::stream_type istream{this->verbose_reads_input};

    constexpr auto fid   = bio::vtag<bio::field::rname, bio::field::pos>;
    constexpr auto ftype = bio::ttag<std::string_view, int64_t>;

    using record_t = bio::record<std::remove_cvref_t<decltype(fid)>, std::remove_cvref_t<decltype(ftype)>>;
    bio::format_input_handler<bio::sam> input_handler{istream, default_options};
    record_t                            rec;

    for (unsigned i = 0; i < 3; ++i)
    {
        input_handler.parse_next_record_into(rec);
        EXPECT_EQ(rec.rname(), this->ref_id);
        EXPECT_EQ(rec.pos(), this->positions[i]);
    }
}

TYPED_TEST_P(sam_file_read, warning_rname_not_in_header)
{
    typename TestFixture::stream_type istream{this->verbose_reads_input};

    constexpr auto fid   = bio::vtag<bio::field::rname, bio::field::pos>;
    constexpr auto ftype = bio::ttag<std::string_view, int64_t>;

    using record_t = bio::record<std::remove_cvref_t<decltype(fid)>, std::remove_cvref_t<decltype(ftype)>>;
    bio::format_input_handler<bio::sam> input_handler{istream, default_options};
    record_t                            rec;

    std::string warning_msg{warning_prefix.data() + std::string{"1] The reference sequence name \""} + this->ref_id +
                            "\" is not present in the header.\n"};

    for (unsigned i = 0; i < 3; ++i)
    {
        testing::internal::CaptureStderr();
        input_handler.parse_next_record_into(rec);
        EXPECT_EQ(testing::internal::GetCapturedStderr(), ((i == 0) ? warning_msg : std::string{}));

        EXPECT_EQ(rec.rname(), this->ref_id);
        EXPECT_EQ(rec.pos(), this->positions[i]);
    }
}

TYPED_TEST_P(sam_file_read, warning_rnext_not_in_header)
{
    typename TestFixture::stream_type istream{this->verbose_reads_input};

    constexpr auto fid   = bio::vtag<bio::field::rnext, bio::field::pnext>;
    constexpr auto ftype = bio::ttag<std::string_view, int64_t>;

    using record_t = bio::record<std::remove_cvref_t<decltype(fid)>, std::remove_cvref_t<decltype(ftype)>>;
    bio::format_input_handler<bio::sam> input_handler{istream, default_options};
    record_t                            rec;

    std::string warning_msg{warning_prefix.data() + std::string{"1] The mates reference sequence name \""} +
                            this->ref_id + "\" is not present in the header.\n"};

    for (unsigned i = 0; i < 3; ++i)
    {
        testing::internal::CaptureStderr();
        input_handler.parse_next_record_into(rec);
        EXPECT_EQ(testing::internal::GetCapturedStderr(), ((i == 0) ? warning_msg : std::string{}));

        EXPECT_EQ(rec.rnext(), this->ref_id);
        EXPECT_EQ(rec.pnext(), std::get<1>(this->mates[i]));
    }
}

TYPED_TEST_P(sam_file_read, format_error_uneven_hexadecimal_tag)
{
    typename TestFixture::stream_type istream{this->wrong_hexadecimal_tag};

    constexpr auto fid   = bio::vtag<bio::field::tags>;
    constexpr auto ftype = bio::ttag<bio::map_io::sam_tag_dictionary>;

    using record_t = bio::record<std::remove_cvref_t<decltype(fid)>, std::remove_cvref_t<decltype(ftype)>>;
    bio::format_input_handler<bio::sam> input_handler{istream, default_options};
    record_t                            rec;

    EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::format_error);
}

REGISTER_TYPED_TEST_SUITE_P(sam_file_read,
                            full_data_set,
                            all_missing_data,
                            select_fields,
                            warning_rname_not_in_header,
                            warning_rnext_not_in_header,
                            format_error_uneven_hexadecimal_tag);
