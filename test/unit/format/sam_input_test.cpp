// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <bio/format/sam.hpp>
#include <bio/format/sam_input_handler.hpp>

#include "sam_input_test_template.hpp"

// ----------------------------------------------------------------------------
// fixture
// ----------------------------------------------------------------------------

template <>
struct sam_file_read<seqan3::format_sam> : public sam_file_data
{
    // -----------------------------------------------------------------------------------------------------------------
    // formatted input
    // -----------------------------------------------------------------------------------------------------------------

    using stream_type = std::istringstream;

    std::string big_header_input{
      R"(@HD	VN:1.6	SO:coordinate	SS:coordinate:queryname	GO:none
@PG	ID:qc	PN:quality_control	CL:qc -f file1	DS:trim reads with low qual	VN:1.0.0
@PG	ID:novoalign	PN:novoalign	VN:V3.02.07	CL:novoalign -d /path/hs37d5.ndx -f /path/file.fastq.gz	PP:qc
@SQ	SN:ref	LN:249250621
@SQ	SN:ref2	LN:243199373	AS:hs37d5
@RG	ID:U0a_A2_L1	PL:illumina	PU:1	LB:1	SM:NA12878
@RG	ID:U0a_A2_L2	PL:illumina	SM:NA12878	PU:1	LB:1
@CO	Tralalalalalala this is a comment
)"};

    std::string simple_three_reads_input{
      R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	1H7M1D1M1S2H	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1P1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)"};

    std::string verbose_reads_input{
      // "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:34\n" // todo: check why adding this segfaults, adding this removes warnings
      "read1\t41\tref\t1\t61\t1S1M1D1M1I\t=\t10\t300\tACGT\t!##$\taa:A:c"
      "\tNM:i:-7"
      "\tAS:i:2"
      "\tff:f:3.1"
      "\tzz:Z:str"
      "\tCC:i:300"
      "\tcc:i:-300\n"
      "read2\t42\tref\t2\t62\t1H7M1D1M1S2H\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\tbc:B:c,-3"
      "\tbC:B:C,3,200"
      "\tbs:B:s,-3,200,-300"
      "\tbS:B:S,300,40,500"
      "\tbi:B:i,-3,200,-66000"
      "\tbI:B:I,294967296"
      "\tbf:B:f,3.5,0.1,43.8"
      "\tbH:H:1AE301\n"
      "read3\t43\tref\t3\t63\t1S1M1P1M1I1M1I1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"};

    std::string empty_input{"*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"};

    std::string empty_cigar{"read1\t41\tref\t1\t61\t*\tref\t10\t300\tACGT\t!##$\n"};

    std::string unknown_ref{
      "read1\t41\traf\t1\t61\t1S1M1D1M1I\t=\t10\t300\tACGT\t!##$\taa:A:c\tAS:i:2\tff:f:3.1\tzz:Z:str\n"};

    std::string unknown_ref_header{"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:34\n*\t0\tunknown_ref\t1\t0\t4M\t*\t0\t0\tAAAA\t*\n"};

    // -----------------------------------------------------------------------------------------------------------------
    // formatted output
    // -----------------------------------------------------------------------------------------------------------------

    std::string simple_three_reads_output{// compared to simple_three_reads_input this has no hard clipping
                                          R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1P1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)"};

    std::string verbose_output{
      R"(@HD	VN:1.6	SO:unknown	GO:none
@SQ	SN:ref	LN:34	AN:other_name
@RG	ID:group1	DS:more info
@PG	ID:prog1	PN:cool_program	CL:./prog1	PP:a	DS:b	VN:c
@CO	This is a comment.
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	CC:i:300	NM:i:-7	aa:A:c	cc:i:-300	ff:f:3.1	zz:Z:str
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	bC:B:C,3,200	bI:B:I,294967296	bS:B:S,300,40,500	bc:B:c,-3	bf:B:f,3.5,0.1,43.8	bi:B:i,-3,200,-66000	bs:B:s,-3,200,-300
read3	43	ref	3	63	1S1M1P1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)"};

    std::string special_output{
      R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	*	1	61	1S1M1D1M1I	*	0	0	ACGT	!##$
)"};

    std::string wrong_hexadecimal_tag{
      R"(read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	bH:H:1AE30
)"};
};

// ---------------------------------------------------------------------------------------------------------------------
// parametrized tests
// ---------------------------------------------------------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_SUITE_P(sam, sam_file_read, seqan3::format_sam, );
// INSTANTIATE_TYPED_TEST_SUITE_P(sam, sam_file_write, seqan3::format_sam, );

// ---------------------------------------------------------------------------------------------------------------------
// SAM specifics
// ---------------------------------------------------------------------------------------------------------------------

struct sam_input_test : public sam_file_data
{};

TEST_F(sam_input_test, no_hd_line_in_header)
{
    // the header line (@HD) is optional
    std::istringstream istream{std::string{"@SQ\tSN:ref\tLN:34\nread1\t41\tref\t1\t61\t*\tref\t10\t300\tACGT\t!##$\n"}};

    using record_t = bio::record<std::remove_cvref_t<decltype(bio::vtag<bio::field::id>)>,
                                 std::remove_cvref_t<decltype(bio::ttag<std::string_view>)>>;
    bio::format_input_handler<bio::sam> input_handler{istream, default_options /*defined in sam_input_test_template*/};
    record_t                            rec;

    EXPECT_NO_THROW(input_handler.parse_next_record_into(rec));

    EXPECT_EQ(rec.id(), std::string{"read1"});
}

TEST_F(sam_input_test, windows_file)
{
    std::istringstream istream(std::string{"read1\t41\tref\t1\t61\t*\tref\t10\t300\tACGT\t!##$\r\n"});

    using record_t = bio::record<std::remove_cvref_t<decltype(bio::vtag<bio::field::id>)>,
                                 std::remove_cvref_t<decltype(bio::ttag<std::string_view>)>>;
    bio::format_input_handler<bio::sam> input_handler{istream, default_options /*defined in sam_input_test_template*/};
    record_t                            rec;

    EXPECT_NO_THROW(input_handler.parse_next_record_into(rec));

    EXPECT_EQ(rec.id(), std::string{"read1"});
}

TEST_F(sam_input_test, seqan3_issue2195)
{ // see issue https://github.com/seqan/seqan3/issues/2195
    using seqan3::operator""_phred42;

    auto ftype = bio::ttag<std::string_view, decltype(std::string_view{} | seqan3::views::char_to<seqan3::phred42>)>;
    using record_t = bio::record<std::remove_cvref_t<decltype(bio::vtag<bio::field::id, bio::field::qual>)>,
                                 std::remove_cvref_t<decltype(ftype)>>;
    record_t rec;

    {
        std::istringstream istream{"*r1\t4\t1\t10\t0\t5M\t=\t136097\t-121\tACTGA\t*9<9;\tNM:i:1\tMQ:i:0\n"};
        bio::format_input_handler<bio::sam> input_handler{istream, default_options};
        EXPECT_NO_THROW(input_handler.parse_next_record_into(rec));

        EXPECT_RANGE_EQ(rec.id(), std::string{"*r1"});
        EXPECT_RANGE_EQ(rec.qual(), "*9<9;"_phred42);
    }

    {
        std::istringstream                  istream{"*\t4\t1\t10\t0\t2M\t=\t136097\t-121\tAC\t*1\tNM:i:1\tMQ:i:0\n"};
        bio::format_input_handler<bio::sam> input_handler{istream, default_options};
        EXPECT_NO_THROW(input_handler.parse_next_record_into(rec));

        EXPECT_RANGE_EQ(rec.id(), std::string{""});
        EXPECT_RANGE_EQ(rec.qual(), "*1"_phred42);
    }
}

struct sam_input_format_error : public sam_file_data
{};

TEST_F(sam_input_format_error, illegal_character_in_seq)
{
    std::istringstream istream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\tAC!T\t*\n"));

    using record_t = bio::record<std::remove_cvref_t<decltype(bio::vtag<bio::field::seq>)>,
                                 std::remove_cvref_t<decltype(bio::ttag<seqan3::dna5_vector>)>>;
    bio::format_input_handler<bio::sam> input_handler{istream, default_options /*defined in sam_input_test_template*/};
    record_t                            rec;

    EXPECT_THROW(input_handler.parse_next_record_into(rec), seqan3::invalid_char_assignment);
}

TEST_F(sam_input_format_error, invalid_arithmetic_value)
{
    using record_t = bio::record<std::remove_cvref_t<decltype(bio::vtag<bio::field::pos, bio::field::pnext>)>,
                                 std::remove_cvref_t<decltype(bio::ttag<int32_t, int32_t>)>>;
    record_t rec;

    // invalid value
    {
        std::istringstream                  istream(std::string("*\t0\t*\t1abc\t0\t*\t*\t0\t0\t*\t*\n"));
        bio::format_input_handler<bio::sam> input_handler{istream, default_options};
        EXPECT_THROW(input_handler.parse_next_record_into(rec), std::runtime_error);
    }
    // overflow error
    {
        std::istringstream istream = std::istringstream(std::string("*\t0\t*\t2147483650\t0\t*\t*\t0\t0\t*\t*\n"));
        bio::format_input_handler<bio::sam> input_handler{istream, default_options};
        EXPECT_THROW(input_handler.parse_next_record_into(rec), std::runtime_error);
    }
    // negative value as pos - do we want this semantic check?
    // {
    //     std::istringstream istream = std::istringstream(std::string("*\t0\t*\t-3\t0\t*\t*\t0\t0\t*\t*\n"));
    //     bio::format_input_handler<bio::sam> input_handler{istream, default_options};
    //     EXPECT_THROW(input_handler.parse_next_record_into(rec), std::runtime_error);
    // }

    // negative value as mate mapping position - do we want this semantic check?
    // {
    //     std::istringstream istream = std::istringstream(std::string("*\t0\t*\t0\t0\t*\t*\t-3\t0\t*\t*\n"));
    //     bio::format_input_handler<bio::sam> input_handler{istream, default_options};
    //     EXPECT_THROW(input_handler.parse_next_record_into(rec), std::runtime_error);
    // }
}

TEST_F(sam_input_format_error, invalid_cigar)
{
    using record_t = bio::record<std::remove_cvref_t<decltype(bio::vtag<bio::field::cigar>)>,
                                 std::remove_cvref_t<decltype(bio::ttag<std::vector<seqan3::cigar>>)>>;
    record_t rec;

    // unkown operation
    {
        std::istringstream                  istream(std::string("*\t0\t*\t0\t0\t5Z\t*\t0\t0\t*\t*\n"));
        bio::format_input_handler<bio::sam> input_handler{istream, default_options};
        EXPECT_THROW(input_handler.parse_next_record_into(rec), seqan3::invalid_char_assignment);
    }
    // negative number as operation count
    {
        std::istringstream istream = std::istringstream(std::string("*\t0\t*\t0\t0\t-5M\t*\t0\t0\t*\t*\n"));
        bio::format_input_handler<bio::sam> input_handler{istream, default_options};
        EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::format_error);
    }
    // negative number as operation count - long cigar
    {
        std::istringstream istream = std::istringstream(std::string("*\t0\t*\t0\t0\t3S4M1I-5M2D2M\t*\t0\t0\t*\t*\n"));
        bio::format_input_handler<bio::sam> input_handler{istream, default_options};
        EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::format_error);
    }
}

TEST_F(sam_input_format_error, invalid_sam_tag_format)
{
    using record_t = bio::record<std::remove_cvref_t<decltype(bio::vtag<bio::field::tags>)>,
                                 std::remove_cvref_t<decltype(bio::ttag<bio::map_io::sam_tag_dictionary>)>>;
    record_t rec;

    // type identifier is wrong
    {
        std::istringstream                  istream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\tNM:X:3\n"));
        bio::format_input_handler<bio::sam> input_handler{istream, default_options};
        EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::format_error);
    }
    // Array subtype identifier is wrong
    {
        std::istringstream istream = std::istringstream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\tNM:B:x3,4\n"));
        bio::format_input_handler<bio::sam> input_handler{istream, default_options};
        EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::format_error);
    }
}
