// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <bio/map_io/header.hpp>

TEST(map_io_header, ctors)
{
    // todo
}

TEST(map_io_header, read_full_header)
{
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

    bio::map_io::header header{};
    EXPECT_NO_THROW(header.read(big_header_input));

    EXPECT_EQ(header.format_version, "1.6");
    EXPECT_EQ(header.sorting, "coordinate");
    EXPECT_EQ(header.subsorting, "coordinate:queryname");
    EXPECT_EQ(header.grouping, "none");

    EXPECT_EQ(header.program_infos[0].id, "qc");
    EXPECT_EQ(header.program_infos[0].name, "quality_control");
    EXPECT_EQ(header.program_infos[0].version, "1.0.0");
    EXPECT_EQ(header.program_infos[0].description, "trim reads with low qual");
    EXPECT_EQ(header.program_infos[0].previous, "");
    EXPECT_EQ(header.program_infos[0].command_line_call, "qc -f file1");
    EXPECT_EQ(header.program_infos[1].id, "novoalign");
    EXPECT_EQ(header.program_infos[1].name, "novoalign");
    EXPECT_EQ(header.program_infos[1].version, "V3.02.07");
    EXPECT_EQ(header.program_infos[1].description, "");
    EXPECT_EQ(header.program_infos[1].previous, "qc");
    EXPECT_EQ(header.program_infos[1].command_line_call, "novoalign -d /path/hs37d5.ndx -f /path/file.fastq.gz");

    std::string id1{"ref"};
    std::string id2{"ref2"};

    EXPECT_EQ(header.rnames_info().size(), 2u);
    EXPECT_EQ(header.rnames_info()[(header.rname_to_pos().at(id1))],
              (std::tuple<uint32_t, std::string>{249250621u, ""}));
    EXPECT_EQ(header.rnames_info()[(header.rname_to_pos().at(id2))],
              (std::tuple<uint32_t, std::string>{243199373u, "AS:hs37d5"}));

    EXPECT_EQ(header.read_groups[0],
              (std::pair<std::string, std::string>{"U0a_A2_L1", "PL:illumina\tPU:1\tLB:1\tSM:NA12878"}));
    EXPECT_EQ(header.read_groups[1],
              (std::pair<std::string, std::string>{"U0a_A2_L2", "PL:illumina\tSM:NA12878\tPU:1\tLB:1"}));

    EXPECT_EQ(header.comments[0], "Tralalalalalala this is a comment");
}

TEST(map_io_header, independent_order)
{
    // order of tags does not matter
    std::string header_str{
      "@HD\tGO:none\tSO:coordinate\tVN:1.6\tSS:coordinate:queryname\n"
      "@CO\tTralalalalalala this is a comment\n"
      "@PG\tPN:novoalign\tPP:qc\tID:novoalign\tVN:V3.02.07\tCL:novoalign -d /hs37d5.ndx -f /file.fastq.gz\n"
      "@SQ\tAS:hs37d5\tSN:ref2\tLN:243199373\n"
      "@RG\tLB:1\tSM:NA12878\tPL:illumina\tPU:1\tID:U0a_A2_L1\n"};

    bio::map_io::header header{};
    EXPECT_NO_THROW(header.read(header_str));

    EXPECT_EQ(header.format_version, "1.6");
    EXPECT_EQ(header.sorting, "coordinate");
    EXPECT_EQ(header.subsorting, "coordinate:queryname");
    EXPECT_EQ(header.grouping, "none");

    EXPECT_EQ(header.program_infos[0].id, "novoalign");
    EXPECT_EQ(header.program_infos[0].name, "novoalign");
    EXPECT_EQ(header.program_infos[0].version, "V3.02.07");
    EXPECT_EQ(header.program_infos[0].description, "");
    EXPECT_EQ(header.program_infos[0].previous, "qc");
    EXPECT_EQ(header.program_infos[0].command_line_call, "novoalign -d /hs37d5.ndx -f /file.fastq.gz");

    std::string id2{"ref2"};
    EXPECT_EQ(header.rnames_info().size(), 1u);
    EXPECT_EQ(header.rnames_info()[(header.rname_to_pos().at(id2))],
              (std::tuple<uint32_t, std::string>{243199373u, "AS:hs37d5"}));

    EXPECT_EQ(header.read_groups[0],
              (std::pair<std::string, std::string>{"U0a_A2_L1", "LB:1\tSM:NA12878\tPL:illumina\tPU:1"}));

    EXPECT_EQ(header.comments[0], "Tralalalalalala this is a comment");
}

TEST(map_io_header, issue2423)
{
    // former seqan3 issue https://github.com/seqan/seqan3/pull/2423
    std::string many_refs_input{[]()
                                {
                                    std::string result{"@HD\tVN:1.6\n"};
                                    for (size_t i = 0; i < 64; ++i)
                                        result += "@SQ\tSN:ref_" + std::to_string(i) + "\tLN:100\n";
                                    return result;
                                }()};

    bio::map_io::header header{};
    EXPECT_NO_THROW(header.read(many_refs_input));

    EXPECT_EQ(header.rnames().size(), 64u);
    EXPECT_EQ(header.rnames_info().size(), 64u);
    EXPECT_EQ(header.rname_to_pos().size(), 64u);
}

TEST(map_io_header, tag_LN_maximum_value)
{
    // maximum LN value is 2^31-1
    std::string header_str{"@SQ\tSN:ref\tLN:2147483647\n"};

    bio::map_io::header header{};
    EXPECT_NO_THROW(header.read(header_str));

    EXPECT_EQ(header.rnames().size(), 1u);
    EXPECT_EQ(header.rnames_info()[(header.rname_to_pos().at("ref"))],
              (std::tuple<uint32_t, std::string>{2147483647u, ""}));
}

TEST(map_io_header_warnings, user_defined_tags)
{
    // user defined tags should not trigger errors, but print warnings to cerr
    std::string header_str{
      "@HD\tVN:1.6\tVB:user_tag\tSB:user_tag\tGB:user_tag\tpb:user_tag\n"
      "@SQ\tSN:ref2\tLN:243199373\tSB:user_tag\tLB:user_tag\tpb:user_tag\n"
      "@RG\tID:U0a_A2_L1\tIB:user_tag\tpb:user_tag\n"
      "@PG\tID:qc\tIB:user_tag\tPB:user_tag\tCB:user_tag\tDB:user_tag\tVB:user_tag\tpb:user_tag\n"};

    std::string expected_cerr{
      "Unsupported SAM header tag in @HD: VB\n"
      "Unsupported SAM header tag in @HD: SB\n"
      "Unsupported SAM header tag in @HD: GB\n"
      "Unsupported SAM header tag in @HD: pb\n"
      "Unsupported SAM header tag in @PG: IB\n"
      "Unsupported SAM header tag in @PG: PB\n"
      "Unsupported SAM header tag in @PG: CB\n"
      "Unsupported SAM header tag in @PG: DB\n"
      "Unsupported SAM header tag in @PG: VB\n"
      "Unsupported SAM header tag in @PG: pb\n"};

    bio::map_io::header header{};
    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(header.read(header_str));
    EXPECT_EQ(testing::internal::GetCapturedStderr(), expected_cerr);
}

TEST(map_io_header_errors, invalid_tag_HA)
{
    // invalid header record type: @HA
    std::string header_str{"@HA\tthis is not a valid tag\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, invalid_tag_SA)
{
    // invalid header record type: @SA
    std::string header_str{"@SA\tthis is not a valid tag\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, invalid_tag_PA)
{
    // invalid header record type: @PA
    std::string header_str{"@PA\tthis is not a valid tag\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, invalid_tag_RA)
{
    // invalid header record type: @RA
    std::string header_str{"@RA\tthis is not a valid tag\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, invalid_tag_CA)
{
    // invalid header record type: @CA
    std::string header_str{"@CA\tthis is not a valid tag\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, invalid_tag_TT)
{
    // invalid header record type: @TT
    std::string header_str{"@TT\tthis is not a valid tag\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, missing_tag_VN)
{
    // missing VN tag in @HD
    std::string header_str{"@HD\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, missing_tag_SN)
{
    // missing SN tag in @SQ
    std::string header_str{"@SQ\tLN:1\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, missing_tag_LN)
{
    // missing LN tag in @SQ
    std::string header_str{"@SQ\tSN:ref\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

// do we want this semantic check?
// TEST(map_io_header_errors, tag_LN_cannot_be_0)
// {
//     // LN cannot be 0
//     std::string header_str
//     {
//         "@SQ\tSN:ref\tLN:0\n"
//     };

//     bio::map_io::header header{};
//     EXPECT_THROW(header.read(header_str), bio::format_error);
// }

// do we want this semantic check?
// TEST(map_io_header_errors, tag_LN_cannot_be_negative)
// {
//     // LN cannot be negative
//     std::string header_str
//     {
//         "@SQ\tSN:ref\tLN:-1\n"
//     };

//     bio::map_io::header header{};
//     EXPECT_THROW(header.read(header_str), bio::format_error);
// }

TEST(map_io_header_errors, tag_LN_overflow)
{
    // LN exceeds maximum value
    std::string header_str{"@SQ\tSN:ref\tLN:2147483648\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, missing_tag_RG_ID)
{
    // missing ID tag in @RG
    std::string header_str{"@RG\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

TEST(map_io_header_errors, missing_tag_PG_ID)
{
    // missing ID tag in @PG
    std::string header_str{"@PG\n"};

    bio::map_io::header header{};
    EXPECT_THROW(header.read(header_str), bio::format_error);
}

// write

// TEST_F(sam_input_test, write_different_header)
// {
//     std::ostringstream ostream;

//     auto write_header = [&] ()
//     {
//         seqan3::sam_file_output fout{ostream, seqan3::format_sam{}, seqan3::fields<seqan3::field::header_ptr,
//                                                                                    seqan3::field::ref_id,
//                                                                                    seqan3::field::ref_offset>{}};
//         ASSERT_NO_THROW(fout.emplace_back(&header, this->ref_id, 0));
//     };

//     header.sorting = "unsorted";
//     header.grouping = "query";

//     write_header();
//     ostream.flush();
//     EXPECT_EQ(ostream.str(),
//               "@HD\tVN:1.6\tSO:unsorted\tGO:query\n@SQ\tSN:ref\tLN:34\n*\t0\tref\t1\t0\t*\t*\t0\t0\t*\t*\n");

//     ostream = std::ostringstream{};
//     header.sorting = "queryname";
//     header.grouping = "reference";

//     write_header();
//     ostream.flush();
//     EXPECT_EQ(ostream.str(),
//               "@HD\tVN:1.6\tSO:queryname\tGO:reference\n@SQ\tSN:ref\tLN:34\n*\t0\tref\t1\t0\t*\t*\t0\t0\t*\t*\n");

//     ostream = std::ostringstream{};
//     header.sorting = "coordinate";
//     header.subsorting = "query";

//     write_header();
//     ostream.flush();
//     EXPECT_EQ(ostream.str(),
//               "@HD\tVN:1.6\tSO:coordinate\tSS:query\tGO:reference\n@SQ\tSN:ref\tLN:34\n*\t0\tref\t1\t0\t*\t*\t0\t0\t*\t*\n");
// }
