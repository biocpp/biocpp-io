// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>

#include <bio/format/fasta_input_handler.hpp>

using seqan3::                             operator""_dna5;
using seqan3::                             operator""_phred42;
using std::literals::string_view_literals::operator""sv;

// ----------------------------------------------------------------------------
// fixture
// ----------------------------------------------------------------------------

struct read : public ::testing::Test
{
    using default_rec_t = bio::record<seqan3::vtag_t<bio::field::id, bio::field::seq>,
                                      std::string_view,
                                      decltype(std::string_view{} | seqan3::views::char_to<seqan3::dna5>)>;

    std::vector<std::string> ids{
      {"ID1"},
      {"ID2"},
      {"ID3 lala"},
    };

    std::vector<std::vector<seqan3::dna5>> seqs{
      {"ACGTTTTTTTTTTTTTTT"_dna5},
      {"ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5},
      {"ACGTTTA"_dna5},
    };

    std::vector<std::vector<seqan3::phred42>> quals{
      {"!##$%&'()*+,-./++-"_phred42},
      {"!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE"_phred42},
      {"!!!!!!!"_phred42},
    };

    template <typename id_t, typename seq_t>
    void do_read_test_impl(std::string const & input)
    {
        std::istringstream istream{input};

        bio::format_input_handler<bio::fasta> input_handler{istream};

        bio::record<seqan3::vtag_t<bio::field::id, bio::field::seq>, id_t, seq_t> rec;

        for (unsigned i = 0; i < 3; ++i)
        {
            input_handler.parse_next_record_into(rec);
            if constexpr (std::same_as<std::ranges::range_value_t<seq_t>, char>)
            {
                EXPECT_RANGE_EQ(rec.seq() | seqan3::views::char_to<seqan3::dna5>, seqs[i]);
            }
            else
            {
                EXPECT_RANGE_EQ(rec.seq(), seqs[i]);
            }
            EXPECT_RANGE_EQ(rec.id(), ids[i]);
        }
    }

    void do_read_test(std::string const & input)
    {
        /* containers */
        do_read_test_impl<std::string, std::string>(input);
        do_read_test_impl<std::string, std::vector<seqan3::dna5>>(input);

        /* views */
        do_read_test_impl<std::string_view, std::string_view>(input);
        do_read_test_impl<std::string_view, decltype(std::string_view{} | seqan3::views::char_to<seqan3::dna5>)>(input);
    }
};

// ----------------------------------------------------------------------------
// simple tests
// ----------------------------------------------------------------------------

TEST_F(read, simple)
{
    std::string input =
      R"raw(>ID1
ACGTTTTTTTTTTTTTTT
>ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>ID3 lala
ACGTTTA
)raw";

    do_read_test(input);
}

TEST_F(read, no_trailing_newline)
{
    std::string input =
      R"raw(>ID1
ACGTTTTTTTTTTTTTTT
>ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>ID3 lala
ACGTTTA)raw";

    do_read_test(input);
}

TEST_F(read, empty_lines)
{
    std::string input =
      R"raw(>ID1
ACGTTTTTTTTTTTTTTT

>ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

>ID3 lala
ACGTTTA
)raw";

    do_read_test(input);
}

TEST_F(read, whitespace_in_seq)
{
    std::string input =
      R"raw(>ID1
ACGTTTT
TTTTTTTT
TTT
>ID2
ACGTTTTTT TTTTTTTTT TTTTTTTT TTTTTTTTTT
TTTTTTTTT TTTTTTTTT TTTTTTTT TTTTTTTTTT
TTTTTTTTT T
>ID3 lala
ACGT	TTA
)raw";

    do_read_test(input);
}

TEST_F(read, digits_in_seq)
{
    std::string input =
      R"raw(>ID1
10  ACGTTTTTTTTTTTTTTT
>ID2
  80 ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT  900
1000 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>ID3 lala
ACGT9T5T2A
)raw";

    do_read_test(input);
}

TEST_F(read, old_id_style)
{
    std::string input =
      R"raw(;ID1
ACGTTTTTTTTTTTTTTT
;ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
;ID3 lala
ACGTTTA
)raw";

    do_read_test(input);
}

// ----------------------------------------------------------------------------
// custom tests
// ----------------------------------------------------------------------------

struct options_t
{
    bool truncate_ids = false;
};

TEST_F(read, truncate_ids_off)
{
    std::string input =
      R"raw(>ID1 foo
ACGTTTTTTTTTTTTTTT
>ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>ID3 lala
ACGTTTA
)raw";

    std::istringstream                    istream{input};
    bio::format_input_handler<bio::fasta> input_handler{istream, options_t{}};
    default_rec_t                         rec;

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID1 foo"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[0]);
    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID2"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[1]);
    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID3 lala"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[2]);
}

TEST_F(read, truncate_ids_on)
{
    std::string input =
      R"raw(>ID1 foo
ACGTTTTTTTTTTTTTTT
>ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>ID3 lala
ACGTTTA
)raw";

    std::istringstream                    istream{input};
    bio::format_input_handler<bio::fasta> input_handler{istream, options_t{.truncate_ids = true}};
    default_rec_t                         rec;

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID1"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[0]);
    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID2"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[1]);
    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID3"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[2]);
}

// ----------------------------------------------------------------------------
// failure
// ----------------------------------------------------------------------------

TEST_F(read, fail_no_input)
{
    std::string input{};

    std::istringstream istream{input};
    EXPECT_THROW(bio::format_input_handler<bio::fasta>{istream}, bio::file_open_error);
}

TEST_F(read, fail_no_id)
{
    std::string input{"foo\nACGT"};

    std::istringstream                    istream{input};
    bio::format_input_handler<bio::fasta> input_handler{istream};
    default_rec_t                         rec;

    EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::parse_error);
}

TEST_F(read, fail_only_id)
{
    std::string input{">foo"};

    std::istringstream                    istream{input};
    bio::format_input_handler<bio::fasta> input_handler{istream};
    default_rec_t                         rec;

    EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::parse_error);
}

TEST_F(read, fail_no_seq)
{
    std::string input{">foo\n>bar"};

    std::istringstream                    istream{input};
    bio::format_input_handler<bio::fasta> input_handler{istream};
    default_rec_t                         rec;

    EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::parse_error);
}

// TODO activate after switch to seqan3::views::char_strictly_to
// TEST_F(read, fail_illegal_alphabet)
// {
//     std::string input{">foo\nFOOBAR\n"};
//
//     std::istringstream istream{input};
//     bio::format_input_handler<bio::fasta> input_handler{istream};
//     using rec_t = bio::record<seqan3::vtag_t<bio::field::id, bio::field::seq>, std::string_view,
//     std::vector<seqan3::dna5>>; rec_t rec;
//
//     EXPECT_THROW(input_handler.parse_next_record_into(rec), seqan3::invalid_char_assignment);
// }
