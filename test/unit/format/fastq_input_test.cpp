// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <sstream>

#include <gtest/gtest.h>

#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/quality/phred42.hpp>
#include <bio/test/expect_range_eq.hpp>
#include <bio/test/expect_same_type.hpp>

#include <bio/io/format/fastq_input_handler.hpp>

using namespace bio::alphabet::literals;
using std::literals::string_view_literals::operator""sv;

// ----------------------------------------------------------------------------
// fixture
// ----------------------------------------------------------------------------

struct read : public ::testing::Test
{
    using default_rec_t = bio::io::record<
      bio::meta::vtag_t<bio::io::field::id, bio::io::field::seq, bio::io::field::qual>,
      bio::meta::type_list<std::string_view,
                           decltype(std::string_view{} | bio::views::char_strictly_to<bio::alphabet::dna5>),
                           decltype(std::string_view{} | bio::views::char_strictly_to<bio::alphabet::phred42>)>>;

    std::string default_input =
      R"raw(@ID1
ACGTTTTTTTTTTTTTTT
+
!##$%&'()*+,-./++-
@ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
@ID3 lala
ACGTTTA
+
!!!!!!!
)raw";

    std::vector<std::string> ids{
      {"ID1"},
      {"ID2"},
      {"ID3 lala"},
    };

    std::vector<std::vector<bio::alphabet::dna5>> seqs{
      {"ACGTTTTTTTTTTTTTTT"_dna5},
      {"ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5},
      {"ACGTTTA"_dna5},
    };

    std::vector<std::vector<bio::alphabet::phred42>> quals{
      {"!##$%&'()*+,-./++-"_phred42},
      {"!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE"_phred42},
      {"!!!!!!!"_phred42},
    };

    template <typename id_t, typename seq_t, typename qual_t>
    void do_read_test_impl(std::string const & input)
    {
        std::istringstream istream{input};

        bio::io::format_input_handler<bio::io::fastq> input_handler{istream};

        bio::io::record<bio::meta::vtag_t<bio::io::field::id, bio::io::field::seq, bio::io::field::qual>,
                        bio::meta::type_list<id_t, seq_t, qual_t>>
          rec;

        for (unsigned i = 0; i < 3; ++i)
        {
            input_handler.parse_next_record_into(rec);
            EXPECT_RANGE_EQ(rec.id(), ids[i]);
            if constexpr (std::same_as<std::ranges::range_value_t<seq_t>, char>)
            {
                EXPECT_RANGE_EQ(rec.seq() | bio::views::char_strictly_to<bio::alphabet::dna5>, seqs[i]);
            }
            else
            {
                EXPECT_RANGE_EQ(rec.seq(), seqs[i]);
            }
            if constexpr (std::same_as<std::ranges::range_value_t<qual_t>, char>)
            {
                EXPECT_RANGE_EQ(rec.qual() | bio::views::char_strictly_to<bio::alphabet::phred42>, quals[i]);
            }
            else
            {
                EXPECT_RANGE_EQ(rec.qual(), quals[i]);
            }
        }
    }

    void do_read_test(std::string const & input)
    {
        /* containers */
        do_read_test_impl<std::string, std::string, std::string>(input);
        do_read_test_impl<std::string, std::vector<bio::alphabet::dna5>, std::vector<bio::alphabet::phred42>>(input);

        /* views */
        do_read_test_impl<std::string_view, std::string_view, std::string_view>(input);
        do_read_test_impl<std::string_view,
                          decltype(std::string_view{} | bio::views::char_strictly_to<bio::alphabet::dna5>),
                          decltype(std::string_view{} | bio::views::char_strictly_to<bio::alphabet::phred42>)>(input);
    }
};

// ----------------------------------------------------------------------------
// simple tests
// ----------------------------------------------------------------------------

TEST_F(read, simple)
{
    do_read_test(default_input);
}

TEST_F(read, no_trailing_newline)
{
    std::string input =
      R"raw(@ID1
ACGTTTTTTTTTTTTTTT
+
!##$%&'()*+,-./++-
@ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
@ID3 lala
ACGTTTA
+
!!!!!!!)raw";

    do_read_test(input);
}

TEST_F(read, double_id)
{
    std::string input =
      R"raw(@ID1
ACGTTTTTTTTTTTTTTT
+ID1
!##$%&'()*+,-./++-
@ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+ID2
!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
@ID3 lala
ACGTTTA
+ID3 lala
!!!!!!!)raw";

    do_read_test(input);
}

// ----------------------------------------------------------------------------
// custom tests
// ----------------------------------------------------------------------------

TEST_F(read, empty_seq)
{
    std::string const input =
      R"raw(@ID1

+

@ID2

+

@ID3 lala

+

)raw";

    std::istringstream                            istream{input};
    bio::io::format_input_handler<bio::io::fastq> input_handler{istream};
    default_rec_t                                 rec;

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID1"sv);
    EXPECT_TRUE(std::ranges::empty(rec.seq()));
    EXPECT_TRUE(std::ranges::empty(rec.qual()));

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID2"sv);
    EXPECT_TRUE(std::ranges::empty(rec.seq()));
    EXPECT_TRUE(std::ranges::empty(rec.qual()));

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID3 lala"sv);
    EXPECT_TRUE(std::ranges::empty(rec.seq()));
    EXPECT_TRUE(std::ranges::empty(rec.qual()));
}

struct options_t
{
    bool truncate_ids = false;
};

TEST_F(read, truncate_ids_off)
{
    std::istringstream                            istream{default_input};
    bio::io::format_input_handler<bio::io::fastq> input_handler{istream, options_t{}};
    default_rec_t                                 rec;

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID1"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[0]);
    EXPECT_RANGE_EQ(rec.qual(), quals[0]);

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID2"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[1]);
    EXPECT_RANGE_EQ(rec.qual(), quals[1]);

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID3 lala"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[2]);
    EXPECT_RANGE_EQ(rec.qual(), quals[2]);
}

TEST_F(read, truncate_ids_on)
{
    std::istringstream                            istream{default_input};
    bio::io::format_input_handler<bio::io::fastq> input_handler{istream, options_t{.truncate_ids = true}};
    default_rec_t                                 rec;

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID1"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[0]);
    EXPECT_RANGE_EQ(rec.qual(), quals[0]);

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID2"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[1]);
    EXPECT_RANGE_EQ(rec.qual(), quals[1]);

    input_handler.parse_next_record_into(rec);
    EXPECT_RANGE_EQ(rec.id(), "ID3"sv);
    EXPECT_RANGE_EQ(rec.seq(), seqs[2]);
    EXPECT_RANGE_EQ(rec.qual(), quals[2]);
}

// ----------------------------------------------------------------------------
// failure
// ----------------------------------------------------------------------------

TEST_F(read, fail_no_input)
{
    std::string const input{};

    std::istringstream istream{input};
    EXPECT_THROW(bio::io::format_input_handler<bio::io::fastq>{istream}, bio::io::file_open_error);
}

TEST_F(read, fail_no_id)
{
    std::string const input{"foo\nACGT"};

    std::istringstream                            istream{input};
    bio::io::format_input_handler<bio::io::fastq> input_handler{istream};
    default_rec_t                                 rec;

    EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::io::parse_error);
}

TEST_F(read, fail_no_plus)
{
    std::string const input{"@foo\nACGT\nbar"};

    std::istringstream                            istream{input};
    bio::io::format_input_handler<bio::io::fastq> input_handler{istream};
    default_rec_t                                 rec;

    EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::io::parse_error);
}

TEST_F(read, fail_size_mismatch)
{
    std::string const input =
      R"raw(@ID1
ACGTTTTTTTTTTTTTTT
+ID1
!##$%&'()*+,-./++
)raw";

    std::istringstream                            istream{input};
    bio::io::format_input_handler<bio::io::fastq> input_handler{istream};
    default_rec_t                                 rec;

    EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::io::parse_error);
}

TEST_F(read, fail_illegal_alphabet)
{
    std::string input{"@foo\nFOOBAR\n+\n!!!!!!\n"};

    std::istringstream                            istream{input};
    bio::io::format_input_handler<bio::io::fastq> input_handler{istream};
    using rec_t = bio::io::record<bio::meta::vtag_t<bio::io::field::id, bio::io::field::seq>,
                                  bio::meta::type_list<std::string_view, std::vector<bio::alphabet::dna5>>>;

    rec_t rec;

    EXPECT_THROW(input_handler.parse_next_record_into(rec), bio::alphabet::invalid_char_assignment);
}
