// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2020-2022, deCODE Genetics
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <bio/io/format/fasta_output_handler.hpp>

#include "seq_output_detail.hpp"

inline std::string_view fasta_default_output = R"(>ID1
ACGTTTTTTTTTTTTTTT

>ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTT

>ID3 lala
ACGTTTA
)";

inline std::string_view fasta_default_output_no_linebreak = R"(>ID1
ACGTTTTTTTTTTTTTTT

>ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

>ID3 lala
ACGTTTA
)";

template <typename T>
class fasta_write : public ::testing::Test
{};

TYPED_TEST_SUITE(fasta_write, numbers, );

TYPED_TEST(fasta_write, do_it)
{
    using salph_t       = bio::meta::list_traits::at<TypeParam::first_value, salphs>;
    using qalph_t       = bio::meta::list_traits::at<TypeParam::first_value, qalphs>;
    constexpr bool deep = TypeParam::first_value % 2;

    EXPECT_EQ(writer_options{}.max_seq_line_length, 70ull);

    std::string s = do_test<bio::io::fasta, deep, salph_t, qalph_t>(writer_options{});

    EXPECT_EQ(s, fasta_default_output);
}

TYPED_TEST(fasta_write, no_linebreak)
{
    using salph_t       = bio::meta::list_traits::at<TypeParam::first_value, salphs>;
    using qalph_t       = bio::meta::list_traits::at<TypeParam::first_value, qalphs>;
    constexpr bool deep = TypeParam::first_value % 2;

    std::string s = do_test<bio::io::fasta, deep, salph_t, qalph_t>(writer_options{.max_seq_line_length = 0});

    EXPECT_EQ(s, fasta_default_output_no_linebreak);
}
