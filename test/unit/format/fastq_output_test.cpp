// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2020-2022, deCODE Genetics
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <bio/io/format/fastq_output_handler.hpp>

#include "seq_output_detail.hpp"

inline std::string_view fastq_default_output = R"(@ID1
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
)";

inline std::string_view fastq_default_output_double_id = R"(@ID1
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
!!!!!!!
)";

template <typename T>
class fastq_write : public ::testing::Test
{};

TYPED_TEST_SUITE(fastq_write, numbers, );

TYPED_TEST(fastq_write, do_it)
{
    using salph_t       = bio::meta::list_traits::at<TypeParam::value, salphs>;
    using qalph_t       = bio::meta::list_traits::at<TypeParam::value, qalphs>;
    constexpr bool deep = TypeParam::value % 2;

    std::string s = do_test<bio::io::fastq, deep, salph_t, qalph_t>(writer_options{});

    EXPECT_EQ(s, fastq_default_output);
}

TYPED_TEST(fastq_write, double_id)
{
    using salph_t       = bio::meta::list_traits::at<TypeParam::value, salphs>;
    using qalph_t       = bio::meta::list_traits::at<TypeParam::value, qalphs>;
    constexpr bool deep = TypeParam::value % 2;

    std::string s = do_test<bio::io::fastq, deep, salph_t, qalph_t>(writer_options{.double_id = true});

    EXPECT_EQ(s, fastq_default_output_double_id);
}

// ----------------------------------------------------------------------------
// failure
// ----------------------------------------------------------------------------

TEST(fastq_write_manual, size_mismatch)
{
    std::ostringstream                             ostr{};
    bio::io::format_output_handler<bio::io::fastq> handler{ostr};
    auto recs = example_records<1, bio::alphabet::dna5, bio::alphabet::phred42>();
    recs[0].qual.resize(2);
    EXPECT_THROW(handler.write_record(recs[0]), bio::io::format_error);
}
