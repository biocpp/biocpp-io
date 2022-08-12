// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/algorithm>
#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/tuple/concept.hpp>

#include <bio/io/misc.hpp>
#include <bio/io/record.hpp>

using seqan3::operator""_dna4;

using default_fields = bio::io::vtag_t<bio::io::field::seq, bio::io::field::id, bio::io::field::qual>;

// This is needed for EXPECT_RANGE_EQ:
namespace seqan3
{
template <typename char_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & stream, bio::io::field f)
{
    stream << "<field: " << static_cast<size_t>(f) << ">";
    return stream;
}
} // namespace seqan3

// ----------------------------------------------------------------------------
// fields
// ----------------------------------------------------------------------------

TEST(fields, usage)
{
    EXPECT_TRUE(default_fields::contains(bio::io::field::seq));
    EXPECT_TRUE(default_fields::contains(bio::io::field::id));
    EXPECT_TRUE(default_fields::contains(bio::io::field::qual));
    EXPECT_EQ(default_fields::index_of(bio::io::field::seq), 0ul);
    EXPECT_EQ(default_fields::index_of(bio::io::field::id), 1ul);
    EXPECT_EQ(default_fields::index_of(bio::io::field::qual), 2ul);
}

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

struct record : public ::testing::Test
{
    using ids         = bio::io::vtag_t<bio::io::field::id, bio::io::field::seq>;
    using record_type = bio::io::record<ids, seqan3::type_list<std::string, seqan3::dna4_vector>>;
};

TEST_F(record, definition_tuple_traits)
{
    EXPECT_TRUE((std::is_same_v<typename record_type::base_type, std::tuple<std::string, seqan3::dna4_vector>>));

    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, record_type>, std::string>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, record_type>, seqan3::dna4_vector>));
    EXPECT_EQ(std::tuple_size_v<record_type>, 2ul);

    EXPECT_TRUE(seqan3::tuple_like<record_type>);
}

TEST_F(record, record_element)
{
    EXPECT_TRUE((std::is_same_v<bio::io::record_element_t<bio::io::field::id, record_type>, std::string>));
    EXPECT_TRUE((std::is_same_v<bio::io::record_element_t<bio::io::field::seq, record_type>, seqan3::dna4_vector>));
}

TEST_F(record, construction)
{
    [[maybe_unused]] record_type r{"MY ID", "ACGT"_dna4};
}

TEST_F(record, get_by_index)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(std::get<0>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<1>(r), "ACGT"_dna4);
}

TEST_F(record, get_by_type)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(std::get<std::string>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<seqan3::dna4_vector>(r), "ACGT"_dna4);
}

TEST_F(record, get_by_field)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(bio::io::get<bio::io::field::id>(r), "MY ID");
    EXPECT_RANGE_EQ(bio::io::get<bio::io::field::seq>(r), "ACGT"_dna4);
}

TEST_F(record, get_by_member)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(r.id(), "MY ID");
    EXPECT_RANGE_EQ(r.seq(), "ACGT"_dna4);
}

TEST_F(record, make_record)
{
    std::string s   = "MY ID";
    auto        vec = "ACGT"_dna4;

    auto r = bio::io::make_record(bio::io::vtag<bio::io::field::id, bio::io::field::seq>, s, vec);
    EXPECT_TRUE((std::same_as<decltype(r), record::record_type>));
}

TEST_F(record, tie_record)
{
    std::string s   = "MY ID";
    auto        vec = "ACGT"_dna4;

    auto r = bio::io::tie_record(bio::io::vtag<bio::io::field::id, bio::io::field::seq>, s, vec);
    EXPECT_TRUE(
      (std::same_as<decltype(r),
                    bio::io::record<record::ids, seqan3::type_list<std::string &, std::vector<seqan3::dna4> &>>>));
}
