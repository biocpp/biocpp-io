// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <sstream>

#include <gtest/gtest.h>

#include <bio/alphabet/nucleotide/dna4.hpp>
#include <bio/test/expect_range_eq.hpp>
#include <bio/test/expect_same_type.hpp>

#include <bio/io/detail/tuple_record.hpp>
#include <bio/io/misc.hpp>

using namespace bio::alphabet::literals;

using default_fields =
  bio::meta::vtag_t<bio::io::detail::field::seq, bio::io::detail::field::id, bio::io::detail::field::qual>;

// ----------------------------------------------------------------------------
// fields
// ----------------------------------------------------------------------------

TEST(fields, usage)
{
    EXPECT_TRUE(default_fields::contains(bio::io::detail::field::seq));
    EXPECT_TRUE(default_fields::contains(bio::io::detail::field::id));
    EXPECT_TRUE(default_fields::contains(bio::io::detail::field::qual));
    EXPECT_EQ(default_fields::index_of(bio::io::detail::field::seq), 0ul);
    EXPECT_EQ(default_fields::index_of(bio::io::detail::field::id), 1ul);
    EXPECT_EQ(default_fields::index_of(bio::io::detail::field::qual), 2ul);
}

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

struct record : public ::testing::Test
{
    using ids = bio::meta::vtag_t<bio::io::detail::field::id, bio::io::detail::field::seq>;
    using record_type =
      bio::io::detail::tuple_record<ids, bio::meta::type_list<std::string, std::vector<bio::alphabet::dna4>>>;
};

TEST_F(record, definition_tuple_traits)
{
    EXPECT_TRUE(
      (std::is_same_v<typename record_type::base_type, std::tuple<std::string, std::vector<bio::alphabet::dna4>>>));

    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, record_type>, std::string>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, record_type>, std::vector<bio::alphabet::dna4>>));
    EXPECT_EQ(std::tuple_size_v<record_type>, 2ul);

    // TODO(bio): reactivate once in biocpp_core
    //      EXPECT_TRUE(bio::meta::tuple_like<record_type>);
}

TEST_F(record, record_element)
{
    EXPECT_TRUE(
      (std::is_same_v<bio::io::detail::tuple_record_element_t<bio::io::detail::field::id, record_type>, std::string>));
    EXPECT_TRUE((std::is_same_v<bio::io::detail::tuple_record_element_t<bio::io::detail::field::seq, record_type>,
                                std::vector<bio::alphabet::dna4>>));
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
    EXPECT_RANGE_EQ(std::get<std::vector<bio::alphabet::dna4>>(r), "ACGT"_dna4);
}

TEST_F(record, get_by_field)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(bio::io::detail::get<bio::io::detail::field::id>(r), "MY ID");
    EXPECT_RANGE_EQ(bio::io::detail::get<bio::io::detail::field::seq>(r), "ACGT"_dna4);
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

    auto r =
      bio::io::detail::make_tuple_record(bio::meta::vtag<bio::io::detail::field::id, bio::io::detail::field::seq>,
                                         s,
                                         vec);
    EXPECT_TRUE((std::same_as<decltype(r), record::record_type>));
}

TEST_F(record, tie_record)
{
    std::string s   = "MY ID";
    auto        vec = "ACGT"_dna4;

    auto r = bio::io::detail::tie_tuple_record(bio::meta::vtag<bio::io::detail::field::id, bio::io::detail::field::seq>,
                                               s,
                                               vec);
    EXPECT_TRUE(
      (std::same_as<
        decltype(r),
        bio::io::detail::tuple_record<record::ids,
                                      bio::meta::type_list<std::string &, std::vector<bio::alphabet::dna4> &>>>));
}
