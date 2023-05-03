// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <ranges>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include <bio/alphabet/custom/char.hpp>
#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/quality/phred42.hpp>
#include <bio/meta/tag/vtag.hpp>
#include <bio/ranges/to.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>
#include <bio/test/expect_range_eq.hpp>
#include <bio/test/expect_same_type.hpp>

#include <bio/io/format/format_output_handler.hpp>
#include <bio/io/seq/record.hpp>

using namespace bio::alphabet::literals;
using namespace bio::meta::literals;
using std::literals::string_view_literals::operator""sv;

struct writer_options
{
    size_t max_seq_line_length = 70;
    bool   double_id           = false;
    bool   windows_eol         = false;
};

inline std::vector<std::string_view> ids   = {"ID1", "ID2", "ID3 lala"};
inline std::vector<std::string_view> seqs  = {"ACGTTTTTTTTTTTTTTT",
                                              "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                               "TTTTTTTTTT",
                                              "ACGTTTA"};
inline std::vector<std::string_view> quals = {"!##$%&'()*+,-./++-",
                                              "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBD"
                                              "EBDEBDEBDE",
                                              "!!!!!!!"};

template <typename t>
constexpr bool operator==(bio::views::char_conversion_view_t<t> const & lhs,
                          bio::views::char_conversion_view_t<t> const & rhs)
{
    return std::ranges::equal(lhs, rhs);
}

template <bool deep, typename salph_t, typename qalph_t>
auto example_records()
{
    if constexpr (deep)
    {
        using rec_t = bio::io::seq::record<std::string, std::vector<salph_t>, std::vector<qalph_t>>;
        std::vector<rec_t> recs;
        recs.resize(3);

        for (size_t i = 0; i < 3; ++i)
        {
            recs[i] = rec_t{static_cast<std::string>(ids[i]),
                            seqs[i] | bio::views::char_strictly_to<salph_t> | bio::ranges::to<std::vector>(),
                            quals[i] | bio::views::char_strictly_to<qalph_t> | bio::ranges::to<std::vector>()};
        }

        return recs;
    }
    else
    {
        using rec_t = bio::io::seq::record<std::string_view,
                                           bio::views::char_conversion_view_t<salph_t>,
                                           bio::views::char_conversion_view_t<qalph_t>>;
        std::vector<rec_t> recs;
        recs.resize(3);

        for (size_t i = 0; i < 3; ++i)
        {
            recs[i] = rec_t{ids[i],
                            seqs[i] | bio::views::char_strictly_to<salph_t>,
                            quals[i] | bio::views::char_strictly_to<qalph_t>};
        }

        return recs;
    }
}

template <typename format_t, bool deep, typename salph_t, typename qalph_t>
std::string do_test(writer_options opt)
{
    std::ostringstream ostr{};

    {
        bio::io::format_output_handler<format_t> handler{ostr, opt};

        auto recs = example_records<deep, salph_t, qalph_t>();

        for (auto & rec : recs)
            handler.write_record(rec);
    }

    return ostr.str();
}

using salphs  = bio::meta::type_list<char, bio::alphabet::dna5>;
using qalphs  = bio::meta::type_list<char, bio::alphabet::phred42>;
using numbers = ::testing::Types<bio::meta::vtag_t<0>, bio::meta::vtag_t<1>>;
