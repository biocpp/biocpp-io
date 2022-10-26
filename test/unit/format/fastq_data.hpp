// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <ranges>
#include <string>

#include <bio/alphabet/custom/char.hpp>
#include <bio/ranges/to.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/detail/magic_get.hpp>
#include <bio/io/seq/record.hpp>

struct writer_options
{
    size_t max_seq_line_length = 70;
    bool   double_id           = false;
    bool   windows_eol         = false;
};

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

inline std::vector<std::string_view> ids   = {"ID1", "ID2", "ID3 lala"};
inline std::vector<std::string_view> seqs  = {"ACGTTTTTTTTTTTTTTT",
                                              "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                               "TTTTTTTTT",
                                              "ACGTTTA"};
inline std::vector<std::string_view> quals = {"!##$%&'()*+,-./++-",
                                              "!##$&'()*+,-./"
                                              "+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE",
                                              "!!!!!!!"};

//=============================================================================
// records
//=============================================================================

/* auxiliary stuff */
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
