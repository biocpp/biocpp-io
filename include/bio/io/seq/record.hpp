// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::seq::record.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/alphabet/aminoacid/aa27.hpp>
#include <bio/alphabet/concept.hpp>
#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/quality/phred42.hpp>
#include <bio/meta/tag/ttag.hpp>
#include <bio/meta/type_list/traits.hpp>
#include <bio/ranges/concept.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/detail/concept.hpp>
#include <bio/io/detail/misc.hpp>
#include <bio/io/detail/range.hpp>
#include <bio/io/detail/tuple_record.hpp>
#include <bio/io/format/fasta.hpp>
#include <bio/io/format/fastq.hpp>
#include <bio/io/misc.hpp>
#include <bio/io/stream/transparent_istream.hpp>

namespace bio::io::seq
{

//-----------------------------------------------------------------------------
// the record
//-----------------------------------------------------------------------------

/*!\brief Record type for sequence I/O.
 * \ingroup seq
 * \tparam id_t     Type of the ID member. See the member for type requirements.
 * \tparam seq_t    Type of the SEQ member. See the member for type requirements.
 * \tparam qual_t   Type of the QUAL member. See the member for type requirements.
 *
 * \details
 *
 * This is the record template for sequence I/O. It encompasses three members.
 *
 * See the \ref record_faq for more information on record-based reading.
 *
 * ### Example 1
 *
 * Simple usage of the record in combination with a reader:
 *
 * \snippet test/snippet/seq/seq_reader.cpp simple_usage_file
 *
 * Note that the record type is hidden behind `auto & rec`.
 *
 * ### Example 2
 *
 * When creating record variables from existing data, the template
 * arguments can be deduced:
 *
 * \snippet test/snippet/seq/seq_record.cpp make_record
 *
 * To avoid copying existing data, you may want to use
 * bio::io::seq::tie_record() instead.
 */
template <typename _id_t   = std::string,
          typename _seq_t  = std::vector<alphabet::dna5>,
          typename _qual_t = std::vector<alphabet::phred42>>
struct record
{
    using id_t   = _id_t;   //!< The type of the ID member.
    using seq_t  = _seq_t;  //!< The type of the SEQ member.
    using qual_t = _qual_t; //!< The type of the QUAL member.

    /*!\brief The ID.
     * \details
     *
     * The type of this member can be changed to read data into a different structure (on input)
     * or write from a different type (on output).
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     * 1. any back-insertable range over the `char` alphabet, e.g. a std::string (**deep**)
     * 2. std::string_view (**shallow**)
     * 3. bio::meta::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over `char`.
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     *
     */
    id_t id{};

    /*!\brief The sequence.
     * \details
     *
     * The type of this member can be changed to read data into a different structure (on input)
     * or write from a different type (on output).
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over the `char` alphabet (**deep**)
     *   2. any back-insertable range over a bio::alphabet::alphabet (**deep**, converted elements)
     *   3. std::string_view (**shallow**)
     *   4. a std::ranges::transform_view defined on a std::string_view (**shallow**, converted elements)
     *   5. bio::meta::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over an bio::alphabet::alphabet.
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     *
     */
    seq_t seq{};

    /*!\brief The qualities.
     * \details
     *
     * The type of this member can be changed to read data into a different structure (on input)
     * or write from a different type (on output).
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over the `char` alphabet (**deep**)
     *   2. any back-insertable range over a bio::alphabet::quality (**deep**, converted elements)
     *   3. std::string_view (**shallow**)
     *   4. a std::ranges::transform_view defined on a std::string_view (**shallow**, converted elements)
     *   5. bio::meta::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over `char` or a bio::alphabet::quality.
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     */
    qual_t qual{};
};

//-----------------------------------------------------------------------------
// tie_record()
//-----------------------------------------------------------------------------

/*!\brief Helper function for easily creating a record of references.
 * \ingroup var
 * \tparam id_t         Type of the ID parameter.
 * \tparam seq_t        Type of the SEQ parameter.
 * \tparam qual_t       Type of the QUAL parameter.
 * \param[in,out] id    The ID parameter.
 * \param[in,out] seq   The SEQ parameter.
 * \param[in,out] qual  The QUAL parameter.
 * \returns A bio::io::seq::record where all the members are references to this function's parameters.
 * \details
 *
 * ### Example
 *
 * \snippet test/snippet/seq/seq_record.cpp tie_record
 */
template <typename id_t, typename seq_t, typename qual_t>
auto tie_record(id_t & id, seq_t & seq, qual_t & qual)

{
    return record<id_t &, seq_t &, qual_t &>{id, seq, qual};
}

//-----------------------------------------------------------------------------
// Aliases
//-----------------------------------------------------------------------------

/*!\name Pre-defined record aliases
 * \ingroup seq
 * \brief These can be used to configure the behaviour of the bio::io::seq::reader via
 * bio::io::seq::reader_options.
 * \{
 */
//!\brief Record type that can hold any kind of sequence (generic `char` alphabet).
//!\ingroup seq
using record_char_deep = record<std::string, std::string, std::string>;

//!\brief Record type that can hold any kind of sequence (generic `char` alphabet); shallow version.
//!\ingroup seq
using record_char_shallow = record<std::string_view, std::string_view, std::string_view>;

//!\brief Record type that reads DNA sequences (bio::alphabet::dna5)
// and corresponding qualities ( bio::alphabet::phred42).
//!\ingroup seq
using record_dna_deep = record<std::string, std::vector<alphabet::dna5>, std::vector<alphabet::phred42>>;

//!\brief Record type that reads DNA sequences (bio::alphabet::dna5)
// and corresponding qualities ( bio::alphabet::phred42); shallow version.
//!\ingroup seq
using record_dna_shallow = record<std::string_view,
                                  views::char_conversion_view_t<alphabet::dna5>,
                                  views::char_conversion_view_t<alphabet::phred42>>;

//!\brief Record type that reads Protein sequences (bio::alphabet::aa27) and ignores qualities.
//!\ingroup seq
using record_protein_deep = record<std::string, std::vector<alphabet::aa27>, meta::ignore_t>;

//!\brief Record type that reads Protein sequences (bio::alphabet::aa27) and ignores qualities; shallow version.
//!\ingroup seq
using record_protein_shallow = record<std::string_view, views::char_conversion_view_t<alphabet::aa27>, meta::ignore_t>;

//!\}

} // namespace bio::io::seq

//-----------------------------------------------------------------------------
// requirement checkers
//-----------------------------------------------------------------------------

namespace bio::io::seq::detail
{

//!\brief Validates the concepts that the record type needs to satisfy when being passed to a reader.
template <typename id_t, typename seq_t, typename qual_t>
constexpr bool record_read_concept_checker(std::type_identity<seq::record<id_t, seq_t, qual_t>>)
{
    // TODO(GCC11): once GCC10 is dropped, remove the "<typename t = seq_t>"
    static_assert(ranges::back_insertable_with<id_t, char> ||
                    meta::one_of<id_t, std::string_view, meta::ignore_t, meta::ignore_t const>,
                  "Requirements for the type of the ID-field not met. See documentation for bio::io::seq::record.");
    static_assert(io::detail::lazy_concept_checker([]<typename t = seq_t>(auto) requires(
                    meta::one_of<t, std::string_view, meta::ignore_t, meta::ignore_t const> ||
                    (ranges::back_insertable<t> && alphabet::alphabet<std::ranges::range_reference_t<t>>) ||
                    io::detail::transform_view_on_string_view<t>) { return std::true_type{}; }),
                  "Requirements for the type of the SEQ-field not met. See documentation for bio::io::seq::record.");
    static_assert(io::detail::lazy_concept_checker([]<typename t = qual_t>(auto) requires(
                    meta::one_of<t, std::string_view, meta::ignore_t, meta::ignore_t const> ||
                    (ranges::back_insertable<t> && alphabet::alphabet<std::ranges::range_reference_t<t>>) ||
                    io::detail::transform_view_on_string_view<t>) { return std::true_type{}; }),
                  "Requirements for the type of the QUAL-field not met. See documentation for "
                  "bio::io::seq::record.");
    return true;
}

//!\brief Validates the concepts that the record type needs to satisfy when being passed to a reader.
template <typename id_t, typename seq_t, typename qual_t>
constexpr bool record_write_concept_checker(std::type_identity<seq::record<id_t, seq_t, qual_t>>)
{
    // TODO(GCC11): once GCC10 is dropped, remove the "<typename t = seq_t>"
    static_assert(io::detail::char_range<id_t>,
                  "Requirements for the type of the ID-field not met. See documentation for bio::io::seq::record.");
    static_assert(io::detail::lazy_concept_checker([]<typename t = seq_t>(auto) requires(
                    std::ranges::forward_range<t> && alphabet::alphabet<std::ranges::range_reference_t<t>>) {
                      return std::true_type{};
                  }),
                  "Requirements for the type of the SEQ-field not met. See documentation for bio::io::seq::record.");
    static_assert(io::detail::lazy_concept_checker([]<typename t = qual_t>(auto) requires(
                    std::ranges::forward_range<t> && (alphabet::quality<std::ranges::range_reference_t<t>> ||
                                                      std::same_as<char, std::ranges::range_value_t<t>>)) {
                      return std::true_type{};
                  }),
                  "Requirements for the type of the QUAL-field not met. See documentation for "
                  "bio::io::seq::record.");
    return true;
}

//-----------------------------------------------------------------------------
// format_handler_mixin
//-----------------------------------------------------------------------------

//!\brief The field_ids used in this domain.
//!\ingroup seq
static constexpr auto field_ids = meta::vtag<io::detail::field::id, io::detail::field::seq, io::detail::field::qual>;

} // namespace bio::io::seq::detail

namespace bio::io::seq
{

template <typename id_t, typename seq_t, typename qual_t>
struct record;

//!\brief Mixin for format handlers that helps converts to the tuple-record.
//!\ingroup seq
struct format_handler_mixin
{
    //!\brief Convert from domain-specific record to tuple-record.
    template <typename... arg_ts>
    static auto record2tuple_record(seq::record<arg_ts...> & in_record)
    {
        return io::detail::tie_tuple_record(detail::field_ids, in_record.id, in_record.seq, in_record.qual);
    }
};

} // namespace bio::io::seq
