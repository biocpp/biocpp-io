// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::seq_io::record.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/alphabet/adaptation/char.hpp>
#include <bio/alphabet/aminoacid/aa27.hpp>
#include <bio/alphabet/concept.hpp>
#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/quality/phred42.hpp>
#include <bio/meta/tag/ttag.hpp>
#include <bio/meta/type_list/traits.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/detail/concept.hpp>
#include <bio/io/detail/misc.hpp>
#include <bio/io/detail/range.hpp>
#include <bio/io/detail/tuple_record.hpp>
#include <bio/io/format/fasta.hpp>
#include <bio/io/format/fastq.hpp>
#include <bio/io/misc.hpp>
#include <bio/io/stream/transparent_istream.hpp>

namespace bio::io::seq_io
{

/*!\brief Record type for sequence I/O.
 * \ingroup seq_io
 * \tparam id_t     Type of the ID member. See the member for type requirements.
 * \tparam seq_t    Type of the SEQ member. See the member for type requirements.
 * \tparam qual_t   Type of the QUAL member. See the member for type requirements.
 *
 * \details
 *
 * This is the record template for sequence I/O. It encompasses three members.
 *
 * See the \ref record_faq for more information on record-based reading.
 */
template <typename id_t   = std::string,
          typename seq_t  = std::vector<alphabet::dna5>,
          typename qual_t = std::vector<alphabet::phred42>>
struct record
{
    /*!\brief The ID.
     * \details
     *
     * The type of this member can be changed to read data into a different structure (on input)
     * or write from a different type (on output).
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
     *
     * 1. any back-insertable range over the `char` alphabet, e.g. a std::string (**deep**)
     * 2. std::string_view (**shallow**)
     * 3. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over an alphabet that is convertible to `char`.
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
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over the `char` alphabet (**deep**)
     *   2. any back-insertable range over a bio::alphabet::alphabet (**deep**, converted elements)
     *   3. std::string_view (**shallow**)
     *   4. a std::ranges::transform_view defined on a std::string_view (**shallow**, converted elements)
     *   5. bio::meta::ignore_t (**ignored**)
     *
     * ### Type requirements when writing
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
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
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over the `char` alphabet (**deep**)
     *   2. any back-insertable range over a bio::alphabet::quality (**deep**, converted elements)
     *   3. std::string_view (**shallow**)
     *   4. a std::ranges::transform_view defined on a std::string_view (**shallow**, converted elements)
     *   5. bio::meta::ignore_t (**ignored**)
     *
     * ### Type requirements when writing
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over an bio::alphabet::alphabet.
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     */
    qual_t qual{};
};

/*!\name Pre-defined record aliases
 * \ingroup seq_io
 * \brief These can be used to configure the behaviour of the bio::io::seq_io::reader via
 * bio::io::seq_io::reader_options.
 * \{
 */
//!\brief Record type that can hold any kind of sequence (generic `char` alphabet).
//!\ingroup seq_io
using record_char = record<std::string, std::string, std::string>;

//!\brief Record type that can hold any kind of sequence (generic `char` alphabet); shallow version.
//!\ingroup seq_io
using record_char_shallow = record<std::string_view, std::string_view, std::string_view>;

//!\brief Record type that reads DNA sequences (bio::alphabet::dna5)
// and corresponding qualities ( bio::alphabet::phred42).
//!\ingroup seq_io
using record_dna = record<std::string, std::vector<alphabet::dna5>, std::vector<alphabet::phred42>>;

//!\brief Record type that reads DNA sequences (bio::alphabet::dna5)
// and corresponding qualities ( bio::alphabet::phred42); shallow version.
//!\ingroup seq_io
using record_dna_shallow =
  record<std::string_view, conversion_view_t<alphabet::dna5>, conversion_view_t<alphabet::phred42>>;

//!\brief Record type that reads DNA sequences (bio::alphabet::aa27) and ignores qualities.
//!\ingroup seq_io
using record_protein = record<std::string, std::vector<alphabet::aa27>, ignore_t>;

//!\brief Record type that reads DNA sequences (bio::alphabet::aa27) and ignores qualities; shallow version.
//!\ingroup seq_io
using record_protein_shallow = record<std::string_view, conversion_view_t<alphabet::aa27>, ignore_t>;

//!\}

} // namespace bio::io::seq_io

namespace bio::io::detail // TODO move this to seq_io::detail?
{

//!\brief Validates the concepts that the record type needs to satisfy when being passed to a reader.
template <typename id_t, typename seq_t, typename qual_t>
constexpr bool record_read_concept_checker(std::type_identity<seq_io::record<id_t, seq_t, qual_t>>)
{
    // TODO(GCC11): once GCC10 is dropped, remove the "<typename t = seq_t>"
    static_assert(io::detail::lazy_concept_checker([]<typename t = id_t>(auto) requires(
                    io::detail::back_insertable_with<t, char> ||
                    io::detail::one_of<t, std::string_view, ignore_t, ignore_t const>) { return std::true_type{}; }),
                  "Requirements for the type of the ID-field not met. See documentation for bio::io::seq_io::record.");
    static_assert(io::detail::lazy_concept_checker([]<typename t = seq_t>(auto) requires(
                    io::detail::one_of<t, std::string_view, ignore_t, ignore_t const> ||
                    (io::detail::back_insertable<t> && alphabet::alphabet<std::ranges::range_reference_t<t>>) ||
                    io::detail::transform_view_on_string_view<t>) { return std::true_type{}; }),
                  "Requirements for the type of the SEQ-field not met. See documentation for bio::io::seq_io::record.");
    static_assert(
      io::detail::lazy_concept_checker([]<typename t = qual_t>(auto) requires(
        io::detail::one_of<t, std::string_view, ignore_t, ignore_t const> ||
        (io::detail::back_insertable<t> && alphabet::alphabet<std::ranges::range_reference_t<t>>) ||
        io::detail::transform_view_on_string_view<t>) { return std::true_type{}; }),
      "Requirements for the type of the QUAL-field not met. See documentation for bio::io::seq_io::record.");
    return true;
}

} // namespace bio::io::detail

namespace bio::io::seq_io
{

//!\brief The field_ids used in this domain.
//!\ingroup seq_io
static constexpr auto field_ids = meta::vtag<detail::field::id, detail::field::seq, detail::field::qual>;

template <typename id_t, typename seq_t, typename qual_t>
struct record;

//!\brief Mixin for format handlers that helps converts to the tuple-record.
//!\ingroup seq_io
struct format_handler_mixin
{
    //!\brief Convert from domain-specific record to tuple-record.
    template <typename... arg_ts>
    static auto record2tuple_record(seq_io::record<arg_ts...> & in_record)
    {
        return io::detail::tie_tuple_record(field_ids, in_record.id, in_record.seq, in_record.qual);
    }
};

} // namespace bio::io::seq_io
