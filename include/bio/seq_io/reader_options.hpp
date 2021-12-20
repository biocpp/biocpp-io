// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::seq_io::reader_options and various field_types definitions.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/misc.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/views/to.hpp>

#include <bio/detail/concept.hpp>
#include <bio/detail/misc.hpp>
#include <bio/detail/range.hpp>
#include <bio/format/fasta.hpp>
#include <bio/misc.hpp>
#include <bio/seq_io/misc.hpp>
#include <bio/stream/transparent_istream.hpp>

// TODO replace seqan3::views::char_to with seqan3::views::char_strictly_to
namespace bio::seq_io
{

/*!\addtogroup seq_io
 * \{
 */

/*!\name Pre-defined field types
 * \brief These can be used to configure the behaviour of the bio::seq_io::reader via bio::seq_io::reader_options.
 * \{
 */

/*!\brief The generic field types template.
 * \tparam ownership Return shallow or deep types.
 * \details
 *
 * You can use this to configure a deep record or one with custom alphabets.
 *
 * ### Example
 *
 * \snippet test/snippet/seq_io/seq_io_reader_options.cpp example_simple
 *
 * Type of the ID will be std::string, type of the sequence will be std::vector<seqan3::dna4> and
 * type of the qualities will be std::vector<seqan3::phred42>.
 */
template <bio::ownership   ownership   = bio::ownership::shallow,
          seqan3::alphabet seq_alph_t  = seqan3::dna5,
          seqan3::alphabet qual_alph_t = seqan3::phred63>
inline constinit auto field_types = []()
{
    if constexpr (ownership == bio::ownership::deep)
    {
        return ttag<std::string,
                    std::conditional_t<std::same_as<seq_alph_t, char>, std::string, std::vector<seq_alph_t>>,
                    std::conditional_t<std::same_as<qual_alph_t, char>, std::string, std::vector<qual_alph_t>>>;
    }
    else
    {
        return ttag<std::string_view,
                    std::conditional_t<std::same_as<seq_alph_t, char>,
                                       std::string_view,
                                       decltype(std::string_view{} | seqan3::views::char_to<seq_alph_t>)>,
                    std::conditional_t<std::same_as<qual_alph_t, char>,
                                       std::string_view,
                                       decltype(std::string_view{} | seqan3::views::char_to<qual_alph_t>)>>;
    }
}();

/*!\brief The field types for reading DNA data.
 * \details
 *
 * This is the default.
 *
 * Configures a shallow record where sequence data is seqan3::dna5 and quality data is seqan3::phred63.
 */
inline constinit auto field_types_dna = field_types<>;

/*!\brief The field types for reading protein data.
 * \tparam ownership Return shallow or deep types.
 * \details
 *
 * Configures a shallow record where sequence data is seqan3::aa27 and quality data is seqan3::phred63.
 */
inline constinit auto field_types_protein = field_types<ownership::shallow, seqan3::aa27>;

/*!\brief The field types for reading any data.
 * \details
 *
 * Configures a shallow record where sequence and quality data are plain characters.
 */
inline constinit auto field_types_char = field_types<ownership::shallow, char, char>;

/*!\brief The field types for raw I/O.
 * \details
 *
 * Every field is configured as a std::span of std::byte (this enables "raw" io).
 *
 * ATTENTION: The exact content of this byte-span depends on the format and is likely not
 * compatible between formats. Use at your own risk!
 */
inline constinit auto field_types_raw =
  ttag<std::span<std::byte const>, std::span<std::byte const>, std::span<std::byte const>>;
// TODO use seqan3::list_traits::repeat as soon as available

//!\}
//!\}

/*!\brief Options that can be used to configure the behaviour of seqan3::seq_io::reader.
 * \tparam field_ids_t   Type of the field_ids member (usually deduced).
 * \tparam field_types_t Type of the field_types member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup seq_io
 * \details
 *
 * By default, the reader options assume DNA data. You can select bio::seq_io::field_types_protein to
 * read protein data or bio::seq_io::field_types_char to store in an agnostic type.
 *
 * ### Example
 *
 * The reader options can be easily set via [designated
 * initialisers](https://en.cppreference.com/w/cpp/language/aggregate_initialization).
 *
 * To switch from DNA reading (the default) to protein reading and activate truncating of IDs, do
 * the following:
 *
 * \snippet test/snippet/seq_io/seq_io_reader_options.cpp example_simple
 *
 * It is not required to specify all options; default values are documented.
 *
 * Please be aware that those options that you modify need to be set in the correct order -- **which is
 * alphabetical order** for all option classes in this library.
 *
 * Typically, the options are set as part of the bio::seq_io::reader construction. See the respective
 * documentation page for more information.
 *
 * ### Example (advanced)
 *
 * This code switches from seqan3::dna5 to seqan3::dna4 alphabet, from seqan3::phred63 to
 * seqan3::phred42, and reduces the amount of threads used:
 *
 * \snippet test/snippet/seq_io/seq_io_reader_options.cpp example_advanced1
 *
 * This code makes FASTA the only legal format and creates records with only the sequence field as
 * a std::string:
 *
 * \snippet test/snippet/seq_io/seq_io_reader_options.cpp example_advanced2
 *
 * ### Field type specific restrictions
 *
 * This section is only relevant if you specify the #field_types member manually, i.e. if you
 * change the field_types but do not use one of the predefined tags.
 *
 * 1. bio::field::id
 *   * any back-insertable range over the `char` alphabet (copy of elements returned)
 *   * std::string_view (view into a buffer is returned)
 * 2. bio::field::seq
 *   * any back-insertable range over the `char` alphabet (copy of elements returned)
 *   * any back-insertable range over a seqan3::alphabet (elements are transformed via seqan3::views::char_strictly_to)
 *   * std::string_view (view into a buffer is returned)
 *   * a std::ranges::transform_view defined on a std::string_view (transformation view is returned)
 * 3. bio::field::qual
 *   * any back-insertable range over the `char` alphabet (copy of elements returned)
 *   * any back-insertable range over a seqan3::alphabet (elements are transformed via seqan3::views::char_strictly_to)
 *   * std::string_view (view into a buffer is returned)
 *   * a std::ranges::transform_view defined on a std::string_view (transformation view is returned)
 */
template <typename field_ids_t   = decltype(default_field_ids),
          typename field_types_t = decltype(field_types_dna),
          typename formats_t     = seqan3::type_list<fasta>>
struct reader_options
{
    /*!\brief The fields that shall be contained in each record; a bio::vtag over bio::field.
     * \details
     *
     * It is usually not necessary to change this.
     * **If you do, you need to adapt field_types, as well!**
     */
    field_ids_t field_ids = default_field_ids;

    /*!\brief The types corresponding to each field; a bio::ttag over the types.
     *
     * \details
     *
     * See seqan3::seq_io::reader for an overview of the supported field/type combinations.
     */
    field_types_t field_types = field_types_dna;

    /*!\brief The formats that input files can take; a bio::ttag over the types.
     *
     * \details
     *
     * See seqan3::seq_io::reader for an overview of the the supported formats.
     */
    formats_t formats = ttag<fasta>;

    //!\brief Options that are passed on to the internal stream oject.
    transparent_istream_options stream_options{};

    //!\brief Truncate IDs at first whitespace.
    bool truncate_ids = false;

private:
    static_assert(detail::is_fields_tag<field_ids_t>, "field_ids must be a bio::vtag over bio::field.");

    static_assert(detail::is_type_list<field_types_t>, "field_types must be a bio::ttag / seqan3::type_list.");

    static_assert(detail::is_type_list<formats_t>, "formats must be a bio::ttag / seqan3::type_list.");

    static_assert(field_ids_t::size == seqan3::list_traits::size<field_types_t>,
                  "field_ids and field_types must have the same size.");

    //!\brief Type of the record.
    using record_t = record<field_ids_t, field_types_t>;

    static_assert(
      detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
        !field_ids_t::contains(field::id) || detail::back_insertable_with<record_element_t<field::id, rec_t>, char> ||
        std::same_as<std::string_view, record_element_t<field::id, rec_t>>) { return std::true_type{}; }),
      "Requirements for the field-type of the ID-field not met. See documentation for bio::seq_io::reader_options.");

    static_assert(
      detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
        !field_ids_t::contains(field::seq) ||
        (detail::back_insertable<record_element_t<field::seq, rec_t>> &&
         seqan3::alphabet<std::ranges::range_reference_t<record_element_t<field::seq, rec_t>>>) ||
        std::same_as<std::string_view, record_element_t<field::seq, rec_t>> ||
        detail::transform_view_on_string_view<record_element_t<field::seq, rec_t>>) { return std::true_type{}; }),
      "Requirements for the field-type of the SEQ-field not met. See documentation for bio::seq_io::reader_options.");

    static_assert(
      detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
        !field_ids_t::contains(field::qual) ||
        (detail::back_insertable<record_element_t<field::qual, rec_t>> &&
         seqan3::alphabet<std::ranges::range_reference_t<record_element_t<field::qual, rec_t>>>) ||
        std::same_as<std::string_view, record_element_t<field::qual, rec_t>> ||
        detail::transform_view_on_string_view<record_element_t<field::qual, rec_t>>) { return std::true_type{}; }),
      "Requirements for the field-type of the QUAL-field not met. See documentation for bio::seq_io::reader_options.");
};

} // namespace bio::seq_io
