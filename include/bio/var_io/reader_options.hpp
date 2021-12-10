// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::var_io::reader_options and various pre-defined field_types.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/utility/type_list/traits.hpp>

#include <bio/format/bcf.hpp>
#include <bio/format/vcf.hpp>
#include <bio/stream/transparent_istream.hpp>
#include <bio/stream/transparent_ostream.hpp>
#include <bio/var_io/dynamic_type.hpp>
#include <bio/var_io/header.hpp>
#include <bio/var_io/misc.hpp>

/* TODO
 *
 * field_types (default) with shallow/deep switch and all VCF except GENOTYPES
 */

namespace bio::detail
{

template <typename t>
concept info_element_concept = detail::decomposable_into_two<t> &&
  (detail::char_range<detail::first_elem_t<t>> ||
   std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::is_dynamic_type<detail::second_elem_t<t>>;

template <typename t>
concept genotype_bcf_style_concept = detail::decomposable_into_two<t> &&
  (detail::char_range<detail::first_elem_t<t>> ||
   std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::is_dynamic_vector_type<detail::second_elem_t<t>>;

template <typename t>
concept genotypes_vcf_style_concept =
  detail::decomposable_into_two<t> && detail::back_insertable<detail::first_elem_t<t>> &&
  detail::char_range<std::ranges::range_reference_t<detail::first_elem_t<t>>> &&
  detail::vector_like<detail::second_elem_t<t>> &&
  detail::vector_like<std::ranges::range_reference_t<detail::second_elem_t<t>>> &&
  detail::is_dynamic_type<std::ranges::range_value_t<std::ranges::range_reference_t<detail::second_elem_t<t>>>>;

} // namespace bio::detail

namespace bio::var_io
{

/*!\name Pre-defined field types
 * \brief These can be used to configure the behaviour of the bio::var_io::reader via bio::var_io::reader_options.
 * \{
 */

/*!\brief The default field types for variant io.
 *!\ingroup var_io
 *
 * \details
 *
 * These traits define a record type with minimal memory allocations for all input formats.
 * It is the recommended record type when iterating ("streaming") over files that ca be any variant IO format.
 *
 * The "style" of the record resembles the BCF specification, i.e. contigs, FILTERs and INFO identifiers are
 * represented as numbers (not strings); and the genotypes are encoded by-genotype (not by-sample).
 * See bio::var_io::genotypes_bcf_style for more information on the latter.
 *
 * \warning Shallow types
 *
 * These records are not self-contained, i.e. they depend on caches and will become invalid when the reader moves to
 * the next record.
 * Since some elements in the record are views, it may not be possible and/or safe to change all values.
 */
template <ownership own = ownership::shallow>
inline constexpr auto field_types_bcf_style =
  ttag<int32_t,                                                             // field::chrom,
       int32_t,                                                             // field::pos,
       std::string_view,                                                    // field::id,
       decltype(std::string_view{} | seqan3::views::char_to<seqan3::dna5>), // field::ref,
       std::vector<std::string_view>,                                       // field::alt,
       float,                                                               // field::qual,
       std::vector<int32_t>,                                                // field::filter,
       std::vector<info_element_bcf<ownership::shallow>>,                   // field::info,
       std::vector<genotype_element_bcf<ownership::shallow>>,               // field::genotypes,
       record_private_data>;                                                // field::_private

/*!\brief Deep field types for variant io.
 *!\ingroup var_io
 *
 * \details
 *
 * These field types result in a record that is self-contained, i.e. it does not depend on internal caches and the
 * state of the reader.
 *
 * Use these field types, if you intend to store individual records or if you need to change fields in the record
 * that are otherwise not modifiable (e.g. views).
 */
template <>
inline constexpr auto field_types_bcf_style<ownership::deep> =
  ttag<int32_t,                                            // field::chrom,
       int32_t,                                            // field::pos,
       std::string,                                        // field::id,
       std::vector<seqan3::dna5>,                          // field::ref,
       std::vector<std::string>,                           // field::alt,
       float,                                              // field::qual,
       std::vector<int32_t>,                               // field::filter,
       std::vector<info_element_bcf<ownership::deep>>,     // field::info,
       std::vector<genotype_element_bcf<ownership::deep>>, // field::genotypes,
       record_private_data>;                               // field::_private

/*!\brief Field types for variant IO that represent VCF more closely (text IDs etc).
 *!\ingroup var_io
 *
 * \details
 *
 * In contrast tp bio::var_io::field_types_bcf_style, these field types encode IDs as text and use
 * bio::var_io::genotypes_vcf_style to encode format/samples.
 *
 * If you know that you will be reading/writing almost exclusively VCF (and not BCF), using these field types
 * might lead to a better performance.
 *
 * \warning Shallow types
 *
 * These records are not self-contained, i.e. they depend on caches and will become invalid when the reader moves to
 * the next record.
 * Since some elements in the record are views, it may not be possible and/or safe to change all values.
 */
template <ownership own = ownership::shallow>
inline constexpr auto field_types_vcf_style =
  ttag<std::string_view,                                                    // field::chrom,
       int32_t,                                                             // field::pos,
       std::string_view,                                                    // field::id,
       decltype(std::string_view{} | seqan3::views::char_to<seqan3::dna5>), // field::ref,
       std::vector<std::string_view>,                                       // field::alt,
       float,                                                               // field::qual,
       std::vector<std::string_view>,                                       // field::filter,
       std::vector<info_element<ownership::shallow>>,                       // field::info,
       genotypes_vcf<ownership::shallow>,                                   // field::genotypes,
       record_private_data>;                                                // field::_private>;

/*!\brief Field types for variant IO that represent VCF more closely (text IDs etc); deep variant.
 *!\ingroup var_io
 *
 * \details
 *
 * The same as bio::var_io::field_types_vcf_style, but with self-contained records.
 */
template <>
inline constexpr auto field_types_vcf_style<ownership::deep> =
  ttag<std::string,                                // field::chrom
       int32_t,                                    // field::pos
       std::string,                                // field::id
       std::vector<seqan3::dna5>,                  // field::ref
       std::vector<std::string>,                   // field::alt
       float,                                      // field::qual
       std::vector<std::string>,                   // field::filter
       std::vector<info_element<ownership::deep>>, // field::info,
       genotypes_vcf<ownership::deep>,             // field::genotypes
       record_private_data>;                       // field::_private

//!\brief Every field is configured as a std::span of std::byte (this enables "raw" io).
//!\ingroup var_io
inline constexpr auto field_types_raw =
  seqan3::list_traits::concat<seqan3::list_traits::repeat<default_field_ids.size - 1, std::span<std::byte const>>,
                              seqan3::type_list<var_io::record_private_data>>{};

//!\}

/*!\brief Options that can be used to configure the behaviour of bio::var_io::reader.
 * \tparam field_ids_t   Type of the field_ids member (usually deduced).
 * \tparam field_types_t Type of the field_types member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup var_io
 *
 * \details
 *
 * TODO snippets
 *
 * ### Field type specific restrictions
 *
 * This section is only relevant if you specify the #field_types member manually, i.e. if you
 * change the field_types but do not use one of the predefined tags.
 *
 * 1. bio::field::chrom
 *   * string or string_view: The chromosome string is returned.
 *   * `int32_t`: The IDX value for the chromosome is returned.
 * 2. bio::field::pos
 *   * any integral: The position is returned as a number (`int32_t` recommended).
 * 3. bio::field::id
 *   * string or string_view: The ID as a string.
 * 4. bio::field::ref
 *   * string or string_view: plaintext.
 *   * `decltype(std::string_view{} | seqan3::views::char_strictly_to<seqan3::dna5>)`: A view
 * over a SeqAn3 alphabet. Other alphabets and/or transform views are also possible.
 * 5. bio::field::alt
 *   * back-insertable range of string or string_view: The ALTs as plaintext.
 *   * back-insertable range over views: similar views as for field::ref are supported but only
 * use this if you are sure there are no breakpoint strings etc. in the file!
 * 6. bio::field::qual
 *   * any arithmetic type: The quality as a number.
 * 7. bio::field::filter
 *   * back-insertable range of string or string_view: The filters as strings.
 *   * back-insertable range of `int32_t`: The IDX values of the filters.
 * 8. bio::field::info
 *   * back-insertable range of elements similar to bio::var_io::info_element
 *   * *similar* means any type decomposable into two elements (`struct` or tuple) where the
 * first is either a string[_view] or `int32_t` (IDX) and the second is bio::var_io::dynamic_type.
 * 9. field::genotypes
 *   1. A range (that supports back-insertion) over elements that are "similar" to
 * bio::var_io::genotype_element:
 *     * The elements must be decomposable into exactly two sub-elements (either `struct` or tuple).
 *     * The first subelement must be a string[_view] or `int32_t`.
 *     * The second subelement must bio::var_io::dynamic_vector_type.
 *   2. Or: A type similar to bio::var_io::genotypes_vcf :
 *     * It must be decomposable into exactly two sub-elements (either `struct` or tuple).
 *     * The first subelement must be a range over string[_views] that supports back-insertion.
 *     * The second subelement must range-of-range over bio::var_io::dynamic_type and both
 * range-dimensions need to support back-insertion.
 */
template <typename field_ids_t   = std::remove_cvref_t<decltype(default_field_ids)>,
          typename field_types_t = std::remove_cvref_t<decltype(field_types_bcf_style<ownership::shallow>)>,
          typename formats_t     = seqan3::type_list<vcf, bcf>>
struct reader_options
{
    //!\brief The fields that shall be contained in each record; a seqan3::tag over seqan3::field.
    field_ids_t field_ids = default_field_ids;

    /*!\brief The types corresponding to each field; a seqan3::type_tag over the types.
     *
     * \details
     *
     * See bio::var_io::reader_options for an overview of the supported field/type combinations.
     */
    field_types_t field_types = field_types_bcf_style<ownership::shallow>;

    /*!\brief The formats that input files can take; a seqan3::type_tag over the types.
     *
     * \details
     *
     * See bio::var_io::reader for an overview of the the supported formats.
     */
    formats_t formats = ttag<vcf, bcf>;

    //!\brief Whether to print non-critical file format warnings.
    bool print_warnings = true;

    //!\brief Options that are passed on to the internal stream oject.
    transparent_istream_options stream_options{};

    // TODO static_assert
};

} // namespace bio::var_io
