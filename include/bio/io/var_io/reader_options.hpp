// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::var_io::reader_options.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include <bio/meta/tag/ttag.hpp>
#include <bio/meta/type_list/traits.hpp>

#include <bio/io/detail/misc.hpp>
#include <bio/io/format/bcf.hpp>
#include <bio/io/format/vcf.hpp>
#include <bio/io/genomic_region.hpp>
#include <bio/io/stream/transparent_istream.hpp>
#include <bio/io/stream/transparent_ostream.hpp>
#include <bio/io/var_io/header.hpp>
#include <bio/io/var_io/misc.hpp>

//-----------------------------------------------------------------------------
// field concepts
//-----------------------------------------------------------------------------

namespace bio::io::detail
{

/*!\interface bio::io::detail::info_element_reader_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var_io::info_element / bio::io::var_io::info_element_bcf.
 */
//!\cond CONCEPT_DEF
template <typename t>
concept info_element_reader_concept = detail::decomposable_into_two<t> &&
  (detail::out_string<detail::first_elem_t<t>> ||
   std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::is_info_element_value_type<detail::second_elem_t<t>>;
//!\endcond

/*!\interface bio::io::detail::genotype_reader_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var_io::genotype_element / bio::io::var_io::genotype_element_bcf.
 */
//!\cond CONCEPT_DEF
template <typename t>
concept genotype_reader_concept = detail::decomposable_into_two<t> &&
  (detail::out_string<detail::first_elem_t<t>> ||
   std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::is_genotype_element_value_type<detail::second_elem_t<t>>;
//!\endcond
} // namespace bio::io::detail

namespace bio::io::var_io
{

//-----------------------------------------------------------------------------
// reader_options
//-----------------------------------------------------------------------------

/*!\brief Options that can be used to configure the behaviour of bio::io::var_io::reader.
 * \tparam field_ids_t   Type of the field_ids member (usually deduced).
 * \tparam field_types_t Type of the field_types member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup var_io
 *
 * \details
 *
 * This object can be configured in a similar way as bio::io::seq_io::reader_options.
 * If you are new to the way options are set in this library, have a look bio::io::seq_io::reader
 * and bio::io::seq_io::reader_options first, as those are much simpler.
 *
 * ## Field types
 *
 * The internal representation of VCF and BCF are different. To be able to freely
 * interchange between these formats, this library needs to choose one representation that
 * everything is converted to when being read.
 *
 * Changing the field_types member configures the reader to return data in different types.
 * One thing that is fixed for all configurations in this library is the layout of the GENOTYPES field
 * which is always grouped "by-field" (BCF-style) and not "by-sample" (VCF-style).
 * Another important choice is that **numbers are always 1-based,** because this is the default in VCF
 * and all other tools that deal with VCF/BCF.
 *
 * Beyond that, a wide variety of types are supported per field (see below), but most users will be happy
 * with one of the predefined tags.
 *
 * ### Pre-defined tags
 *
 * Two "styles" of field types are predefined:
 *
 * 1. bio::io::var_io::field_types (the default)
 *   * All "strings" are represented as strings.
 * 2. bio::io::var_io::field_types_bcf_style (BCF-style)
 *   * Most "strings" are represented by their in-header IDX value (see the BCF spec for more details).
 *   * When reading and writing, you need to make sure that the IDX values in the output header are the same as in the
 * input header, otherwise your record fields might change meaning or even become invalid.
 *
 * Both styles are *shallow* by default, but can be configured to be *deep*.
 * For more details, see \ref shallow_vs_deep
 *
 * ### Manual configuration
 *
 * This section is only relevant if you specify the #field_types member manually via
 * a bio::meta::ttag, i.e. if you change the field_types but do not use one of the predefined tags
 * (see above).
 *
 * The following types are valid for the respective fields and you can mix-and-match shallow/deep and integral/text IDs:
 *
 * 1. bio::io::field::chrom
 *   * string or string_view: The chromosome string is returned.
 *   * `int32_t`: The IDX value for the chromosome is returned.
 * 2. bio::io::field::pos
 *   * any integral: The position is returned as a number (`int32_t` recommended).
 * 3. bio::io::field::id
 *   * string or string_view: The ID as a string.
 * 4. bio::io::field::ref
 *   * string or string_view: plaintext.
 *   * back-insertable range over bio::alphabet::alphabet (a container with converted elements).
 *   * `decltype(std::string_view{} | bio::views::char_strictly_to<alphabet::dna5>)`: A view
 * over a SeqAn3 alphabet. Other alphabets and/or transform views are also possible.
 * 5. bio::io::field::alt
 *   * back-insertable range of string or string_view: The ALTs as plaintext.
 *   * back-insertable range over views: similar views as for field::ref are supported but only
 * use this if you are sure there are no breakpoint strings etc. in the file!
 * 6. bio::io::field::qual
 *   * any arithmetic type: The quality as a number.
 * 7. bio::io::field::filter
 *   * back-insertable range of string or string_view: The filters as strings.
 *   * back-insertable range of `int32_t`: The IDX values of the filters.
 * 8. bio::io::field::info
 *   * back-insertable range of elements "similar" to bio::io::var_io::info_element:
 *     * The elements must be decomposable into two subelements (`struct` or tuple).
 *     * The first subelement must be either a string[_view] (ID) or `int32_t` (IDX).
 *     * The second subelement must be bio::io::var_io::info_element_value_type.
 * 9. field::genotypes
 *   * back-insertable range of elements "similar" to bio::io::var_io::genotype_element:
 *     * The elements must be decomposable into exactly two sub-elements (either `struct` or tuple).
 *     * The first subelement must be either a string[_view] (ID) or `int32_t` (IDX).
 *     * The second subelement must bio::io::var_io::genotype_element_value_type.
 *
 * This example shows how to read only a subset of the available fields and manually specify their type:
 *
 * \snippet test/snippet/var_io/var_io_reader_options.cpp field_types_expert
 *
 * Reading fewer fields than available may provide a noticeable speed-up since only the
 * requested fields are actually parsed.
 */
template <typename field_ids_t   = decltype(default_field_ids),
          typename field_types_t = decltype(field_types<ownership::shallow>),
          typename formats_t     = meta::type_list<vcf, bcf>>
struct reader_options
{
    /*!\name Field configuration
     * \brief These options allow configuring the record's field order and types.
     * \{
     */
    //!\brief The fields that shall be contained in each record; a bio::meta::vtag over bio::io::field.
    field_ids_t field_ids = default_field_ids;

    /*!\brief The types corresponding to each field; a bio::meta::ttag over the types.
     *
     * \details
     *
     * See bio::io::var_io::reader_options for an overview of the supported field/type combinations.
     */
    field_types_t field_types = bio::io::var_io::field_types<ownership::shallow>;
    //!\}

    /*!\brief The formats that input files can take; a bio::meta::ttag over the types.
     *
     * \details
     *
     * See bio::io::var_io::reader for an overview of the the supported formats.
     */
    formats_t formats = meta::ttag<vcf, bcf>;

    //!\brief Whether to print non-critical file format warnings.
    bool print_warnings = true;

    /*!\name Region filtering
     * \brief These options allow filtering the file for a sub-region.
     * \{
     */
    /*!\brief Only display records that overlap with the given region (ignored if default-initialised).
     * \details
     *
     * The region must be specified as a **0-based, half-open interval** (even though VCF is 1-based).
     *
     * This option works with and without an index, but reading without an index is disabled by default.
     * Indexes are not usable when the format is uncompressed VCF or if the input is read from STDIN.
     *
     * ### Region filter with an index
     *
     * An index file allows skipping regions of the file on-disk. This is usually a lot faster than filtering
     * without an index. The filename is detected automatically but can also be specified manually (see
     * #region_index_file).
     *
     * ### Region filter without an index
     *
     * The file will be scanned linearly for records that overlap the region.
     * Only minimal parts of the record are parsed to compare the positions, so this is faster than
     * using a `std::views::filter` on the file.
     * This option is disabled by default, so an exception is thrown when a region is given but no
     * index file is found. To allow region filtering without an index, set #region_index_optional to
     * `true`.
     */
    genomic_region<ownership::deep> region{};

    /*!\brief Path to the index file [optional, auto-detected if not specified].
     * \details
     *
     * This option is ignored if no #region is specified, if the file format is uncompressed VCF or the data
     * are read from standard input.
     *
     * If no path is given, the library will look for `FILE.tbi` where `FILE` is the filename of the input
     * file. CSI-indexes are not (yet) supported.
     */
    std::filesystem::path region_index_file{};

    //!\brief Allow linear-time region filtering when no index file is available.
    bool region_index_optional = false;
    //!\}

    //!\brief Options that are passed on to the internal stream oject.
    transparent_istream_options stream_options{};

private:
    static_assert(detail::is_fields_tag<field_ids_t>, "field_ids must be a bio::meta::vtag over bio::io::field.");

    static_assert(meta::detail::is_type_list<field_types_t>,
                  "field_types must be a bio::meta::ttag / bio::meta::type_list.");

    static_assert(meta::detail::is_type_list<formats_t>, "formats must be a bio::meta::ttag / bio::meta::type_list.");

    static_assert(field_ids_t::size == field_types_t::size(), "field_ids and field_types must have the same size.");

    //!\brief Type of the record.
    using record_t = record<field_ids_t, field_types_t>;

    static_assert(detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
                    !field_ids_t::contains(field::chrom) ||
                    detail::back_insertable_with<record_element_t<field::chrom, rec_t>, char> ||
                    detail::one_of<record_element_t<field::chrom, rec_t>, std::string_view, int32_t>) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the CHROM-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    static_assert(detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
                    !field_ids_t::contains(field::pos) ||
                    std::integral<std::remove_reference_t<record_element_t<field::pos, rec_t>>>) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the POS-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    static_assert(detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
                    !field_ids_t::contains(field::id) ||
                    detail::back_insertable_with<record_element_t<field::id, rec_t>, char> ||
                    detail::one_of<std::remove_reference_t<record_element_t<field::id, rec_t>>, std::string_view>) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the ID-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");
    // TODO re-activate me later
    //      static_assert(detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
    //                      !field_ids_t::contains(field::ref) ||
    //                      (detail::back_insertable<record_element_t<field::ref, rec_t>> &&
    //                       alphabet::alphabet<std::ranges::range_reference_t<record_element_t<field::ref, rec_t>>>) ||
    //                      std::same_as<std::remove_reference_t<record_element_t<field::ref, rec_t>>, std::string_view>
    //                      || detail::transform_view_on_string_view<record_element_t<field::ref, rec_t>>) {
    //                        return std::true_type{};
    //                    }),
    //                    "Requirements for the field-type of the REF-field not met. See documentation for "
    //                    "bio::io::var_io::reader_options.");

    // TODO re-activate me later
    //      static_assert(
    //        detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
    //          !field_ids_t::contains(field::alt) ||
    //          (detail::back_insertable<record_element_t<field::alt, rec_t>> &&
    //             (detail::back_insertable<std::ranges::range_reference_t<record_element_t<field::alt, rec_t>>> &&
    //              alphabet::alphabet<
    //                std::ranges::range_reference_t<std::ranges::range_reference_t<record_element_t<field::alt,
    //                rec_t>>>>) ||
    //           std::same_as<std::remove_reference_t<std::ranges::range_reference_t<record_element_t<field::alt,
    //           rec_t>>>,
    //                        std::string_view> ||
    //           detail::transform_view_on_string_view<std::ranges::range_reference_t<record_element_t<field::alt,
    //           rec_t>>>)) {
    //            return std::true_type{};
    //        }),
    //        "Requirements for the field-type of the ALT-field not met. See documentation for "
    //        "bio::io::var_io::reader_options.");

    static_assert(detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
                    !field_ids_t::contains(field::qual) ||
                    seqan3::arithmetic<std::remove_reference_t<record_element_t<field::qual, rec_t>>>) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the QUAL-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    static_assert(
      detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
        !field_ids_t::contains(field::filter) ||
        (detail::back_insertable<record_element_t<field::filter, rec_t>> &&
         (detail::back_insertable_with<std::ranges::range_reference_t<record_element_t<field::filter, rec_t>>, char> ||
          detail::one_of<
            std::remove_reference_t<std::ranges::range_reference_t<record_element_t<field::filter, rec_t>>>,
            std::string_view,
            int32_t>))) { return std::true_type{}; }),
      "Requirements for the field-type of the FILTER-field not met. See documentation for "
      "bio::io::var_io::reader_options.");

    static_assert(
      detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
        !field_ids_t::contains(field::info) ||
        (detail::back_insertable<record_element_t<field::info, rec_t>> &&
         detail::info_element_reader_concept<
           std::remove_reference_t<std::ranges::range_reference_t<record_element_t<field::info, rec_t>>>>)) {
          return std::true_type{};
      }),
      "Requirements for the field-type of the INFO-field not met. See documentation for "
      "bio::io::var_io::reader_options.");

    static_assert(
      detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
        !field_ids_t::contains(field::genotypes) ||
        (detail::back_insertable<record_element_t<field::genotypes, rec_t>> &&
         detail::genotype_reader_concept<
           std::remove_reference_t<std::ranges::range_reference_t<record_element_t<field::genotypes, rec_t>>>>)) {
          return std::true_type{};
      }),
      "Requirements for the field-type of the GENOTYPES-field not met. See documentation for "
      "bio::io::var_io::reader_options.");
};

} // namespace bio::io::var_io
