// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::ann_io::reader_options.
 * \author Joshua Kim <kim_j AT molgen.mpg.de>
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include <seqan3/utility/type_list/traits.hpp>

#include <bio/detail/misc.hpp>
#include <bio/format/bed.hpp>
#include <bio/stream/transparent_istream.hpp>
#include <bio/stream/transparent_ostream.hpp>

#include <bio/ann_io/misc.hpp>

namespace bio::ann_io
{

//-----------------------------------------------------------------------------
// reader_options
//-----------------------------------------------------------------------------

/*!\brief Options that can be used to configure the behaviour of bio::ann_io::reader.
 * \tparam field_ids_t   Type of the field_ids member (usually deduced).
 * \tparam field_types_t Type of the field_types member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup ann_io
 *
 * \details
 *
 * This object can be configured in a similar way as bio::seq_io::reader_options.
 * If you are new to the way options are set in this library, have a look bio::seq_io::reader
 * and bio::seq_io::reader_options first, as those are much simpler.
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
 * with one of the predefined sets.
 *
 * ### Pre-defined tags
 *
 * Two "styles" of field types are predefined:
 *
 * 1. bio::ann_io::field_types (the default)
 *   * All "strings" are represented as strings.
 * 2. bio::ann_io::field_types_bcf_style (BCF-style)
 *   * Most "strings" are represented by their in-header IDX value (see the BCF spec for more details).
 *   * When reading and writing, you need to make sure that the IDX values in the output header are the same as in the
 * input header, otherwise your record fields might change meaning or even become invalid.
 *
 * Both styles are "shallow" by default, but can be configured to be "deep":
 *
 * 1. shallow (bio::ownership::shallow)
 *   * The record contains light-weight data structures like views.
 *   * Record cannot be "stored"; it depends on internal caches and buffers, and it becomes invalid
 * as soon as the next record is read from the file.
 * 2. deep (bio::ownership::deep)
 *   * The record is self-contained; sequences and strings are stored in containers.
 *   * Record can be copied or stored and can "live on" independently of the reader.
 *
 * This example shows how to use deep records:
 *
 * \snippet test/snippet/ann_io/ann_io_reader_options.cpp field_types_deep
 *
 * Performance implications:
 *   * Shallow records imply fewer allocations and lower overhead during reading.
 *   * If you know that you need to copy your fields anyway, using a deep record can be faster than using a shallow
 * record and copying the data "manually" out of that (because certain internal caches are re-used to create deep
 * records).
 *   * field_types_bcf_style<ownership::deep> is faster than field_types<ownership::deep>, but for the shallow variants
 * there is almost no difference.
 *
 * TODO some of this should be moved to a general documentation page on configuring records; shallow vs deep; etc
 *
 * ### Manual configuration
 *
 * This section is only relevant if you specify the #field_types member manually via
 * a bio::ttag, i.e. if you change the field_types but do not use one of the predefined tags
 * (see above).
 *
 * The following types are valid for the respective fields and you can mix-and-match shallow/deep and integral/text IDs:
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
 *   * back-insertable range over seqan3::alphabet (a container with converted elements).
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
 *   * back-insertable range of elements "similar" to bio::ann_io::info_element:
 *     * The elements must be decomposable into two subelements (`struct` or tuple).
 *     * The first subelement must be either a string[_view] (ID) or `int32_t` (IDX).
 *     * The second subelement must be bio::ann_io::info_element_value_type.
 * 9. field::genotypes
 *   * back-insertable range of elements "similar" to bio::ann_io::genotype_element:
 *     * The elements must be decomposable into exactly two sub-elements (either `struct` or tuple).
 *     * The first subelement must be either a string[_view] (ID) or `int32_t` (IDX).
 *     * The second subelement must bio::ann_io::genotype_element_value_type.
 *
 * This example shows how to read only a subset of the available fields and manually specify their type:
 *
 * \snippet test/snippet/ann_io/ann_io_reader_options.cpp field_types_expert
 *
 * Reading fewer fields than available may provide a noticeable speed-up since only the
 * requested fields are actually parsed. Any field may also be set to `std::span<std::byte const>` which
 * results in no parsing happening for that field.
 *
 */
template <typename field_ids_t   = decltype(default_field_ids),
          typename field_types_t = decltype(field_types<ownership::shallow>),
          typename formats_t     = seqan3::type_list<bed>>
struct reader_options
{
    //!\brief The fields that shall be contained in each record; a seqan3::tag over seqan3::field.
    field_ids_t field_ids = default_field_ids;

    /*!\brief The types corresponding to each field; a bio::ttag over the types.
     *
     * \details
     *
     * See bio::ann_io::reader_options for an overview of the supported field/type combinations.
     */
    field_types_t field_types = bio::ann_io::field_types<ownership::shallow>;

    /*!\brief The formats that input files can take; a bio::ttag over the types.
     *
     * \details
     *
     * See bio::ann_io::reader for an overview of the the supported formats.
     */
    formats_t formats = ttag<bed>;

    //!\brief Whether to print non-critical file format warnings.
    bool print_warnings = true;

    //!\brief Options that are passed on to the internal stream oject.
    transparent_istream_options stream_options{};

private:
    static_assert(detail::is_fields_tag<field_ids_t>, "field_ids must be a bio::vtag over bio::field.");

    static_assert(detail::is_type_list<field_types_t>, "field_types must be a bio::ttag / seqan3::type_list.");

    static_assert(detail::is_type_list<formats_t>, "formats must be a bio::ttag / seqan3::type_list.");

    static_assert(field_ids_t::size == field_types_t::size(), "field_ids and field_types must have the same size.");

    //!\brief Type of the record.
    using record_t = record<field_ids_t, field_types_t>;

    static_assert(
      detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
        !field_ids_t::contains(field::chrom) ||
        detail::back_insertable_with<record_element_t<field::chrom, rec_t>, char> ||
        detail::one_of<record_element_t<field::chrom, rec_t>, std::string_view>) { return std::true_type{}; }),
      "Requirements for the field-type of the CHROM-field not met. See documentation for bio::ann_io::reader_options.");

    static_assert(
      detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
        !field_ids_t::contains(field::chromStart) ||
        std::integral<std::remove_reference_t<record_element_t<field::chromStart, rec_t>>>) { return std::true_type{}; }),
      "Requirements for the field-type of the chromStart-field not met. See documentation for bio::ann_io::reader_options.");

    static_assert(
      detail::lazy_concept_checker([]<typename rec_t = record_t>(auto) requires(
        !field_ids_t::contains(field::chromEnd) ||
        std::integral<std::remove_reference_t<record_element_t<field::chromEnd, rec_t>>>) { return std::true_type{}; }),
      "Requirements for the field-type of the chromEnd-field not met. See documentation for bio::ann_io::reader_options.");
};

} // namespace bio::ann_io
