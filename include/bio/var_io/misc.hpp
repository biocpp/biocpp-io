// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::var_io::tag_dictionary class and auxiliaries.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_strictly_to.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>

#include <bio/detail/magic_get.hpp>
#include <bio/detail/range.hpp>
#include <bio/misc.hpp>
#include <bio/record.hpp>
#include <bio/var_io/dynamic_type.hpp>
#include <bio/var_io/header.hpp>

namespace bio::var_io
{

/*!\addtogroup var_io
 * \{
 */
/*!\name Special values in variant I/O.
 * \{
 */
//!\brief Default implementation.
template <typename t>
//     requires false
inline t missing_value = t{};

//!\brief Specialisation for char.
template <>
inline constexpr char missing_value<char> = char{0x07};

//!\brief Specialisation for integral types.
template <typename int_t>
    requires(std::same_as<int_t, int8_t> || std::same_as<int_t, int16_t> || std::same_as<int_t, int32_t>)
inline constexpr int_t missing_value<int_t> = std::numeric_limits<int_t>::lowest();

//!\brief Specialisation for float.
//!\hideinitializer
template <>
inline float missing_value<float> = []()
{
    static_assert(sizeof(uint32_t) == 4, "uint32_t is not four bytes big.");
    static_assert(sizeof(float) == 4, "float is not four bytes big.");
    union
    {
        uint32_t i;
        float    f;
    } u;
    u.i = 0x7F800001U;
    return u.f;
}();
//!\}
//!\}
} // namespace bio::var_io

namespace bio::detail
{
/*!\addtogroup var_io
 * \{
 */
/*!\name Special values in variant I/O.
 * \{
 */
//!\brief Default implementation. [not used]
template <typename t>
//     requires false
inline t end_of_vector = t{};

//!\brief Specialisation for char.
template <>
inline constexpr char end_of_vector<char> = '\0';

//!\brief Specialisation for integral types.
template <typename int_t>
    requires(std::same_as<int_t, int8_t> || std::same_as<int_t, int16_t> || std::same_as<int_t, int32_t>)
inline constexpr int_t end_of_vector<int_t> = std::numeric_limits<int_t>::lowest() + 1;

//!\brief Specialisation for float.
//!\hideinitializer
template <>
inline float end_of_vector<float> = []()
{
    static_assert(sizeof(uint32_t) == 4, "uint32_t is not four bytes big.");
    static_assert(sizeof(float) == 4, "float is not four bytes big.");
    union
    {
        uint32_t i;
        float    f;
    } u;
    u.i = 0x7F800002U;
    return u.f;
}();
//!\}
//!\}
} // namespace bio::detail

namespace bio::var_io
{

//-----------------------------------------------------------------------------
// Helper classes for field types
//-----------------------------------------------------------------------------

/*!\brief The type of elements in an INFO field. [default]
 * \ingroup var_io
 * \tparam own Ownership of the type; see bio::ownership.
 */
template <ownership own = ownership::shallow>
struct info_element
{
    //!\brief Type of the ID.
    using string_t = std::conditional_t<own == ownership::shallow, std::string_view, std::string>;

    //!\brief The ID of the element (as a string or string_view).
    string_t          id;
    //!\brief The value of the element.
    dynamic_type<own> value;

    //!\brief Defaulted three-way comparisons.
    auto operator<=>(info_element const &) const = default;
};

/*!\brief The type of elements in an INFO field. [full BCF-style]
 * \ingroup var_io
 * \tparam own Ownership of the type; see bio::ownership.
 */
template <ownership own = ownership::shallow>
struct info_element_bcf
{
    //!\brief The IDX of the element (index of that descriptor in the header).
    int32_t           idx;
    //!\brief The value of the element.
    dynamic_type<own> value;

    //!\brief Defaulted three-way comparisons.
    auto operator<=>(info_element_bcf const &) const = default;
};

/*!\brief A type representing an element in the GENOTYPES field.
 * \ingroup var_io
 *
 * \details
 *
 * Genotypes are represented as decribed in the BCF specification by default, i.e. information is grouped by
 * FORMAT identifier, not by sample.
 *
 * This element consists of the FORMAT ID given as a string and a vector of values inside a variant.
 * The size of the vector is:
 *
 *   * equal to the number of samples; or
 *   * 0 -- if the field is missing from all samples.
 *
 * The variant vector is guaranteed to be over the type defined in the header. Note that this is a vector over such
 * types (one element per sample!), so bio::var_io::dynamic_type_id::vector_of_int32 corresponds to
 * std::vector<std::vector<int32_t>>. See bio::var_io::dynamic_vector_type for more details.
 *
 * If fields are missing from some samples but not others, the vector will have full size but the respective values
 * will be set to the missing value (see bio::var_io::missing_value) or be the empty vector (in case the element type
 * is a vector).
 */
template <ownership own = ownership::shallow>
struct genotype_element
{
    //!\brief Type of the ID.
    using string_t = std::conditional_t<own == ownership::shallow, std::string_view, std::string>;

    //!\brief The ID of the element (as a string or string_view).
    string_t                 id;
    //!\brief The value of the element.
    dynamic_vector_type<own> value;

    //!\brief Defaulted three-way comparisons.
    auto operator<=>(genotype_element const &) const = default;
};

/*!\brief A type representing an element in the GENOTYPES field. [full BCF-style]
 * \ingroup var_io
 *
 * \details
 *
 * The same as bio::var_io::genotype_element except that a numeric IDX is used instead of the ID string.
 */
template <ownership own = ownership::shallow>
struct genotype_element_bcf
{
    //!\brief The IDX of the element (index of that descriptor in the header).
    int32_t                  idx;
    //!\brief The value of the element.
    dynamic_vector_type<own> value;

    //!\brief Defaulted three-way comparisons.
    auto operator<=>(genotype_element_bcf const &) const = default;
};

/*!\brief A type representing the FORMATS column and all sample columns in VCF-style.
 * \ingroup var_io
 *
 * \details
 *
 * This type can be used as the field-type for the GENOTYPES field as an alternative to a range of
 * bio::var_io::genotype_element.
 *
 * It uses the data layout as it appears in a VCF file with the FORMAT strings in one member and a vector of "samples".
 * Each element of that vector represents a single sample column and is implemented as a vector of values of
 * dynamic type (see bio::var_io::dynamic_type).
 *
 * **This data layout is not recommended, because it is almost always slower.**
 * Use it only, if you know that the user will never read or write BCF and if you do very little processing of the
 * sample values.
 */
template <ownership own = ownership::shallow>
struct genotypes_vcf
{
    //!\brief Type of the format strings.
    using string_t = std::conditional_t<own == ownership::shallow, std::string_view, std::string>;

    //!\brief The FORMAT strings.
    std::vector<string_t>                       format_strings;
    //!\brief The sample columns.
    std::vector<std::vector<dynamic_type<own>>> samples;

    //!\brief Defaulted three-way comparisons.
    auto operator<=>(genotypes_vcf const &) const = default;
};

//!\brief A datastructure that contains private data of variant IO records.
//!\ingroup var_io
struct record_private_data
{
    //!\privatesection
    //!\brief Pointer to the header
    header const * header_ptr = nullptr;

    //!\brief Defaulted three-way comparison.
    friend bool operator==(record_private_data const &, record_private_data const &) = default;
    // TODO pointer to bcf-record
};

//-----------------------------------------------------------------------------
// default_field_ids
//-----------------------------------------------------------------------------

//!\brief Default fields for bio::var_io::reader_options.
//!\ingroup var_io
inline constinit auto default_field_ids = vtag<field::chrom,
                                               field::pos,
                                               field::id,
                                               field::ref,
                                               field::alt,
                                               field::qual,
                                               field::filter,
                                               field::info,
                                               field::genotypes,
                                               field::_private>;

//-----------------------------------------------------------------------------
// Pre-defined field types (reader)
//-----------------------------------------------------------------------------

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
 * The "style" of the record resembles the VCF specification, i.e. contigs, FILTERs and INFO identifiers are
 * represented as string/string_views. **However,**  the genotypes are encoded by-genotype (BCF-style) and not by-sample
 *(VCF-style) for performance reasons.
 *
 * See bio::var_io::genotypes_bcf_style for more information on the latter.
 */
template <ownership own = ownership::shallow>
inline constinit auto field_types =
  ttag<std::string_view,                                                             // field::chrom,
       int32_t,                                                                      // field::pos,
       std::string_view,                                                             // field::id,
       decltype(std::string_view{} | seqan3::views::char_strictly_to<seqan3::dna5>), // field::ref,
       std::vector<std::string_view>,                                                // field::alt,
       float,                                                                        // field::qual,
       std::vector<std::string_view>,                                                // field::filter,
       std::vector<info_element<ownership::shallow>>,                                // field::info,
       std::vector<genotype_element<ownership::shallow>>,                            // field::genotypes,
       record_private_data>;                                                         // field::_private

//!\brief Deep version of bio::var_io::field_types.
//!\ingroup var_io
template <>
inline constinit auto field_types<ownership::deep> =
  ttag<std::string,                                    // field::chrom,
       int32_t,                                        // field::pos,
       std::string,                                    // field::id,
       std::vector<seqan3::dna5>,                      // field::ref,
       std::vector<std::string>,                       // field::alt,
       float,                                          // field::qual,
       std::vector<std::string>,                       // field::filter,
       std::vector<info_element<ownership::deep>>,     // field::info,
       std::vector<genotype_element<ownership::deep>>, // field::genotypes,
       record_private_data>;                           // field::_private

/*!\brief Alternative set of field types (BCF-style, shallow).
 *!\ingroup var_io
 *
 * \details
 *
 * See bio::var_io::reader_options for when and why to choose these field types.
 */
template <ownership own = ownership::shallow>
inline constinit auto field_types_bcf_style =
  ttag<int32_t,                                                                      // field::chrom,
       int32_t,                                                                      // field::pos,
       std::string_view,                                                             // field::id,
       decltype(std::string_view{} | seqan3::views::char_strictly_to<seqan3::dna5>), // field::ref,
       std::vector<std::string_view>,                                                // field::alt,
       float,                                                                        // field::qual,
       std::vector<int32_t>,                                                         // field::filter,
       std::vector<info_element_bcf<ownership::shallow>>,                            // field::info,
       std::vector<genotype_element_bcf<ownership::shallow>>,                        // field::genotypes,
       record_private_data>;                                                         // field::_private

/*!\brief Alternative set of field types (BCF-style, deep).
 *!\ingroup var_io
 *
 * \details
 *
 * See bio::var_io::reader_options for when and why to choose these field types.
 */
template <>
inline constinit auto field_types_bcf_style<ownership::deep> =
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

/*!\brief Alternative set of field types (VCF-style, shallow).
 *!\ingroup var_io
 *
 * \details
 *
 * See bio::var_io::reader_options for when and why to choose these field types.
 */
template <ownership own = ownership::shallow>
inline constinit auto field_types_vcf_style =
  ttag<std::string_view,                                                             // field::chrom,
       int32_t,                                                                      // field::pos,
       std::string_view,                                                             // field::id,
       decltype(std::string_view{} | seqan3::views::char_strictly_to<seqan3::dna5>), // field::ref,
       std::vector<std::string_view>,                                                // field::alt,
       float,                                                                        // field::qual,
       std::vector<std::string_view>,                                                // field::filter,
       std::vector<info_element<ownership::shallow>>,                                // field::info,
       genotypes_vcf<ownership::shallow>,                                            // field::genotypes,
       record_private_data>;                                                         // field::_private>;

/*!\brief Alternative set of field types (BCF-style, deep).
 *!\ingroup var_io
 *
 * \details
 *
 * See bio::var_io::reader_options for when and why to choose these field types.
 */
template <>
inline constinit auto field_types_vcf_style<ownership::deep> =
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
inline constinit auto field_types_raw =
  seqan3::list_traits::concat<seqan3::list_traits::repeat<default_field_ids.size - 1, std::span<std::byte const>>,
                              seqan3::type_list<var_io::record_private_data>>{};

//!\}

/*!\brief A alias for bio::record that is usable with variant IO.
 * \ingroup var_io
 * \details
 *
 * This alias is provided purely for convenience. See the documentation for
 * bio::var_io::writer for an example of how to use it.
 */
template <ownership own = ownership::deep>
using default_record = record<decltype(default_field_ids), decltype(field_types<own>)>;

} // namespace bio::var_io

namespace bio::detail
{

//-----------------------------------------------------------------------------
// BCF record core
//-----------------------------------------------------------------------------

//!\brief The "core" of a BCF record in bit-compatible representation to the on-disk format.
//!\ingroup var_io
struct bcf_record_core
{
    int32_t  chrom         = -1; //!< CHROM as IDX.
    int32_t  pos           = -1; //!< POS.
    int32_t  rlen          = -1; //!< Not used by this implementation.
    float    qual          = bio::var_io::missing_value<float>; //!< QUAL.
    uint16_t n_info        = 0;  //!< Number of INFOS values.
    uint16_t n_allele      = 0;  //!< Number of alleles.
    uint32_t n_sample : 24 = 0;  //!< Number of samples.
    uint8_t  n_fmt         = 0;  //!< Number of FORMAT values.
};

static_assert(sizeof(bcf_record_core) == 24, "Bit alignment problem in declaration of bcf_record_core.");

//-----------------------------------------------------------------------------
// bcf_type_descriptor and utilities
//-----------------------------------------------------------------------------

/*!\name Pre-defined field types
 * \brief These can be used to configure the behaviour of the bio::var_io::reader via bio::var_io::reader_options.
 * \{
 */

//!\brief The BCF type descriptor with values as described in the specification.
//!\ingroup var_io
enum class bcf_type_descriptor : uint8_t
{
    missing = 0,
    int8    = 1,
    int16   = 2,
    int32   = 3,
    float32 = 5,
    char8   = 7
};

//TODO this needs standalone tests
//!\brief Compute the smallest possible integral type descriptor able to represent the value.
//!\relates bcf_type_descriptor
detail::bcf_type_descriptor smallest_int_desc(std::unsigned_integral auto const num)
{
    // bits required to represent number (the +1 because signed integral has smaller range)
    switch (std::bit_ceil(std::bit_width(num) + 1))
    {
        case 128:
        case 64:
            throw std::runtime_error{std::string{"Could not write number '"} + detail::to_string(num) + "'. Value out of range (only int32 supported)."};
            return {};
        case 32:
            return detail::bcf_type_descriptor::int32;
        case 16:
            return detail::bcf_type_descriptor::int16;
        default:
            return detail::bcf_type_descriptor::int8;
    }
}

//!\overload
detail::bcf_type_descriptor smallest_int_desc(std::signed_integral auto const num)
{
    return smallest_int_desc(static_cast<uint64_t>(std::abs(num)));
}

//!\overload
detail::bcf_type_descriptor smallest_int_desc(std::ranges::forward_range auto && range)
{
//         //TODO check if this is faster:
//         val_t max = std::numeric_limits<val_t>::lowest();
//         for (val_t elem : data)
//         {
//             if (elem > max)
//             {
//                 max = elem;
//                 desc = smallest_int_desc(elem);
//                 if (desc == detail::bcf_type_descriptor::int32) // this is max(type_descriptor)
//                     break;
//             }
//         }
    return smallest_int_desc(std::ranges::empty(range) ? 0 : *std::ranges::max_element(range));
}

//!\overload
template <std::ranges::forward_range rng_t>
    requires std::ranges::range<std::ranges::range_reference_t<rng_t>>
detail::bcf_type_descriptor smallest_int_desc(rng_t & range)
{
    return smallest_int_desc(range | std::views::join);
}

/*!\addtogroup var_io
 * \{
 */
/*!\name A compile-time mapping of types to bio::detail::bcf_type_descriptor
 * \{
 */
//!\brief Default implementation.
template <typename t>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor = bcf_type_descriptor::missing;

//!\brief Specialisation for int8.
template <std::signed_integral t>
    requires (sizeof(t) == 1)
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::int8;

//!\brief Specialisation for int16.
template <std::signed_integral t>
    requires (sizeof(t) == 2)
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::int16;
//!\brief Specialisation for int16.
template <std::unsigned_integral t>
    requires (sizeof(t) == 1)
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::int16;

//!\brief Specialisation for int32.
template <std::signed_integral t>
    requires (sizeof(t) >= 4)
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::int32;

//!\brief Specialisation for int32.
template <std::unsigned_integral t>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::int32;

//!\brief Specialisation for float.
template <>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<float> = bcf_type_descriptor::float32;
//!\brief Specialisation for double.
template <>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<double> = bcf_type_descriptor::float32;

//!\brief Specialisation for char.
template <>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<char> = bcf_type_descriptor::char8;
//!\brief Specialisation for seqan3 alphabets.
template <detail::deliberate_alphabet t>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::char8;
//!\}

//!\}
} // namespace bio::detail
