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
#include <seqan3/core/debug_stream/debug_stream_type.hpp>

#include <bio/detail/magic_get.hpp>
#include <bio/detail/range.hpp>
#include <bio/misc.hpp>
#include <bio/record.hpp>
#include <bio/var_io/dynamic_type.hpp>

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
inline t missing_value = t{};

//!\brief Specialisation for char.
template <>
inline constexpr char missing_value<char> = '.';

//!\brief Specialisation for std::string.
template <>
inline std::string missing_value<std::string> = ".";

//!\brief Specialisation for std::string_view.
template <>
inline constexpr std::string_view missing_value<std::string_view> = ".";

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
inline t end_of_vector = t{};

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

} // namespace bio::var_io

namespace bio::detail
{

//!\brief The "core" of a BCF record in bit-compatible representation to the on-disk format.
//!\ingroup var_io
struct bcf_record_core
{
    int32_t  chrom         = -1; //!< CHROM as IDX.
    int32_t  pos           = -1; //!< POS.
    int32_t  rlen          = -1; //!< Not used by this implementation.
    float    qual          = -1; //!< QUAL.
    uint16_t n_info        = 0;  //!< Number of INFOS values.
    uint16_t n_allele      = 0;  //!< Number of alleles.
    uint32_t n_sample : 24 = 0;  //!< Number of samples.
    uint8_t  n_fmt         = 0;  //!< Number of FORMAT values.
};

static_assert(sizeof(bcf_record_core) == 24, "Bit alignment problem in declaration of bcf_record_core.");

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
template <>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<int8_t> = bcf_type_descriptor::int8;
//!\brief Specialisation for int16.
template <>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<int16_t> = bcf_type_descriptor::int16;
//!\brief Specialisation for int32.
template <>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<int32_t> = bcf_type_descriptor::int32;
//!\brief Specialisation for float.
template <>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<float> = bcf_type_descriptor::float32;
//!\brief Specialisation for char.
template <>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<char> = bcf_type_descriptor::char8;
//!\}
//!\}
} // namespace bio::detail
