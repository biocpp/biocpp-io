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
//!\brief Default implementation. [not used]
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

/*!\brief A type representing an variant file INFO field [index of the INFO in header, value].
 * \ingroup var_io
 */
template <ownership own = ownership::shallow>
using info_element = std::pair<int32_t, dynamic_type<own>>;

/*!\brief A type representing an genotype.
 * \ingroup var_io
 *
 * \details
 *
 * Genotypes / samples are represented as decribed in the BCF specification, i.e. information is grouped by FORMAT
 * identifier, not by sample.
 *
 * This entry consists of the FORMAT index in the file's header and a vector of values. The size of the vector is:
 *
 *   * equal to the number of samples; or
 *   * 0 -- if the field is missing from all samples.
 *
 * The variant vector is guaranteed to be over the type defined in the header. Note that this is a vector over such
 * types (one element per sample!), so bio::var_io::dynamic_type_id::vector_of_int32 corresponds to
 * std::vector<std::vector<int32_t>>. See seqan3::dynamic_vector_type for more details.
 *
 * If fields are missing from some samples but not others, the vector will have full size but the respective values
 * will be set to the missing value (see bio::var_io::is_missing()) or be the empty vector (in case the element type
 * is a vector).
 */
template <ownership own = ownership::shallow>
using genotype_bcf_style = std::pair<int32_t, dynamic_vector_type<own>>;

//!\brief TODO this will get a rework
//!\ingroup var_io
template <ownership own = ownership::shallow>
using genotypes_bcf_style = std::vector<genotype_bcf_style<own>>;

//!\brief TODO this will get a rework
//!\ingroup var_io
template <ownership own = ownership::shallow>
using genotypes_vcf_style =
  std::pair<std::vector<std::conditional_t<own == ownership::shallow, std::string_view, std::string>>,
            std::vector<std::vector<dynamic_type<own>>>>;

//!\brief TODO this will get a rework
//!\ingroup var_io
template <ownership own = ownership::shallow>
using genotypes_as_strings = std::vector<std::conditional_t<own == ownership::shallow, std::string_view, std::string>>;

//!\brief Default fields for bio::var_io::reader_options.
//!\ingroup var_io
inline constexpr auto default_field_ids = vtag<field::chrom,
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
