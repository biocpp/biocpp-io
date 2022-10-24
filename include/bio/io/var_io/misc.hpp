// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::var_io::tag_dictionary class and auxiliaries.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/meta/tag/ttag.hpp>
#include <bio/meta/tag/vtag.hpp>
#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/detail/magic_get.hpp>
#include <bio/io/detail/range.hpp>
#include <bio/io/detail/tuple_record.hpp>
#include <bio/io/misc.hpp>

namespace bio::io::var_io::detail
{} // namespace bio::io::var_io::detail

//-----------------------------------------------------------------------------
// missing_value
//-----------------------------------------------------------------------------

namespace bio::io::var_io
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
inline constexpr char missing_value<char> = char{0x07};

//!\brief Specialisation for integral types.
template <meta::one_of<int8_t, int16_t, int32_t> int_t>
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
} // namespace bio::io::var_io

//-----------------------------------------------------------------------------
// end_of_vector
//-----------------------------------------------------------------------------

namespace bio::io::var_io::detail
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

//!\brief Specialisation for char.
template <>
inline constexpr char end_of_vector<char> = '\0';

//!\brief Specialisation for integral types.
template <meta::one_of<int8_t, int16_t, int32_t> int_t>
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
} // namespace bio::io::var_io::detail

namespace bio::io::var_io
{

//-----------------------------------------------------------------------------
// value_type_id
//-----------------------------------------------------------------------------

//!\brief Enumerator to ease "dynamic typing" in variant IO.
//!\ingroup var_io
enum class value_type_id : size_t
{
    char8,             //!< Used for "Character" fields of size 1.
    int8,              //!< Used for "Integer" fields of size 1 where the value fits in one byte.
    int16,             //!< Used for "Integer" fields of size 1 where the value fits in two bytes.
    int32,             //!< Used for "Integer" fields of size 1 where the value fits in four bytes.
    float32,           //!< Used for "Float" fields of size 1.
    string,            //!< Used for "String" fields of size 1 and "Character" fields of size != 1.
    vector_of_int8,    //!< Used for "Integer" fields of size != 1 where each value fits in one byte.
    vector_of_int16,   //!< Used for "Integer" fields of size != 1 where each value fits in two bytes.
    vector_of_int32,   //!< Used for "Integer" fields of size != 1 where each value fits in four bytes.
    vector_of_float32, //!< Used for "Float" fields of size != 1.
    vector_of_string,  //!< Used for "String" fields of size != 1.
    flag               //!< Used for "Flat" fields (size must be 0).
};

} // namespace bio::io::var_io

namespace bio::io::var_io::detail
{

//!\brief int* and vector_of_int* are each "compatible" with each other; the rest only with self.
//!\ingroup var_io
constexpr bool type_id_is_compatible(var_io::value_type_id const lhs, var_io::value_type_id const rhs)
{
    switch (lhs)
    {
        case var_io::value_type_id::int8:
        case var_io::value_type_id::int16:
        case var_io::value_type_id::int32:
            switch (rhs)
            {
                case var_io::value_type_id::int8:
                case var_io::value_type_id::int16:
                case var_io::value_type_id::int32:
                    return true;
                default:
                    return false;
            };
            break;
        case var_io::value_type_id::vector_of_int8:
        case var_io::value_type_id::vector_of_int16:
        case var_io::value_type_id::vector_of_int32:
            switch (rhs)
            {
                case var_io::value_type_id::vector_of_int8:
                case var_io::value_type_id::vector_of_int16:
                case var_io::value_type_id::vector_of_int32:
                    return true;
                default:
                    return false;
            };
            break;
        default:
            return lhs == rhs;
    }
}

//-----------------------------------------------------------------------------
// bcf_type_descriptor and utilities
//-----------------------------------------------------------------------------

/*!\name BCF Type descriptor and utilities
 * \brief TODO some things here should get stand-alone tests maybe?
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
    //  int64
    float32 = 5,
    //  double
    char8   = 7,
};

/*!\addtogroup var_io
 * \{
 */

//!\brief Compute the smallest possible integral type descriptor able to represent the value.
detail::bcf_type_descriptor smallest_int_desc(std::unsigned_integral auto const num)
{
    // bits required to represent number (the +1 because signed integral has smaller range)
    switch (std::bit_ceil(std::bit_width(num) + 1))
    {
        case 128:
        case 64:
            throw bio_error{"Could not write number '", num, "'. Value out of range (only int32 supported)."};
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
    // If a value is the missing value (lowest), we can always encode it as the one-byte missing value
    return num == var_io::missing_value<decltype(num)> ? detail::bcf_type_descriptor::int8
                                                       : smallest_int_desc(static_cast<uint64_t>(std::abs(num)));
}

//!\overload
detail::bcf_type_descriptor smallest_int_desc(std::ranges::input_range auto && range)
{
    using val_t = bio::ranges::range_innermost_value_t<decltype(range)>;
    int64_t max = 0;

    // get max:
    // we don't use std::ranges::max_element here, so we can abort early if know we get max descriptor type
    for (int64_t elem : range)
    {
        if constexpr (std::signed_integral<val_t>)
            elem = static_cast<val_t>(elem) == var_io::missing_value<val_t> ? 0 : std::abs(elem);

        if (elem > max)
        {
            max = elem;
            if (max >= std::numeric_limits<int16_t>::max()) // this will always lead to bcf_type_descriptor::int32
                break;
        }
    }

    return smallest_int_desc(static_cast<uint64_t>(max));
}

//!\overload
template <std::ranges::forward_range rng_t>
    requires(std::ranges::range<std::ranges::range_reference_t<rng_t>>)
detail::bcf_type_descriptor smallest_int_desc(rng_t & range)
{
    return smallest_int_desc(range | std::views::join);
}

//!\brief Whether the value is any integer value.
inline bool type_descriptor_is_int(bcf_type_descriptor const type_desc)
{
    switch (type_desc)
    {
        case bcf_type_descriptor::int8:
        case bcf_type_descriptor::int16:
        case bcf_type_descriptor::int32:
            return true;
        default:
            return false;
    }
}

//!\brief Convert from bio::io::var_io::value_type_id to bio::io::detail::bcf_type_descriptor.
inline bcf_type_descriptor value_type_id_2_type_descriptor(var_io::value_type_id const type_id)
{
    switch (type_id)
    {
        case var_io::value_type_id::char8:
        case var_io::value_type_id::string:
        case var_io::value_type_id::vector_of_string:
            return bcf_type_descriptor::char8;
        case var_io::value_type_id::int8:
        case var_io::value_type_id::vector_of_int8:
        case var_io::value_type_id::flag:
            return bcf_type_descriptor::int8;
        case var_io::value_type_id::int16:
        case var_io::value_type_id::vector_of_int16:
            return bcf_type_descriptor::int16;
        case var_io::value_type_id::int32:
        case var_io::value_type_id::vector_of_int32:
            return bcf_type_descriptor::int32;
        case var_io::value_type_id::float32:
        case var_io::value_type_id::vector_of_float32:
            return bcf_type_descriptor::float32;
    }
    return bcf_type_descriptor::missing;
}

//!\}
//!\}

/*!\addtogroup var_io
 * \{
 */
/*!\name A compile-time mapping of types to bio::io::detail::bcf_type_descriptor
 * \{
 */
//!\brief Default implementation.
template <typename t>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor = bcf_type_descriptor::missing;

//!\brief Specialisation for int8.
template <std::signed_integral t>
    requires(sizeof(t) == 1)
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::int8;
//!\brief Specialisation for int8.
template <>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<bool> = bcf_type_descriptor::int8;

//!\brief Specialisation for int16.
template <std::signed_integral t>
    requires(sizeof(t) == 2)
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::int16;
//!\brief Specialisation for int16.
template <std::unsigned_integral t>
    requires(sizeof(t) == 1)
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::int16;

//!\brief Specialisation for int32.
template <std::signed_integral t>
    requires(sizeof(t) >= 4)
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
template <io::detail::deliberate_alphabet t>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::char8;

//!\brief Specialisation for cstring.
template <meta::decays_to<char const *> t>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> = bcf_type_descriptor::char8;

//!\brief Specialisation for range.
template <std::ranges::input_range t>
inline constexpr bcf_type_descriptor type_2_bcf_type_descriptor<t> =
  type_2_bcf_type_descriptor<bio::ranges::range_innermost_value_t<t>>;

//!\}

//!\brief Formula for computing indexes in genotype fields with number "G"; see VCF spec for details.
constexpr size_t vcf_gt_formula(size_t const a, size_t const b)
{
    return (b * (b + 1)) / 2 + a;
}

//!\}

//!\brief The field_ids used in this domain.
//!\ingroup var_io
static constexpr auto field_ids = meta::vtag<io::detail::field::chrom,
                                             io::detail::field::pos,
                                             io::detail::field::id,
                                             io::detail::field::ref,
                                             io::detail::field::alt,
                                             io::detail::field::qual,
                                             io::detail::field::filter,
                                             io::detail::field::info,
                                             io::detail::field::genotypes,
                                             io::detail::field::_private>;

} // namespace bio::io::var_io::detail

namespace bio::io::var_io
{

template <typename chrom_t,
          typename pos_t,
          typename id_t,
          typename ref_t,
          typename alt_t,
          typename qual_t,
          typename filter_t,
          typename info_t,
          typename genotypes_t>
struct record;

//!\brief Mixin for format handlers that helps converts to the tuple-record.
//!\ingroup var_io
struct format_handler_mixin
{
    //!\brief Convert from domain-specific record to tuple-record.
    template <typename... arg_ts>
    static auto record2tuple_record(var_io::record<arg_ts...> & in_record)
    {
        return io::detail::tie_tuple_record(detail::field_ids,
                                            in_record.chrom,
                                            in_record.pos,
                                            in_record.id,
                                            in_record.ref,
                                            in_record.alt,
                                            in_record.qual,
                                            in_record.filter,
                                            in_record.info,
                                            in_record.genotypes,
                                            in_record._private);
    }
};

} // namespace bio::io::var_io
