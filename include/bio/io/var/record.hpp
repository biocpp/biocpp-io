// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::var::record and auxilliary data structures.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include <bio/alphabet/custom/char.hpp>
#include <bio/alphabet/nucleotide/concept.hpp>
#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/meta/concept/core_language.hpp>
#include <bio/meta/tag/ttag.hpp>
#include <bio/meta/tag/vtag.hpp>
#include <bio/meta/tuple.hpp>
#include <bio/ranges/concept.hpp>
#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/container/dictionary.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/detail/range.hpp>
#include <bio/io/detail/tuple_record.hpp>
#include <bio/io/misc.hpp>
#include <bio/io/var/header.hpp>
#include <bio/io/var/misc.hpp>

//-----------------------------------------------------------------------------
// forwards
//-----------------------------------------------------------------------------

namespace bio::io::var
{

class header;

} // namespace bio::io::var

//-----------------------------------------------------------------------------
// BCF record core
//-----------------------------------------------------------------------------

namespace bio::io::var::detail
{

//!\brief The "core" of a BCF record in bit-compatible representation to the on-disk format.
//!\ingroup var
struct bcf_record_core
{
    int32_t  chrom         = -1;                                 //!< CHROM as IDX.
    int32_t  pos           = -1;                                 //!< POS.
    int32_t  rlen          = -1;                                 //!< Not used by this implementation.
    float    qual          = bio::io::var::missing_value<float>; //!< QUAL.
    uint16_t n_info        = 0;                                  //!< Number of INFOS values.
    uint16_t n_allele      = 0;                                  //!< Number of alleles.
    uint32_t n_sample : 24 = 0;                                  //!< Number of samples.
    uint8_t  n_fmt         = 0;                                  //!< Number of FORMAT values.
};

static_assert(sizeof(bcf_record_core) == 24, "Bit alignment problem in declaration of bcf_record_core.");

} // namespace bio::io::var::detail

//-----------------------------------------------------------------------------
// The info element
//-----------------------------------------------------------------------------

namespace bio::io::var
{

/*!\brief The base type of bio::io::var::info_variant_shallow .
 * \relates bio::io::var::info_variant_shallow
 */
using info_variant_shallow_base_t = std::variant<char,
                                                 int8_t,
                                                 int16_t,
                                                 int32_t,
                                                 float,
                                                 std::string_view,
                                                 std::vector<int8_t>,
                                                 std::vector<int16_t>,
                                                 std::vector<int32_t>,
                                                 std::vector<float>,
                                                 std::vector<std::string_view>,
                                                 bool>;

/*!\brief std::variant that stores the value of an INFO field [shallow version].
 * \ingroup var
 * \details
 *
 * This is a type for storing the value in an INFO field key-value pair.
 * Since these values can be of different types, this type is derived of std::variant
 * which allows storing values of different types ("variant" refers to the C++ type here, not
 * the biological meaning).
 * See bio::io::var::info_variant_shallow_base_t for the exact base type.
 *
 * To retrieve the contained value from a variable called `val`, use one of the following interfaces:
 *
 *   * `get<std::string_view>(val)` returns the contained value as a `std::string_view`.
 *   * `get<5>(val)` returns the contained value as a `std::string_view`; see
 * bio::io::var::info_variant_shallow_base_t regarding the order.
 *   * `get<"AA">(val)` returns the contained value as a `std::string_view`, because bio::io::var::info_key2type_enum associates that type with the key "AA".
 *
 * In all cases, an exception of std::bad_variant_access is thrown if the variant currently holds a
 * value of a different type.
 */
struct info_variant_shallow : info_variant_shallow_base_t
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    info_variant_shallow()                                         = default; //!< Defaulted.
    info_variant_shallow(info_variant_shallow &&)                  = default; //!< Defaulted.
    info_variant_shallow(info_variant_shallow const &)             = default; //!< Defaulted.
    info_variant_shallow & operator=(info_variant_shallow &&)      = default; //!< Defaulted.
    info_variant_shallow & operator=(info_variant_shallow const &) = default; //!< Defaulted.

    //!\brief Inherit base class's constructors.
    using info_variant_shallow_base_t::info_variant_shallow_base_t;
    //!\brief Inherit base class's assignment operator.
    using info_variant_shallow_base_t::operator=;
    //!\}

    /*!\brief Access the contained value by string.
     * \tparam key The string (literal).
     * \tparam me The variant parameter.
     * \throws std::bad_variant_access If the variant is in a different state.
     *
     * The mapping of string to type is defined by bio::io::var::info_key2type_enum.
     */
    template <ranges::small_string key>
    friend decltype(auto) get(meta::decays_to<info_variant_shallow> auto && me)
    {
        static_assert(std::same_as<decltype(info_key2type_enum<key>), value_type_id const>,
                      "No value_type_id found in bio::io::var::info_key2type_enum for this key.");

        return std::get<static_cast<size_t>(info_key2type_enum<key>)>(std::forward<decltype(me)>(me));
    }
};

/*!\brief The base type of bio::io::var::info_variant_deep .
 * \relates bio::io::var::info_variant_deep
 */
using info_variant_deep_base_t = std::variant<char,
                                              int8_t,
                                              int16_t,
                                              int32_t,
                                              float,
                                              std::string,
                                              std::vector<int8_t>,
                                              std::vector<int16_t>,
                                              std::vector<int32_t>,
                                              std::vector<float>,
                                              std::vector<std::string>,
                                              bool>;

/*!\brief std::variant that stores the value of an INFO field [deep version].
 * \ingroup var
 * \details
 *
 * This is a type for storing the value in an INFO field key-value pair.
 * Since these values can be of different types, this type is derived of std::variant
 * which allows storing values of different types ("variant" refers to the C++ type here, not
 * the biological meaning).
 * See bio::io::var::info_variant_shallow_base_t for the exact base type.
 *
 * To retrieve the contained value from a variable called `val`, use one of the following interfaces:
 *
 *   * `get<std::string>(val)` returns the contained value as a `std::string`.
 *   * `get<5>(val)` returns the contained value as a `std::string`; see
 * bio::io::var::info_variant_deep_base_t regarding the order.
 *   * `get<"AA">(val)` returns the contained value as a `std::string`, because bio::io::var::info_key2type_enum associates that type with the key "AA".
 *
 * In all cases, an exception of std::bad_variant_access is thrown if the variant currently holds a
 * value of a different type.
 */

struct info_variant_deep : info_variant_deep_base_t
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    info_variant_deep()                                      = default; //!< Defaulted.
    info_variant_deep(info_variant_deep &&)                  = default; //!< Defaulted.
    info_variant_deep(info_variant_deep const &)             = default; //!< Defaulted.
    info_variant_deep & operator=(info_variant_deep &&)      = default; //!< Defaulted.
    info_variant_deep & operator=(info_variant_deep const &) = default; //!< Defaulted.

    //!\brief Inherit base class's constructors.
    using info_variant_deep_base_t::info_variant_deep_base_t;
    //!\brief Inherit base class's assignment operator.
    using info_variant_deep_base_t::operator=;
    //!\}

    /*!\brief Access the contained value by key string.
     * \tparam key The string (literal).
     * \tparam me The variant parameter.
     * \throws std::bad_variant_access If the variant is in a different state.
     *
     * The mapping of string to type is defined by bio::io::var::info_key2type_enum.
     */
    template <ranges::small_string key>
    friend decltype(auto) get(meta::decays_to<info_variant_deep> auto && me)
    {
        static_assert(std::same_as<decltype(info_key2type_enum<key>), value_type_id const>,
                      "No value_type_id found in bio::io::var::info_key2type_enum for this key.");

        return std::get<static_cast<size_t>(info_key2type_enum<key>)>(std::forward<decltype(me)>(me));
    }
};

} // namespace bio::io::var

namespace std
{

//!\cond
template <>
struct variant_size<bio::io::var::info_variant_shallow> : variant_size<bio::io::var::info_variant_shallow_base_t>
{};

template <>
struct variant_size<bio::io::var::info_variant_deep> : variant_size<bio::io::var::info_variant_deep_base_t>
{};
//!\endcond
} // namespace std

namespace bio::io::var::detail
{
//!\brief Auxilliary concept that encompasses bio::io::var::info_variant.
//!\ingroup var
template <typename t>
concept is_info_variant = meta::one_of<t, var::info_variant_shallow, var::info_variant_deep>;

} // namespace bio::io::var::detail

namespace bio::io::var
{

//-----------------------------------------------------------------------------
// The genotype element
//-----------------------------------------------------------------------------

/*!\brief The base type of bio::io::var::genotype_variant_shallow.
 * \relates bio::io::var::genotype_variant_shallow
 */
using genotype_variant_shallow_base_t = std::variant<std::vector<char>,
                                                     std::vector<int8_t>,
                                                     std::vector<int16_t>,
                                                     std::vector<int32_t>,
                                                     std::vector<float>,
                                                     std::vector<std::string_view>,
                                                     ranges::concatenated_sequences<std::vector<int8_t>>,
                                                     ranges::concatenated_sequences<std::vector<int16_t>>,
                                                     ranges::concatenated_sequences<std::vector<int32_t>>,
                                                     ranges::concatenated_sequences<std::vector<float>>,
                                                     std::vector<std::vector<std::string>>
                                                     /* no flag here */>;

/*!\brief std::variant that stores the value of a GENOTYPE field [shallow version].
 * \ingroup var
 * \details
 *
 * This is a type for storing the value in a GENOTYPE field key-value pair.
 * Since the value in such field can have different types (depending on the field),
 * this type is derived of std::variant
 * which allows storing values of different types ("variant" refers to the C++ type here, not
 * the biological meaning).
 * See bio::io::var::genotype_variant_shallow for the exact base type.
 *
 * Note that this follows the BCF representation where data is grouped "by-field" and not "by-sample",
 * e.g. the "GT" values of all samples are in one vector.
 * This also means that a bio::io::var::value_type_id::string implies a vector-of-strings and
 * bio::io::var::value_type_id::vector_of_string implies a vector-of-vector-of-strings.
 * **All possible types for this variant are a container** that is either the same size as the number of
 * samples or empty.
 *
 * To retrieve the contained value from a variable called `val`, use one of the following interfaces:
 *
 *   * `get<std::vector<std::string_view>>(val)` returns the contained value as a `std::vector<std::string_view>`.
 *   * `get<5>(val)` returns the contained value as a `std::vector<std::string_view>`; see
 * bio::io::var::genotype_variant_shallow regarding the order.
 *   * `get<"GT">(val)` returns the contained value as a `std::vector<std::string_view>`, because bio::io::var::info_key2type_enum associates that type with the key "AA".
 *
 * In all cases, an exception of std::bad_variant_access is thrown if the variant currently holds a
 * value of a different type.
 */
struct genotype_variant_shallow : genotype_variant_shallow_base_t
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    genotype_variant_shallow()                                             = default; //!< Defaulted.
    genotype_variant_shallow(genotype_variant_shallow &&)                  = default; //!< Defaulted.
    genotype_variant_shallow(genotype_variant_shallow const &)             = default; //!< Defaulted.
    genotype_variant_shallow & operator=(genotype_variant_shallow &&)      = default; //!< Defaulted.
    genotype_variant_shallow & operator=(genotype_variant_shallow const &) = default; //!< Defaulted.

    //!\brief Inherit base class's constructors.
    using genotype_variant_shallow_base_t::genotype_variant_shallow_base_t;
    //!\brief Inherit base class's assignment operator.
    using genotype_variant_shallow_base_t::operator=;
    //!\}

    /*!\brief Access the contained value by key string.
     * \tparam key The string (literal).
     * \tparam me The variant parameter.
     * \throws std::bad_variant_access If the variant is in a different state.
     *
     * The mapping of string to type is defined by bio::io::var::genotype_key2type_enum.
     */
    template <ranges::small_string key>
    friend decltype(auto) get(meta::decays_to<genotype_variant_shallow> auto && me)
    {
        static_assert(std::same_as<decltype(format_key2type_enum<key>), value_type_id const>,
                      "No value_type_id found in bio::io::var::format_key2type_enum for this key.");

        return std::get<static_cast<size_t>(format_key2type_enum<key>)>(std::forward<decltype(me)>(me));
    }
};

/*!\brief The base type of bio::io::var::genotype_variant_deep.
 * \relates bio::io::var::genotype_variant_deep
 */
using genotype_variant_deep_base_t = std::variant<std::vector<char>,
                                                  std::vector<int8_t>,
                                                  std::vector<int16_t>,
                                                  std::vector<int32_t>,
                                                  std::vector<float>,
                                                  std::vector<std::string>,
                                                  ranges::concatenated_sequences<std::vector<int8_t>>,
                                                  ranges::concatenated_sequences<std::vector<int16_t>>,
                                                  ranges::concatenated_sequences<std::vector<int32_t>>,
                                                  ranges::concatenated_sequences<std::vector<float>>,
                                                  std::vector<std::vector<std::string>>
                                                  /* no flag here */>;

/*!\brief std::variant that stores the value of a GENOTYPE field [deep version].
 * \ingroup var
 * \details
 *
 * This is a type for storing the value in a GENOTYPE field key-value pair.
 * Since the value in such field can have different types (depending on the field),
 * this type is derived of std::variant
 * which allows storing values of different types ("variant" refers to the C++ type here, not
 * the biological meaning).
 * See bio::io::var::genotype_variant_deep for the exact base type.
 *
 * Note that this follows the BCF representation where data is grouped "by-field" and not "by-sample",
 * e.g. the "GT" values of all samples are in one vector.
 * This also means that a bio::io::var::value_type_id::string implies a vector-of-strings and
 * bio::io::var::value_type_id::vector_of_string implies a vector-of-vector-of-strings.
 * **All possible types for this variant are a container** that is either the same size as the number of
 * samples or empty.
 *
 * To retrieve the contained value from a variable called `val`, use one of the following interfaces:
 *
 *   * `get<std::vector<std::string>>(val)` returns the contained value as a `std::vector<std::string>`.
 *   * `get<5>(val)` returns the contained value as a `std::vector<std::string>`; see
 * bio::io::var::genotype_variant_deep regarding the order.
 *   * `get<"GT">(val)` returns the contained value as a `std::vector<std::string>`, because bio::io::var::info_key2type_enum associates that type with the key "AA".
 *
 * In all cases, an exception of std::bad_variant_access is thrown if the variant currently holds a
 * value of a different type.
 */
struct genotype_variant_deep : genotype_variant_deep_base_t
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    genotype_variant_deep()                                          = default; //!< Defaulted.
    genotype_variant_deep(genotype_variant_deep &&)                  = default; //!< Defaulted.
    genotype_variant_deep(genotype_variant_deep const &)             = default; //!< Defaulted.
    genotype_variant_deep & operator=(genotype_variant_deep &&)      = default; //!< Defaulted.
    genotype_variant_deep & operator=(genotype_variant_deep const &) = default; //!< Defaulted.

    //!\brief Inherit base class's constructors.
    using genotype_variant_deep_base_t::genotype_variant_deep_base_t;
    //!\brief Inherit base class's assignment operator.
    using genotype_variant_deep_base_t::operator=;
    //!\}

    /*!\brief Access the contained value by key string.
     * \tparam key The string (literal).
     * \tparam me The variant parameter.
     * \throws std::bad_variant_access If the variant is in a different state.
     *
     * The mapping of string to type is defined by bio::io::var::genotype_key2type_enum.
     */
    template <ranges::small_string key>
    friend decltype(auto) get(meta::decays_to<genotype_variant_deep> auto && me)
    {
        static_assert(std::same_as<decltype(format_key2type_enum<key>), value_type_id const>,
                      "No value_type_id found in bio::io::var::format_key2type_enum for this key.");

        return std::get<static_cast<size_t>(format_key2type_enum<key>)>(std::forward<decltype(me)>(me));
    }
};

} // namespace bio::io::var

namespace std
{

//!\cond
template <>
struct variant_size<bio::io::var::genotype_variant_shallow> :
  variant_size<bio::io::var::genotype_variant_shallow_base_t>
{};

template <>
struct variant_size<bio::io::var::genotype_variant_deep> : variant_size<bio::io::var::genotype_variant_deep_base_t>
{};
//!\endcond

} // namespace std

namespace bio::io::var::detail
{

//!\brief Auxilliary concept that encompasses bio::io::var::genotype_variant.
//!\ingroup var
template <typename t>
concept is_genotype_variant = meta::one_of<t, var::genotype_variant_shallow, var::genotype_variant_deep>;

} // namespace bio::io::var::detail

namespace bio::io::var
{

//-----------------------------------------------------------------------------
// record_private_data
//-----------------------------------------------------------------------------

//!\brief A datastructure that contains private data of variant IO records.
//!\ingroup var
struct record_private_data
{
    //!\privatesection
    //!\brief Pointer to the header
    header const * header_ptr = nullptr;

    //!\brief Pointer to record core (if BCF).
    detail::bcf_record_core const * record_core = nullptr;

    //!\brief Raw record type.
    using raw_record_t =
      io::detail::tuple_record<meta::vtag_t<io::detail::field::chrom,
                                            io::detail::field::pos,
                                            io::detail::field::id,
                                            io::detail::field::ref,
                                            io::detail::field::alt,
                                            io::detail::field::qual,
                                            io::detail::field::filter,
                                            io::detail::field::info,
                                            io::detail::field::genotypes,
                                            io::detail::field::_private>,
                               meta::list_traits::concat<meta::list_traits::repeat<9, std::span<std::byte const>>,
                                                         meta::type_list<var::record_private_data>>>;
    //!\brief Pointer to raw record.
    raw_record_t const * raw_record = nullptr;

    //!\brief Defaulted three-way comparison.
    friend bool operator==(record_private_data const &, record_private_data const &) = default;
};

//-----------------------------------------------------------------------------
// The record type
//-----------------------------------------------------------------------------

/*!\brief Record type for variant I/O.
 * \ingroup var
 * \tparam chrom_t     Type of the CHROM member. See the member for type requirements.
 * \tparam pos_t       Type of the POS member. See the member for type requirements.
 * \tparam id_t        Type of the ID member. See the member for type requirements.
 * \tparam ref_t       Type of the REF member. See the member for type requirements.
 * \tparam alt_t       Type of the ALT member. See the member for type requirements.
 * \tparam qual_t      Type of the QUAL member. See the member for type requirements.
 * \tparam filter_t    Type of the FILTER member. See the member for type requirements.
 * \tparam info_t      Type of the INFO member. See the member for type requirements.
 * \tparam genotypes_t Type of the GENOTYPES member. See the member for type requirements.
 * \details
 *
 * This is the record template for variant I/O.
 * The internal representation of VCF and BCF are different. To be able to freely
 * interchange between these formats, this library needs to choose one representation that
 * everything is converted to when being read.
 *
 * One thing that is fixed for all configurations in this library is the layout of the GENOTYPES field
 * which is always grouped "by-field" (BCF-style) and not "by-sample" (VCF-style).
 * Another important choice is that **numbers are always 1-based,** because this is the default in VCF
 * and all other tools that deal with VCF/BCF.
 *
 * Some other things about the reading/writing behaviour can be tweaked by configuring this record
 * to store different types than the default.
 * However, most users will be happy with one of the pre-defined aliases.
 *
 * ### Pre-defined aliases
 *
 * * bio::io::var::record_deep:
 *   * the record is **deep**
 *   * CHROM, FILTER, INFO element IDs and GENOTYPE element IDs are **strings**
 * * bio::io::var::record_shallow:
 *   * the record is **shallow**
 *   * CHROM, FILTER, INFO element IDs and GENOTYPE element IDs are **string_views**
 * * bio::io::var::record_idx:
 *   * the record is **deep**
 *   * CHROM, FILTER, INFO element IDs and GENOTYPE element IDs are IDX values (int32_t).
 * * bio::io::var::record_idx_shallow:
 *   * the record is **shallow**
 *   * CHROM, FILTER, INFO element IDs and GENOTYPE element IDs are IDX values (int32_t).
 *
 * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
 *
 * ### Strings VS IDX
 *
 * By default, all string-like fields are always read/written as strings or string_views.
 * However, the internal representation of BCF uses numeric IDs, so called IDX values,
 * instead.
 * If an application predominanly reads/writes BCF and the headers of input and output
 * files are not modified, it may be faster to use and IDX-based record.
 *
 * ### Manually specifying the memeber types
 *
 * The types for members can be set individually. Among other things, this allows
 * disabling the parsing of unwanted fields:
 *
 * \snippet test/snippet/var/var_reader_options.cpp field_types_expert
 *
 * Designated intitialisers are used here to simplify the
 * definition of the type.
 *
 *
 * See the \ref record_faq for more information on record-based reading.
 */
template <typename _chrom_t     = std::string,
          typename _pos_t       = int32_t,
          typename _id_t        = std::string,
          typename _ref_t       = std::vector<alphabet::dna5>,
          typename _alt_t       = std::vector<std::string>,
          typename _qual_t      = float,
          typename _filter_t    = std::vector<std::string>,
          typename _info_t      = ranges::dictionary<std::string, info_variant_deep>,
          typename _genotypes_t = ranges::dictionary<std::string, genotype_variant_deep>>
struct record
{
    using chrom_t     = _chrom_t;     //!< Type of the chrom member.
    using pos_t       = _pos_t;       //!< Type of the pos member.
    using id_t        = _id_t;        //!< Type of the id member.
    using ref_t       = _ref_t;       //!< Type of the ref member.
    using alt_t       = _alt_t;       //!< Type of the alt member.
    using qual_t      = _qual_t;      //!< Type of the qual member.
    using filter_t    = _filter_t;    //!< Type of the filter member.
    using info_t      = _info_t;      //!< Type of the info member.
    using genotypes_t = _genotypes_t; //!< Type of the genotypes member.

    /*!\brief The CHROM field – an identifier from the reference genome.
     * \details
     *
     * This member can be represented either as a string or an IDX value.
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     * 1. any back-insertable range over the `char` alphabet, e.g. a std::string (**deep**, string)
     * 2. std::string_view (**shallow**, string)
     * 3. `int32_t` (IDX value)
     * 4. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over `char` (string).
     * 2. `int32_t` (IDX value)
     *
     * It may not be ignored when writing.
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     *
     */
    chrom_t chrom{};

    /*!\brief The POS field – the reference position (1-based).
     * \details
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     * 1. any std::integral type; int32_t recommended!
     * 2. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any std::integral type
     *
     * It may not be ignored when writing.
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     */
    pos_t pos{};

    /*!\brief The ID field – a unique identifier.
     * \details
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     * 1. any back-insertable range over the `char` alphabet, e.g. a std::string (**deep**)
     * 2. std::string_view (**shallow**)
     * 3. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over `char`.
     * 2. bio::meta::ignore_t
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     */
    id_t id{};

    /*!\brief The REF field – reference bases.
     * \details
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over the `char` alphabet (**deep**)
     *   2. any back-insertable range over a bio::alphabet::nucleotide (**deep**, converted elements)
     *   3. std::string_view (**shallow**)
     *   4. a std::ranges::transform_view defined on a std::string_view (**shallow**, converted elements)
     *   5. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over `char` or a bio::alphabet::nucleotide.
     *
     * It may not be ignored when writing.
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     *
     */
    ref_t ref{};

    /*!\brief The ALT field – a vector of alternative allele bases.
     * \details
     *
     * This field is always encoded as a range/container of multiple alleles.
     * If there are no alleles in the file (denoted by a single '.'), this
     * range will be empty.
     *
     * Since breakpoint strings are legal, the element type is typically
     * std::string / std::string_view and not a vector/view of bio::alphabet::dna5.
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over a back-insertable of `char` (**deep**)
     *   2. any back-insertable range over a back-insertable of bio::alphabet::nucleotide (**deep**, converted elements)
     *   3. any back-insertable range over std::string_view (**shallow**)
     *   4. any back-insertable range over std::ranges::transform_view defined on a std::string_view (**shallow**,
     * converted elements)
     *   5. bio::meta::ignore_t (**ignored**)
     *
     * For 2. and 4., an exception will be thrown if a breakpoint string is encountered.
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over a std::ranges::forward_range of `char`
     * or a bio::alphabet::nucleotide.
     * 2. bio::meta::ignore_t
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     *
     */
    alt_t alt{};

    /*!\brief The QUAL field – phred-scaled quality score for the assertion made in ALT.
     * \details
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     * 1. any bio::meta::arithmetic type
     * 2. bio::meta::ignore_t (**ignored**)
     *
     * Note that '.' will be read as the respective bio::io::var::missing_value.
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any bio::meta::arithmetic type
     * 2. bio::meta::ignore_t
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     */
    qual_t qual{};

    /*!\brief The FILTER field – a vector of FILTER values.
     * \details
     *
     * This field is always encoded as a range/container of multiple filters.
     * If there are no filters given in the file (denoted by a single '.'), this
     * range will be empty.
     *
     * The individual filter values may be encoded as either string/string_views or
     * IDX values.
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over a back-insertable range of `char` (**deep**, string)
     *   2. any back-insertable range over std::string_view (**shallow**, string)
     *   3. any back-insertable range over int32_t (IDX values)
     *   4. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over a std::ranges::forward_range of `char` (interpreted as strings)
     * 2. any std::ranges::forward_range over int32_t (interpreted as IDX values)
     * 2. bio::meta::ignore_t
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     *
     */
    filter_t filter{};

    /*!\brief The INFO field – a vector of additional information values.
     * \details
     *
     * This field is always encoded as a range/container of "info elements".
     * If there are no information values given in the file (denoted by a single '.'), this
     * range will be empty.
     *
     * "Info elements" are a pair/tuple of size 2, where the first element is a key (e.g. "AD") and
     * the second element is a value (e.g. "20,3").
     * Since fields have values of different types, the value is usually encoded in a std::variant-like
     * type: bio::io::var::info_variant_deep or bio::io::var::info_variant_shallow.
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over a bio::meta::tuple (or std::tuple or std::pair) of
     *      1. std::string (**deep**) or std::string_view (**shallow**) or int32_t (IDX value)
     *      2. bio::io::var::info_variant_deep or bio::io::var::info_variant_shallow
     *   2. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     *   1. any std::ranges::forward_range over a bio::meta::tuple (or std::tuple or std::pair) of
     *      1. std::string (**deep**) or std::string_view (**shallow**) or int32_t (IDX value)
     *      2. bio::io::var::info_variant_deep or bio::io::var::info_variant_shallow
     *   2. bio::meta::ignore_t (**ignored**)
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     *
     */
    info_t info{};

    /*!\brief The GENOTYPES field – a vector of genotype values.
     * \details
     *
     * This field is always encoded as a range/container of "genotype elements".
     * If there are no information values given in the file (no sample columns), this
     * range will be empty.
     *
     * NOTE that this field is encoded "by genotype" (BCF-style) and not
     * "by sample" (VCF-style).
     * Every element represents one piece of genotyping information (e.g. GT or
     * AD), and contains a vector of values for each sample.
     *
     * "Genotype elements" are a pair/tuple size of 2, where the first element is a key (e.g. "GT") and
     * the second element is a value (e.g. "0/1;1/1;0/1;...").
     * Since fields have values of different types, the value is usually encoded in a std::variant-like
     * type: bio::io::var::genotype_variant_deep or bio::io::var::genotype_variant_shallow.
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over a bio::meta::tuple (or std::tuple or std::pair) of
     *      1. std::string (**deep**) or std::string_view (**shallow**) or int32_t (IDX value)
     *      2. bio::io::var::genotype_variant_deep or bio::io::var::genotype_variant_shallow
     *   2. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq::writer), the type can be one of the following:
     *
     *   1. any std::ranges::forward_range over a bio::meta::tuple (or std::tuple or std::pair) of
     *      1. std::string (**deep**) or std::string_view (**shallow**) or int32_t (IDX value)
     *      2. bio::io::var::genotype_variant_deep or bio::io::var::genotype_variant_shallow
     *   2. bio::meta::ignore_t (**ignored**)
     *
     * The default and all pre-defined aliases satisfy the requirements for reading and writing.
     *
     */
    genotypes_t genotypes{};

    //!\brief A member that holds internal use data structures.
    //!\private
    record_private_data _private{};

    //!\brief Clear all members.
    void clear()
    {
        io::detail::clear(chrom);
        io::detail::clear(pos);
        io::detail::clear(id);
        io::detail::clear(ref);
        io::detail::clear(alt);
        io::detail::clear(qual);
        io::detail::clear(filter);
        io::detail::clear(info);
        io::detail::clear(genotypes);
        io::detail::clear(_private);
    }

    //!\brief Defaulted comparison operators.
    friend auto operator<=>(record const & lhs, record const & rhs) = default;
};

//-----------------------------------------------------------------------------
// tie_record()
//-----------------------------------------------------------------------------

/*!\brief Helper function for easily creating a record of references.
 * \ingroup var
 * \tparam chrom_t     Type of the CHROM parameter.
 * \tparam pos_t       Type of the POS parameter.
 * \tparam id_t        Type of the ID parameter.
 * \tparam ref_t       Type of the REF parameter.
 * \tparam alt_t       Type of the ALT parameter.
 * \tparam qual_t      Type of the QUAL parameter.
 * \tparam filter_t    Type of the FILTER parameter.
 * \tparam info_t      Type of the INFO parameter.
 * \tparam genotypes_t Type of the GENOTYPES parameter.
 * \param[in,out] chrom     The CHROM parameter.
 * \param[in,out] pos       The POS parameter.
 * \param[in,out] id        The ID parameter.
 * \param[in,out] ref       The REF parameter.
 * \param[in,out] alt       The ALT parameter.
 * \param[in,out] qual      The QUAL parameter.
 * \param[in,out] filter    The FILTER parameter.
 * \param[in,out] info      The INFO parameter.
 * \param[in,out] genotypes The GENOTYPES parameter.
 * \returns A bio::io::var::record where all the members are references to this function's parameters.
 * \details
 *
 * ### Example
 *
 * TODO
 */
template <typename chrom_t,
          typename pos_t,
          typename id_t,
          typename ref_t,
          typename alt_t,
          typename qual_t,
          typename filter_t,
          typename info_t,
          typename genotypes_t>
auto tie_record(chrom_t &     chrom,
                pos_t &       pos,
                id_t &        id,
                ref_t &       ref,
                alt_t &       alt,
                qual_t &      qual,
                filter_t &    filter,
                info_t &      info,
                genotypes_t & genotypes)
{
    using r_t = record<chrom_t &, pos_t &, id_t &, ref_t &, alt_t &, qual_t &, filter_t &, info_t &, genotypes_t &>;
    return r_t{chrom, pos, id, ref, alt, qual, filter, info, genotypes};
}

//-----------------------------------------------------------------------------
// Pre-defined record aliases
//-----------------------------------------------------------------------------

/*!\brief The deep version of the default record type.
 *!\ingroup var
 * \see bio::io::var::record
 */
using record_deep = record<std::string,                                             // chrom,
                           int32_t,                                                 // pos,
                           std::string,                                             // id,
                           std::vector<alphabet::dna5>,                             // ref,
                           std::vector<std::string>,                                // alt,
                           float,                                                   // qual,
                           std::vector<std::string>,                                // filter,
                           ranges::dictionary<std::string, info_variant_deep>,      // info,
                           ranges::dictionary<std::string, genotype_variant_deep>>; // genotypes,

/*!\brief The record type used by bio::io::var::reader by default.
 *!\ingroup var
 * \see bio::io::var::record
 */
using record_shallow = record<std::string_view,                                                // chrom,
                              int32_t,                                                         // pos,
                              std::string_view,                                                // id,
                              views::char_conversion_view_t<alphabet::dna5>,                   // ref,
                              std::vector<std::string_view>,                                   // alt,
                              float,                                                           // qual,
                              std::vector<std::string_view>,                                   // filter,
                              ranges::dictionary<std::string_view, info_variant_shallow>,      // info,
                              ranges::dictionary<std::string_view, genotype_variant_shallow>>; // genotypes,

/*!\brief A record type with IDX values (shallow).
 *!\ingroup var
 * \see bio::io::var::record
 */
using record_idx_shallow = record<int32_t,                                                    // chrom,
                                  int32_t,                                                    // pos,
                                  std::string_view,                                           // id,
                                  views::char_conversion_view_t<alphabet::dna5>,              // ref,
                                  std::vector<std::string_view>,                              // alt,
                                  float,                                                      // qual,
                                  std::vector<int32_t>,                                       // filter,
                                  std::vector<std::pair<int32_t, info_variant_shallow>>,      // info,
                                  std::vector<std::pair<int32_t, genotype_variant_shallow>>>; // genotypes,

/*!\brief A record type with IDX values (deep).
 *!\ingroup var
 * \see bio::io::var::record
 */
using record_idx_deep = record<int32_t,                                                 // chrom,
                               int32_t,                                                 // pos,
                               std::string,                                             // id,
                               std::vector<alphabet::dna5>,                             // ref,
                               std::vector<std::string>,                                // alt,
                               float,                                                   // qual,
                               std::vector<int32_t>,                                    // filter,
                               std::vector<std::pair<int32_t, info_variant_deep>>,      // info,
                               std::vector<std::pair<int32_t, genotype_variant_deep>>>; // genotypes,

} // namespace bio::io::var

namespace bio::io::var::detail
{

//-----------------------------------------------------------------------------
// Record concept checks for reading
//-----------------------------------------------------------------------------

template <typename t>
concept pair_like = requires
{
    typename std::tuple_size<t>::type;
    requires std::tuple_size_v<t>
    == 2;
    typename std::tuple_element<0, t>::type;
    typename std::tuple_element<1, t>::type;
};

/*!\interface bio::io::detail::info_element_reader_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var::info_element / bio::io::var::info_element_idx.
 */
//!\cond CONCEPT_DEF
template <typename t, typename noref_t = std::remove_reference_t<t>>
concept info_element_reader_concept = pair_like<noref_t> &&
  (io::detail::out_string<std::tuple_element_t<0, noref_t>> ||
   std::same_as<int32_t &, std::tuple_element_t<0, noref_t> &>)&&detail::
    is_info_variant<std::remove_reference_t<std::tuple_element_t<1, noref_t>>>;
//!\endcond

/*!\interface bio::io::detail::genotype_reader_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var::genotype_element / bio::io::var::genotype_element_idx.
 */
//!\cond CONCEPT_DEF
template <typename t, typename noref_t = std::remove_reference_t<t>>
concept genotype_reader_concept = pair_like<noref_t> &&
  (io::detail::out_string<std::tuple_element_t<0, noref_t>> ||
   std::same_as<int32_t &, std::tuple_element_t<0, noref_t> &>)&&detail::
    is_genotype_variant<std::remove_reference_t<std::tuple_element_t<1, noref_t>>>;
//!\endcond

//!\brief Validates the concepts that the record type needs to satisfy when being passed to a reader.
template <typename chrom_t,
          typename pos_t,
          typename id_t,
          typename ref_t,
          typename alt_t,
          typename qual_t,
          typename filter_t,
          typename info_t,
          typename genotypes_t>
constexpr bool record_read_concept_checker(
  std::type_identity<var::record<chrom_t, pos_t, id_t, ref_t, alt_t, qual_t, filter_t, info_t, genotypes_t>>)

{
    // TODO(GCC11): once GCC10 is dropped, remove the "<typename t = seq_t>"
    static_assert(
      ranges::back_insertable_with<chrom_t, char> ||
        meta::one_of<std::remove_reference_t<chrom_t>, std::string_view, int32_t, meta::ignore_t, meta::ignore_t const>,
      "Requirements for the field-type of the CHROM-field not met. See documentation for "
      "bio::io::var::record.");

    static_assert(std::integral<std::remove_reference_t<pos_t>> || meta::decays_to<pos_t, meta::ignore_t>,
                  "Requirements for the field-type of the POS-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(ranges::back_insertable_with<id_t, char> ||
                    meta::one_of<std::remove_reference_t<id_t>, std::string_view, meta::ignore_t, meta::ignore_t const>,
                  "Requirements for the field-type of the ID-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(io::detail::lazy_concept_checker([]<typename t = ref_t>(auto) requires(
                    (ranges::back_insertable<t> && alphabet::alphabet<std::ranges::range_reference_t<t>>) ||
                    meta::one_of<std::remove_reference_t<t>, std::string_view, meta::ignore_t, meta::ignore_t const> ||
                    io::detail::transform_view_on_string_view<t>) { return std::true_type{}; }),
                  "Requirements for the field-type of the REF-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(io::detail::lazy_concept_checker([]<typename t = alt_t>(auto) requires(
                    meta::decays_to<t, meta::ignore_t> ||
                    (ranges::back_insertable<t> &&
                     ((ranges::back_insertable<std::ranges::range_reference_t<t>> &&
                       alphabet::alphabet<std::ranges::range_reference_t<std::ranges::range_reference_t<t>>>) ||
                      meta::decays_to<std::ranges::range_reference_t<t>, std::string_view> ||
                      io::detail::transform_view_on_string_view<std::ranges::range_reference_t<t>>))) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the ALT-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(meta::arithmetic<std::remove_reference_t<qual_t>> || meta::decays_to<qual_t, meta::ignore_t>,
                  "Requirements for the field-type of the QUAL-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(
      io::detail::lazy_concept_checker([]<typename t = filter_t>(auto) requires(
        meta::decays_to<t, meta::ignore_t> ||
        (ranges::back_insertable<t> &&
         (ranges::back_insertable_with<std::ranges::range_reference_t<t>, char> ||
          meta::one_of<std::remove_reference_t<std::ranges::range_reference_t<t>>, std::string_view, int32_t>))) {
          return std::true_type{};
      }),
      "Requirements for the field-type of the FILTER-field not met. See documentation for "
      "bio::io::var::record.");

    static_assert(io::detail::lazy_concept_checker([]<typename t = info_t>(auto) requires(
                    meta::decays_to<t, meta::ignore_t> ||
                    (ranges::back_insertable<t> &&
                     detail::info_element_reader_concept<std::ranges::range_value_t<t>>)) { return std::true_type{}; }),
                  "Requirements for the field-type of the INFO-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(io::detail::lazy_concept_checker([]<typename t = genotypes_t>(auto) requires(
                    meta::decays_to<t, meta::ignore_t> ||
                    (ranges::back_insertable<t> && detail::genotype_reader_concept<std::ranges::range_value_t<t>>)) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the GENOTYPES-field not met. See documentation for "
                  "bio::io::var::record.");

    return true;
}

//-----------------------------------------------------------------------------
// Record concept checks for writing
//-----------------------------------------------------------------------------

//!\brief A char, signed_integral, floating point number or CString.
template <typename t>
concept legal_type_aux =
  std::same_as<t, char> || std::signed_integral<t> || std::floating_point<t> || meta::decays_to<t, char const *>;

/*!\interface bio::io::var::detail::legal_type <>
 * \tparam t The type to check.
 * \brief A type that is similar to one of the alternatives of bio::io::var::info_variant
 */
//!\cond CONCEPT_DEF
template <typename t>
concept legal_type = legal_type_aux<std::remove_cvref_t<t>> || std::same_as<t const &, bool const &> ||
  (std::ranges::forward_range<t> && (legal_type_aux<std::remove_cvref_t<std::ranges::range_reference_t<t>>> ||
                                     io::detail::char_range<std::ranges::range_reference_t<t>>));
//!\endcond

/*!\interface bio::io::var::detail::legal_vector_type <>
 * \tparam t The type to check.
 * \brief A type that is similar to one of the alternatives of bio::io::var::info_variant
 */
//!\cond CONCEPT_DEF
template <typename t>
concept legal_vector_type = std::ranges::forward_range<t> && legal_type<std::ranges::range_reference_t<t>> &&
  !std::same_as<bool const &, std::ranges::range_reference_t<t>>;
//!\endcond

/*!\interface bio::io::var::detail::legal_or_dynamic <>
 * \tparam t The type to check.
 * \brief A type that is similar to one of the alternatives of bio::io::var::info_variant
 */
//!\cond CONCEPT_DEF
template <typename t>
concept legal_or_dynamic = legal_type<t> || is_info_variant<t>;
//!\endcond

/*!\interface bio::io::var::detail::vector_legal_or_dynamic <>
 * \tparam t The type to check.
 * \brief A type that is similar to one of the alternatives of bio::io::var::info_variant
 */
//!\cond CONCEPT_DEF
template <typename t>
concept vector_legal_or_dynamic = legal_vector_type<t> || is_genotype_variant<t>;
//!\endcond

/*!\interface bio::io::var::detail::info_element_writer_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var::info_element / bio::io::var::info_element_idx.
 */
//!\cond CONCEPT_DEF
template <typename t, typename noref_t = std::remove_reference_t<t>>
concept info_element_writer_concept = pair_like<noref_t> &&
  (io::detail::char_range_or_cstring<std::remove_cvref_t<std::tuple_element_t<0, noref_t>>> ||
   std::same_as<int32_t, std::remove_cvref_t<std::tuple_element_t<0, noref_t>>>)&&detail::
    legal_or_dynamic<std::remove_cvref_t<std::tuple_element_t<1, noref_t>>>;
//!\endcond

/*!\interface bio::io::var::detail::genotype_writer_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var::genotype_element / bio::io::var::genotype_element_idx.
 */
//!\cond CONCEPT_DEF
template <typename t, typename noref_t = std::remove_reference_t<t>>
concept genotype_writer_concept = pair_like<noref_t> &&
  (io::detail::char_range_or_cstring<std::remove_cvref_t<std::tuple_element_t<0, noref_t>>> ||
   std::same_as<int32_t, std::remove_cvref_t<std::tuple_element_t<0, noref_t>>>)&&detail::
    vector_legal_or_dynamic<std::remove_cvref_t<std::tuple_element_t<1, noref_t>>>;
//!\endcond

//!\brief Validates the concepts that the record type needs to satisfy when being passed to a writer.
template <typename chrom_t,
          typename pos_t,
          typename id_t,
          typename ref_t,
          typename alt_t,
          typename qual_t,
          typename filter_t,
          typename info_t,
          typename genotypes_t>
constexpr bool record_write_concept_checker(
  std::type_identity<var::record<chrom_t, pos_t, id_t, ref_t, alt_t, qual_t, filter_t, info_t, genotypes_t>>)

{
    static_assert(io::detail::char_range<chrom_t> || meta::decays_to<chrom_t, int32_t>,
                  "Requirements for the field-type of the CHROM-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(std::integral<std::remove_cvref_t<pos_t>>,
                  "Requirements for the field-type of the POS-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(io::detail::char_range<id_t> || meta::decays_to<id_t, meta::ignore_t>,
                  "Requirements for the field-type of the ID-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(io::detail::lazy_concept_checker([]<typename t = ref_t>(auto) requires(
                    std::ranges::forward_range<t> && (alphabet::nucleotide<std::ranges::range_reference_t<t>> ||
                                                      std::convertible_to<std::ranges::range_reference_t<t>, char>)) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the REF-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(io::detail::lazy_concept_checker([]<typename t = alt_t>(auto) requires(
                    meta::decays_to<id_t, meta::ignore_t> ||
                    (std::ranges::forward_range<t> && std::ranges::forward_range<std::ranges::range_reference_t<t>> &&
                     (alphabet::nucleotide<std::ranges::range_reference_t<std::ranges::range_reference_t<t>>> ||
                      std::convertible_to<std::ranges::range_reference_t<std::ranges::range_reference_t<t>>, char>))) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the ALT-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(meta::arithmetic<std::remove_cvref_t<qual_t>> || meta::decays_to<qual_t, meta::ignore_t>,
                  "Requirements for the field-type of the QUAL-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(io::detail::lazy_concept_checker([]<typename t = filter_t>(auto) requires(
                    meta::decays_to<t, meta::ignore_t> ||
                    (std::ranges::forward_range<t> && (io::detail::char_range<std::ranges::range_reference_t<t>> ||
                                                       std::same_as<std::ranges::range_value_t<t>, int32_t>))) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the FILTER-field not met. See documentation for "
                  "bio::io::var::record.");

    static_assert(
      io::detail::lazy_concept_checker([]<typename t = info_t>(auto) requires(
        meta::decays_to<t, meta::ignore_t> ||
        (std::ranges::forward_range<t> && detail::info_element_writer_concept<std::ranges::range_reference_t<t>>)) {
          return std::true_type{};
      }),
      "Requirements for the field-type of the INFO-field not met. See documentation for "
      "bio::io::var::record.");

    static_assert(io::detail::lazy_concept_checker([]<typename t = genotypes_t>(auto) requires(
                    meta::decays_to<t, meta::ignore_t> ||
                    (ranges::back_insertable<t> &&
                     detail::genotype_writer_concept<std::ranges::range_reference_t<t>>)) { return std::true_type{}; }),
                  "Requirements for the field-type of the GENOTYPES-field not met. See documentation for "
                  "bio::io::var::record.");
    return true;
}

} // namespace bio::io::var::detail
