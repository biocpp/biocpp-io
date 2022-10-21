// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::var_io::record and auxilliary data structures.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/meta/concept/core_language.hpp>
#include <bio/meta/tag/ttag.hpp>
#include <bio/meta/tag/vtag.hpp>
#include <bio/ranges/concept.hpp>
#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/detail/magic_get.hpp>
#include <bio/io/detail/range.hpp>
#include <bio/io/detail/tuple_record.hpp>
#include <bio/io/misc.hpp>
#include <bio/io/var_io/misc.hpp>

//-----------------------------------------------------------------------------
// forwards
//-----------------------------------------------------------------------------

namespace bio::io::var_io
{

class header;

} // namespace bio::io::var_io

//-----------------------------------------------------------------------------
// BCF record core
//-----------------------------------------------------------------------------

namespace bio::io::detail
{

//!\brief The "core" of a BCF record in bit-compatible representation to the on-disk format.
//!\ingroup var_io
struct bcf_record_core
{
    int32_t  chrom         = -1;                                    //!< CHROM as IDX.
    int32_t  pos           = -1;                                    //!< POS.
    int32_t  rlen          = -1;                                    //!< Not used by this implementation.
    float    qual          = bio::io::var_io::missing_value<float>; //!< QUAL.
    uint16_t n_info        = 0;                                     //!< Number of INFOS values.
    uint16_t n_allele      = 0;                                     //!< Number of alleles.
    uint32_t n_sample : 24 = 0;                                     //!< Number of samples.
    uint8_t  n_fmt         = 0;                                     //!< Number of FORMAT values.
};

static_assert(sizeof(bcf_record_core) == 24, "Bit alignment problem in declaration of bcf_record_core.");

} // namespace bio::io::detail

//-----------------------------------------------------------------------------
// The info element
//-----------------------------------------------------------------------------

namespace bio::io::var_io
{

/*!\brief Variant to handle "dynamic typing" in Var I/O INFO fields.
 * \ingroup var_io
 * \details
 *
 * This variant can hold values for the INFO field.
 * Since the type of such fields may only be determined at run-time (depends on values in header), variables
 * of this type can be set to different types at run-time.
 */
template <ownership own = ownership::shallow>
using info_element_value_type =
  std::variant<char,
               int8_t,
               int16_t,
               int32_t,
               float,
               std::conditional_t<own == ownership::shallow, std::string_view, std::string>,
               std::vector<int8_t>,
               std::vector<int16_t>,
               std::vector<int32_t>,
               std::vector<float>,
               std::vector<std::conditional_t<own == ownership::shallow, std::string_view, std::string>>,
               bool>;

} // namespace bio::io::var_io

namespace bio::io::detail
{
//!\brief Auxilliary concept that encompasses bio::io::var_io::info_element_value_type.
//!\ingroup var_io
template <typename t>
concept is_info_element_value_type = meta::
  one_of<t, var_io::info_element_value_type<ownership::shallow>, var_io::info_element_value_type<ownership::deep>>;

} // namespace bio::io::detail

namespace bio::io::var_io
{

/*!\brief The type of elements in an INFO field. [default]
 * \ingroup var_io
 * \tparam own Ownership of the type; see bio::io::ownership.
 */
template <ownership own = ownership::shallow>
struct info_element
{
    //!\brief Type of the ID.
    using string_t = std::conditional_t<own == ownership::shallow, std::string_view, std::string>;

    //!\brief The ID of the element (as a string or string_view).
    string_t                     id;
    //!\brief The value of the element.
    info_element_value_type<own> value;

    //!\brief Defaulted three-way comparisons.
    auto operator<=>(info_element const &) const = default;
};

/*!\brief The type of elements in an INFO field. [full BCF-style]
 * \ingroup var_io
 * \tparam own Ownership of the type; see bio::io::ownership.
 */
template <ownership own = ownership::shallow>
struct info_element_idx
{
    //!\brief The IDX of the element (index of that descriptor in the header).
    int32_t                      idx;
    //!\brief The value of the element.
    info_element_value_type<own> value;

    //!\brief Defaulted three-way comparisons.
    auto operator<=>(info_element_idx const &) const = default;
};

//-----------------------------------------------------------------------------
// The genotype element
//-----------------------------------------------------------------------------

/*!\brief Variant to handle "dynamic typing" in Var I/O GENOTYPE fields.
 * \ingroup var_io
 * \details
 *
 * This type is similar to bio::io::var_io::info_element_value_type except that it encodes a range of the respective
 * value.
 *
 * It does not contain an entry for bio::io::var_io::value_type_id::flag, because flags cannot appear in
 * the genotype field.
 */
template <ownership own = ownership::shallow>
using genotype_element_value_type =
  std::variant<std::vector<char>,
               std::vector<int8_t>,
               std::vector<int16_t>,
               std::vector<int32_t>,
               std::vector<float>,
               std::vector<std::conditional_t<own == ownership::shallow, std::string_view, std::string>>,
               ranges::concatenated_sequences<std::vector<int8_t>>,
               ranges::concatenated_sequences<std::vector<int16_t>>,
               ranges::concatenated_sequences<std::vector<int32_t>>,
               ranges::concatenated_sequences<std::vector<float>>,
               std::vector<std::vector<std::conditional_t<own == ownership::shallow, std::string_view, std::string>>>
               /* no flag here */>;

} // namespace bio::io::var_io

namespace bio::io::detail
{

//!\brief Auxilliary concept that encompasses bio::io::var_io::genotype_element_value_type.
//!\ingroup var_io
template <typename t>
concept is_genotype_element_value_type = meta::one_of<t,
                                                      var_io::genotype_element_value_type<ownership::shallow>,
                                                      var_io::genotype_element_value_type<ownership::deep>>;

} // namespace bio::io::detail

namespace bio::io::var_io
{

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
 * types (one element per sample!), so bio::io::var_io::value_type_id::vector_of_int32 corresponds to
 * std::vector<std::vector<int32_t>>. See bio::io::var_io::genotype_element_value_type for more details.
 *
 * If fields are missing from some samples but not others, the vector will have full size but the respective values
 * will be set to the missing value (see bio::io::var_io::missing_value) or be the empty vector (in case the element
 * type is a vector).
 */
template <ownership own = ownership::shallow>
struct genotype_element
{
    //!\brief Type of the ID.
    using string_t = std::conditional_t<own == ownership::shallow, std::string_view, std::string>;

    //!\brief The ID of the element (as a string or string_view).
    string_t                         id;
    //!\brief The value of the element.
    genotype_element_value_type<own> value;

    //!\brief Defaulted three-way comparisons.
    auto operator<=>(genotype_element const &) const = default;
};

/*!\brief A type representing an element in the GENOTYPES field. [full BCF-style]
 * \ingroup var_io
 *
 * \details
 *
 * The same as bio::io::var_io::genotype_element except that a numeric IDX is used instead of the ID string.
 */
template <ownership own = ownership::shallow>
struct genotype_element_idx
{
    //!\brief The IDX of the element (index of that descriptor in the header).
    int32_t                          idx;
    //!\brief The value of the element.
    genotype_element_value_type<own> value;

    //!\brief Defaulted three-way comparisons.
    auto operator<=>(genotype_element_idx const &) const = default;
};

//-----------------------------------------------------------------------------
// record_private_data
//-----------------------------------------------------------------------------

//!\brief A datastructure that contains private data of variant IO records.
//!\ingroup var_io
struct record_private_data
{
    //!\privatesection
    //!\brief Pointer to the header
    header const * header_ptr = nullptr;

    //!\brief Pointer to record core (if BCF).
    detail::bcf_record_core const * record_core = nullptr;

    //!\brief Raw record type.
    using raw_record_t =
      io::detail::tuple_record<meta::vtag_t<detail::field::chrom,
                                            detail::field::pos,
                                            detail::field::id,
                                            detail::field::ref,
                                            detail::field::alt,
                                            detail::field::qual,
                                            detail::field::filter,
                                            detail::field::info,
                                            detail::field::genotypes,
                                            detail::field::_private>,
                               meta::list_traits::concat<meta::list_traits::repeat<9, std::span<std::byte const>>,
                                                         meta::type_list<var_io::record_private_data>>>;
    //!\brief Pointer to raw record.
    raw_record_t const * raw_record = nullptr;

    //!\brief Defaulted three-way comparison.
    friend bool operator==(record_private_data const &, record_private_data const &) = default;
};

//-----------------------------------------------------------------------------
// The record type
//-----------------------------------------------------------------------------

/*!\brief Record type for variant I/O.
 * \ingroup var_io
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
 * * bio::io::var_io::record_default:
 *   * the record is **deep**
 *   * CHROM, FILTER, INFO element IDs and GENOTYPE element IDs are **strings**
 * * bio::io::var_io::record_default_shallow:
 *   * the record is **shallow**
 *   * CHROM, FILTER, INFO element IDs and GENOTYPE element IDs are **string_views**
 * * bio::io::var_io::record_idx:
 *   * the record is **deep**
 *   * CHROM, FILTER, INFO element IDs and GENOTYPE element IDs are IDX values (int32_t).
 * * bio::io::var_io::record_idx_shallow:
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
 * \snippet test/snippet/var_io/var_io_reader_options.cpp field_types_expert
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
          typename _info_t      = std::vector<info_element<ownership::deep>>,
          typename _genotypes_t = std::vector<genotype_element<ownership::deep>>>
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
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
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
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
     *
     * 1. any std::ranges::forward_range over an alphabet that is convertible to `char` (string).
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
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
     *
     * 1. any std::integral type; int32_t recommended!
     * 2. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
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
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
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
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
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
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
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
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
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
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
     *
     * 1. any bio::meta::arithmetic type
     * 2. bio::meta::ignore_t (**ignored**)
     *
     * Note that '.' will be read as the respective bio::io::var_io::missing_value.
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
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
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over a back-insertable of `char` (**deep**, string)
     *   2. any back-insertable range over std::string_view (**shallow**, string)
     *   3. any back-insertable range over int32_t (IDX values)
     *   4. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
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
     * "Info elements" can contain string identifiers or IDX values.
     * It is possible to use tuples instead of bio::io::var_io::info_element,
     * but this is not recommended.
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over a bio::io::var_io::info_element
     * (deep/shallow depends on template paremeter, string)
     *   2. any back-insertable range over a bio::io::var_io::info_element_idx
     * (deep/shallow depends on template paremeter, IDX)
     *   3. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
     *
     *   1. any std::ranges::forward_range over a bio::io::var_io::info_element
     *   2. any std::ranges::forward_range over a bio::io::var_io::info_element_idx
     *   3. bio::meta::ignore_t (**ignored**)
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
     * "Genotype elements" can contain string identifiers or IDX values.
     * It is possible to use tuples instead of bio::io::var_io::genotype_element,
     * but this is not recommended.
     *
     * ### Type requirements when reading
     *
     * When reading (bio::io::seq_io::reader) the type can be one of the following:
     *
     *   1. any back-insertable range over a bio::io::var_io::genotype_element
     * (deep/shallow depends on template paremeter, string)
     *   2. any back-insertable range over a bio::io::var_io::genotype_element_idx
     * (deep/shallow depends on template paremeter, IDX)
     *   3. bio::meta::ignore_t (**ignored**)
     *
     * See \ref shallow_vs_deep for more details on what "deep" and "shallow" mean here.
     *
     * ### Type requirements when writing
     *
     * When writing (bio::io::seq_io::writer), the type can be one of the following:
     *
     *   1. any std::ranges::forward_range over a bio::io::var_io::genotype_element
     *   2. any std::ranges::forward_range over a bio::io::var_io::genotype_element_idx
     *   3. bio::meta::ignore_t (**ignored**)
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
        detail::clear(chrom);
        detail::clear(pos);
        detail::clear(id);
        detail::clear(ref);
        detail::clear(alt);
        detail::clear(qual);
        detail::clear(filter);
        detail::clear(info);
        detail::clear(genotypes);
        detail::clear(_private);
    }

    //!\brief Defaulted comparison operators.
    friend auto operator<=>(record const & lhs, record const & rhs) = default;
};

//-----------------------------------------------------------------------------
// tie_record()
//-----------------------------------------------------------------------------

/*!\brief Helper function for easily creating a record of references.
 * \ingroup var_io
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
 * \returns A bio::io::var_io::record where all the members are references to this function's parameters.
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
 *!\ingroup var_io
 * \see bio::io::var_io::record
 */
using record_default = record<std::string,                                     // chrom,
                              int32_t,                                         // pos,
                              std::string,                                     // id,
                              std::vector<alphabet::dna5>,                     // ref,
                              std::vector<std::string>,                        // alt,
                              float,                                           // qual,
                              std::vector<std::string>,                        // filter,
                              std::vector<info_element<ownership::deep>>,      // info,
                              std::vector<genotype_element<ownership::deep>>>; // genotypes,

/*!\brief The record type used by bio::io::var_io::reader by default.
 *!\ingroup var_io
 * \see bio::io::var_io::record
 */
using record_default_shallow = record<std::string_view,                                   // chrom,
                                      int32_t,                                            // pos,
                                      std::string_view,                                   // id,
                                      conversion_view_t<alphabet::dna5>,                  // ref,
                                      std::vector<std::string_view>,                      // alt,
                                      float,                                              // qual,
                                      std::vector<std::string_view>,                      // filter,
                                      std::vector<info_element<ownership::shallow>>,      // info,
                                      std::vector<genotype_element<ownership::shallow>>>; // genotypes,

/*!\brief A record type with IDX values (shallow).
 *!\ingroup var_io
 * \see bio::io::var_io::record
 */
using record_idx_shallow = record<int32_t,                                                // chrom,
                                  int32_t,                                                // pos,
                                  std::string_view,                                       // id,
                                  conversion_view_t<alphabet::dna5>,                      // ref,
                                  std::vector<std::string_view>,                          // alt,
                                  float,                                                  // qual,
                                  std::vector<int32_t>,                                   // filter,
                                  std::vector<info_element_idx<ownership::shallow>>,      // info,
                                  std::vector<genotype_element_idx<ownership::shallow>>>; // genotypes,

/*!\brief A record type with IDX values (deep).
 *!\ingroup var_io
 * \see bio::io::var_io::record
 */
using record_idx = record<int32_t,                                             // chrom,
                          int32_t,                                             // pos,
                          std::string,                                         // id,
                          std::vector<alphabet::dna5>,                         // ref,
                          std::vector<std::string>,                            // alt,
                          float,                                               // qual,
                          std::vector<int32_t>,                                // filter,
                          std::vector<info_element_idx<ownership::deep>>,      // info,
                          std::vector<genotype_element_idx<ownership::deep>>>; // genotypes,

} // namespace bio::io::var_io

namespace bio::io::detail // TODO move this to var_io::detail?
{

//-----------------------------------------------------------------------------
// Record concept checks for reading
//-----------------------------------------------------------------------------

/*!\interface bio::io::detail::info_element_reader_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var_io::info_element / bio::io::var_io::info_element_idx.
 */
//!\cond CONCEPT_DEF
template <typename t>
concept info_element_reader_concept = detail::decomposable_into_two<t> &&
  (detail::out_string<detail::first_elem_t<t>> ||
   std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::is_info_element_value_type<detail::second_elem_t<t>>;
//!\endcond

/*!\interface bio::io::detail::genotype_reader_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var_io::genotype_element / bio::io::var_io::genotype_element_idx.
 */
//!\cond CONCEPT_DEF
template <typename t>
concept genotype_reader_concept = detail::decomposable_into_two<t> &&
  (detail::out_string<detail::first_elem_t<t>> ||
   std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::is_genotype_element_value_type<detail::second_elem_t<t>>;
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
  std::type_identity<var_io::record<chrom_t, pos_t, id_t, ref_t, alt_t, qual_t, filter_t, info_t, genotypes_t>>)

{
    // TODO(GCC11): once GCC10 is dropped, remove the "<typename t = seq_t>"
    static_assert(detail::lazy_concept_checker([]<typename t = chrom_t>(auto) requires(
                    ranges::back_insertable_with<t, char> ||
                    meta::one_of<std::remove_reference_t<t>, std::string_view, int32_t, ignore_t, ignore_t const>) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the CHROM-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    static_assert(detail::lazy_concept_checker([]<typename t = pos_t>(auto) requires(
                    std::integral<std::remove_reference_t<t>> ||
                    meta::one_of<std::remove_reference_t<t>, ignore_t, ignore_t const>) { return std::true_type{}; }),
                  "Requirements for the field-type of the POS-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    static_assert(detail::lazy_concept_checker([]<typename t = id_t>(auto) requires(
                    ranges::back_insertable_with<t, char> ||
                    meta::one_of<std::remove_reference_t<t>, std::string_view, ignore_t, ignore_t const>) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the ID-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    static_assert(detail::lazy_concept_checker([]<typename t = ref_t>(auto) requires(
                    (ranges::back_insertable<t> && alphabet::alphabet<std::ranges::range_reference_t<t>>) ||
                    meta::one_of<std::remove_reference_t<t>, std::string_view, ignore_t, ignore_t const> ||
                    detail::transform_view_on_string_view<t>) { return std::true_type{}; }),
                  "Requirements for the field-type of the REF-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    static_assert(
      detail::lazy_concept_checker([]<typename t = alt_t>(auto) requires(
        meta::one_of<std::remove_reference_t<t>, ignore_t, ignore_t const> ||
        (ranges::back_insertable<t> &&
         ((ranges::back_insertable<std::ranges::range_reference_t<t>> &&
           alphabet::alphabet<std::ranges::range_reference_t<std::ranges::range_reference_t<t>>>) ||
          meta::decays_to<std::ranges::range_reference_t<t>, std::string_view> ||
          detail::transform_view_on_string_view<std::ranges::range_reference_t<t>>))) { return std::true_type{}; }),
      "Requirements for the field-type of the ALT-field not met. See documentation for "
      "bio::io::var_io::reader_options.");

    static_assert(detail::lazy_concept_checker([]<typename t = qual_t>(auto) requires(
                    meta::arithmetic<std::remove_reference_t<t>> ||
                    meta::one_of<std::remove_reference_t<t>, ignore_t, ignore_t const>) { return std::true_type{}; }),
                  "Requirements for the field-type of the QUAL-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    static_assert(
      detail::lazy_concept_checker([]<typename t = filter_t>(auto) requires(
        meta::one_of<std::remove_reference_t<t>, ignore_t, ignore_t const> ||
        (ranges::back_insertable<t> &&
         (ranges::back_insertable_with<std::ranges::range_reference_t<t>, char> ||
          meta::one_of<std::remove_reference_t<std::ranges::range_reference_t<t>>, std::string_view, int32_t>))) {
          return std::true_type{};
      }),
      "Requirements for the field-type of the FILTER-field not met. See documentation for "
      "bio::io::var_io::reader_options.");

    static_assert(detail::lazy_concept_checker([]<typename t = info_t>(auto) requires(
                    meta::one_of<std::remove_reference_t<t>, ignore_t, ignore_t const> ||
                    (ranges::back_insertable<t> &&
                     detail::info_element_reader_concept<std::remove_reference_t<std::ranges::range_reference_t<t>>>)) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the INFO-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    static_assert(detail::lazy_concept_checker([]<typename t = genotypes_t>(auto) requires(
                    meta::one_of<std::remove_reference_t<t>, ignore_t, ignore_t const> ||
                    (ranges::back_insertable<t> &&
                     detail::genotype_reader_concept<std::remove_reference_t<std::ranges::range_reference_t<t>>>)) {
                      return std::true_type{};
                  }),
                  "Requirements for the field-type of the GENOTYPES-field not met. See documentation for "
                  "bio::io::var_io::reader_options.");

    return true;
}

//-----------------------------------------------------------------------------
// Record concept checks for writing
//-----------------------------------------------------------------------------

//!\brief A char, signed_integral, floating point number or CString.
template <typename t>
concept var_io_legal_type_aux =
  std::same_as<t, char> || std::signed_integral<t> || std::floating_point<t> || std::same_as < std::decay_t<t>,
char const * > ;

/*!\interface bio::io::detail::var_io_legal_type <>
 * \tparam t The type to check.
 * \brief A type that is similar to one of the alternatives of bio::io::var_io::info_element_value_type
 */
//!\cond CONCEPT_DEF
template <typename t>
concept var_io_legal_type = var_io_legal_type_aux<std::remove_cvref_t<t>> || std::same_as<t const &, bool const &> ||
  (std::ranges::forward_range<t> && (var_io_legal_type_aux<std::remove_cvref_t<std::ranges::range_reference_t<t>>> ||
                                     (std::ranges::forward_range<std::ranges::range_reference_t<t>> &&
                                      std::same_as<char const &, std::ranges::range_reference_t<t> const &>)));
//!\endcond

/*!\interface bio::io::detail::var_io_legal_vector_type <>
 * \tparam t The type to check.
 * \brief A type that is similar to one of the alternatives of bio::io::var_io::info_element_value_type
 */
//!\cond CONCEPT_DEF
template <typename t>
concept var_io_legal_vector_type =
  std::ranges::forward_range<t> && var_io_legal_type<std::ranges::range_reference_t<t>> &&
  !std::same_as<bool const &, std::ranges::range_reference_t<t>>;
//!\endcond

/*!\interface bio::io::detail::var_io_legal_or_dynamic <>
 * \tparam t The type to check.
 * \brief A type that is similar to one of the alternatives of bio::io::var_io::info_element_value_type
 */
//!\cond CONCEPT_DEF
template <typename t>
concept var_io_legal_or_dynamic = var_io_legal_type<t> || is_info_element_value_type<t>;
//!\endcond

/*!\interface bio::io::detail::var_io_vector_legal_or_dynamic <>
 * \tparam t The type to check.
 * \brief A type that is similar to one of the alternatives of bio::io::var_io::info_element_value_type
 */
//!\cond CONCEPT_DEF
template <typename t>
concept var_io_vector_legal_or_dynamic = var_io_legal_vector_type<t> || is_genotype_element_value_type<t>;
//!\endcond

/*!\interface bio::io::detail::info_element_writer_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var_io::info_element / bio::io::var_io::info_element_idx.
 */
//!\cond CONCEPT_DEF
template <typename t>
concept info_element_writer_concept = detail::decomposable_into_two<t> &&
  (detail::char_range_or_cstring<detail::first_elem_t<t>> ||
   std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::var_io_legal_or_dynamic<detail::second_elem_t<t>>;
//!\endcond

/*!\interface bio::io::detail::genotype_writer_concept <>
 * \tparam t The type to check.
 * \brief Types "similar" to bio::io::var_io::genotype_element / bio::io::var_io::genotype_element_idx.
 */
//!\cond CONCEPT_DEF
template <typename t>
concept genotype_writer_concept = detail::decomposable_into_two<t> &&
  (detail::char_range_or_cstring<detail::first_elem_t<t>> ||
   std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::var_io_vector_legal_or_dynamic<detail::second_elem_t<t>>;
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
  std::type_identity<var_io::record<chrom_t, pos_t, id_t, ref_t, alt_t, qual_t, filter_t, info_t, genotypes_t>>)

{
    // TODO
    return true;
}

} // namespace bio::io::detail
