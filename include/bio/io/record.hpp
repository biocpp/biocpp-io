// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::record template and the bio::io::field enum.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/core/detail/template_inspection.hpp>

#include <bio/meta/tag/vtag.hpp>
#include <bio/meta/type_list/traits.hpp>
#include <bio/meta/type_list/type_list.hpp>

#include <bio/io/detail/concept.hpp>
#include <bio/io/misc.hpp>

namespace bio::io
{

// ----------------------------------------------------------------------------
// enum field
// ----------------------------------------------------------------------------

/*!\brief An enumerator for the fields used in file formats.
 * \ingroup io
 *
 * \details
 *
 * Some of the fields are shared between formats.
 *
 * TODO improve documentation
 */
enum class field : uint64_t
{
    // Fields used in multiple contexts
    seq,      //!< The "sequence", usually a range of nucleotides or amino acids.
    id,       //!< The identifier, usually a string.
    qual,     //!< The qualities, usually in phred-score notation.
    seq_qual, //!< Sequence and qualities combined in one range.
    offset,   //!< Sequence (SEQ) relative start position (0-based), unsigned value.
    ref_id,
    ref_seq,  //!< The (reference) "sequence" information, usually a range of nucleotides or amino acids.
    pos,      //!< Sequence (REF_SEQ) relative start position (0-based), unsigned value.
    _private, //!< Refers to arbitrary internal datastructures (**never modify this field**).

    // Fields unique to alignment io
    qname = id,
    flag  = _private + 1, //!< The alignment flag (bit information), `uint16_t` value.
    /*ref_id*/            //!< The identifier of the (reference) sequence that SEQ was aligned to.
    /*pos*/
    mapq,  //!< The mapping quality of the SEQ alignment, usually a ohred-scaled score.
    cigar, //!< The cigar vector (std::vector<seqan3::cigar>) representing the alignment in SAM/BAM format.
    next_ref_id,
    next_pos,
    tlen,
    /*seq*/
    /*qual*/
    optionals, //!< The optional fields in the SAM format, stored in a dictionary.
    /*_private*/

    // Fields unique to variant io
    chrom = ref_id, //!< CHROM field in Var I/O.
    /*pos*/
    /*id*/
    ref   = optionals + 1, //!< REF field in Var I/O.
    alt,                   //!< ALT field in Var I/O.
    /*qual*/
    filter,    //!< FILTER field in Var I/O.
    info,      //!< INFO field in Var I/O.
    genotypes, //!< GENOTYPES in Var I/O.
    /*_private*/

    // User defined field aliases
    user_defined = uint64_t{1} << 32, //!< Identifier for user defined file formats and specialisations.
};

} // namespace bio::io

namespace bio::io::detail
{

//!\brief Checks whether a type is a bio::meta::vtag_t over bio::io::field.
template <typename t>
inline constexpr bool is_fields_tag = false;

//!\brief Checks whether a type is a bio::meta::vtag_t over bio::io::field.
template <field... vs>
inline constexpr bool is_fields_tag<meta::vtag_t<vs...>> = true;

} // namespace bio::io::detail

namespace bio::io
{

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

/*!\brief The class template that file records are based on; behaves like an std::tuple.
 * \implements seqan3::tuple_like
 * \ingroup io
 * \tparam field_ids   A meta::vtag_t type with bio::io::field IDs corresponding to field_types.
 * \tparam field_types The types of the fields in this record as a meta::type_list.
 * \details
 *
 * This class template behaves like a std::tuple, with the exception that it provides an additional
 * get-interface that takes a bio::io::field identifier. The traditional get interfaces (via index and
 * via type) are also supported, but discouraged, because accessing via bio::io::field is unambiguous and
 * better readable.
 *
 * In addition to the get()-interfaces, member accessors are provided with the same name as the fields.
 *
 * See bio::io::seq_io::reader for how this data structure is used in practice.
 *
 * See #make_record() and #tie_record() for easy ways to create stand-alone record variables.
 *
 * See the \ref record_faq for more details.
 */
template <typename field_ids_, typename field_types_>
struct record : seqan3::detail::transfer_template_args_onto_t<field_types_, std::tuple>
{
public:
    using field_types = field_types_; //!< The field types as a type_list.
    using field_ids   = field_ids_;   //!< The field ids corresponding to the field_types.

    //!\brief A specialisation of std::tuple.s
    using base_type = seqan3::detail::transfer_template_args_onto_t<field_types, std::tuple>;

private:
    //!\brief Auxiliary functions for clear().
    template <typename t>
        //!\cond REQ
        requires(requires(t & v) { v.clear(); })
    //!\endcond
    static constexpr void clear_element(t & v) noexcept(noexcept(v.clear())) { v.clear(); }

    //!\overload
    template <typename t>
    static constexpr void clear_element(t & v) noexcept(noexcept(std::declval<t &>() = t{}))
    {
        v = {};
    }

    //!\brief A lambda function that expands a pack and calls `clear_element` on every argument in the pack.
    static constexpr auto expander = [](auto &... args) { (clear_element(args), ...); };

    //!\brief Get the base_type.
    base_type & to_base() { return static_cast<base_type &>(*this); }

    //!\brief Get the base_type.
    base_type const & to_base() const { return static_cast<base_type const &>(*this); }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    record()                           = default; //!< Defaulted.
    record(record const &)             = default; //!< Defaulted.
    record(record &&)                  = default; //!< Defaulted.
    ~record()                          = default; //!< Defaulted.
    record & operator=(record const &) = default; //!< Defaulted.
    record & operator=(record &&)      = default; //!< Defaulted.

    //!\brief Inherit tuple's constructors.
    using base_type::base_type;
    //!\}

    static_assert(meta::list_traits::size<field_types> == field_ids::size,
                  "You must give as many IDs as types to bio::io::record.");

    //!\brief Clears containers that provide `.clear()` and (re-)initialises all other elements with `= {}`.
    void clear() noexcept(noexcept(std::apply(expander, std::declval<record &>()))) { std::apply(expander, *this); }

    /*!\name Get accessors
     * \{
     */
    //!\brief Get a specific field by it's position.
    template <size_t i>
    decltype(auto) get() noexcept(noexcept(std::get<i>(to_base())))
    {
        return std::get<i>(to_base());
    }

    //!\brief Get a specific field by it's position.
    template <size_t i>
    decltype(auto) get() const noexcept(noexcept(std::get<i>(std::as_const(to_base()))))
    {
        return std::get<i>(to_base());
    }

    //!\brief Get a specific field by it's field id.
    template <field f>
        //!\cond REQ
        requires(field_ids::contains(f))
    //!\endcond
    decltype(auto) get() noexcept(noexcept(std::get<field_ids::index_of(f)>(to_base())))
    {
        return std::get<field_ids::index_of(f)>(to_base());
    }

    //!\brief Get a specific field by it's field id.
    template <field f>
        //!\cond REQ
        requires(field_ids::contains(f))
    //!\endcond
    decltype(auto) get() const noexcept(noexcept(std::get<field_ids::index_of(f)>(to_base())))
    {
        return std::get<field_ids::index_of(f)>(to_base());
    }
    //!\}

//!\brief A macro that defines all getter functions for fields contained in bio::io::record.
#define BIOCPP_IO_RECORD_MEMBER(F)                                                                                     \
    /*!\brief Return the bio::io::field F if available.*/                                                              \
    decltype(auto) F() noexcept(noexcept(get<field::F>()))                                                             \
    {                                                                                                                  \
        return get<field::F>();                                                                                        \
    }                                                                                                                  \
    /*!\brief Return the bio::io::field F if available. [const-qualified version] */                                   \
    decltype(auto) F() const noexcept(noexcept(get<field::F>()))                                                       \
    {                                                                                                                  \
        return get<field::F>();                                                                                        \
    }

    /*!\name Member accessors
     * \brief This is the same as calling #get<field::X>(); functions are only defined if record has that element.
     * \{
     */
    BIOCPP_IO_RECORD_MEMBER(seq)
    BIOCPP_IO_RECORD_MEMBER(id)
    BIOCPP_IO_RECORD_MEMBER(qual)
    BIOCPP_IO_RECORD_MEMBER(seq_qual)
    BIOCPP_IO_RECORD_MEMBER(offset)
    BIOCPP_IO_RECORD_MEMBER(ref_id)
    BIOCPP_IO_RECORD_MEMBER(ref_seq)
    BIOCPP_IO_RECORD_MEMBER(pos)
    BIOCPP_IO_RECORD_MEMBER(qname)
    BIOCPP_IO_RECORD_MEMBER(flag)
    BIOCPP_IO_RECORD_MEMBER(mapq)
    BIOCPP_IO_RECORD_MEMBER(cigar)
    BIOCPP_IO_RECORD_MEMBER(next_ref_id)
    BIOCPP_IO_RECORD_MEMBER(next_pos)
    BIOCPP_IO_RECORD_MEMBER(tlen)
    BIOCPP_IO_RECORD_MEMBER(optionals)
    BIOCPP_IO_RECORD_MEMBER(chrom)
    BIOCPP_IO_RECORD_MEMBER(ref)
    BIOCPP_IO_RECORD_MEMBER(alt)
    BIOCPP_IO_RECORD_MEMBER(filter)
    BIOCPP_IO_RECORD_MEMBER(info)
    BIOCPP_IO_RECORD_MEMBER(genotypes)
    //!\}
#undef BIOCPP_IO_RECORD_MEMBER
};

} // namespace bio::io

//-------------------------------------------------------------------------------
// tuple traits
//-------------------------------------------------------------------------------

namespace std
{

/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::unary_type_trait
 * \relates bio::io::record
 * \see std::tuple_size_v
 */
template <typename field_ids, typename field_types>
struct tuple_size<bio::io::record<field_ids, field_types>>
{
    //!\brief The value member. Delegates to same value on base_type.
    static constexpr size_t value = tuple_size_v<typename bio::io::record<field_ids, field_types>::base_type>;
};

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::transformation_trait
 * \relates bio::io::record
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template <size_t elem_no, typename field_ids, typename field_types>
struct tuple_element<elem_no, bio::io::record<field_ids, field_types>>
{
    //!\brief The member type. Delegates to same type on base_type.
    using type = std::tuple_element_t<elem_no, typename bio::io::record<field_ids, field_types>::base_type>;
};

} // namespace std

namespace bio::io
{

//-------------------------------------------------------------------------------
// record_element
//-------------------------------------------------------------------------------

/*!\brief Like std::tuple_element but with bio::io::field on bio::io::record. [declaration]
 * \implements seqan3::transformation_trait
 * \relates bio::io::record
 */
template <field f, typename t>
struct record_element;

//!\brief Like std::tuple_element but with bio::io::field on bio::io::record. [implementation]
template <field f, typename field_ids, typename field_types>
    //!\cond REQ
    requires(field_ids::contains(f))
//!\endcond
struct record_element<f, record<field_ids, field_types>> :
  public std::tuple_element<field_ids::index_of(f), record<field_ids, field_types>>
{};

//!\brief Like std::tuple_element but with bio::io::field on bio::io::record. [type trait shortcut]
template <field f, typename t>
    //!\cond REQ
    requires(requires { typename record_element<f, t>::type; })
//!\endcond
using record_element_t = typename record_element<f, t>::type;

//-------------------------------------------------------------------------------
// bio::io::get
//-------------------------------------------------------------------------------

/*!\name Free function get() interface for bio::io::record based on bio::io::field.
 * \brief This is the tuple interface via bio::io::field, e.g. `get<field::seq>(tuple)`.
 * \relates bio::io::record
 * \{
 */
//!\brief Free function get() for bio::io::record based on bio::io::field.
template <field f, typename field_ids, typename field_types>
auto & get(record<field_ids, field_types> & r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(r);
}

//!\overload
template <field f, typename field_ids, typename field_types>
auto const & get(record<field_ids, field_types> const & r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(r);
}

//!\overload
template <field f, typename field_ids, typename field_types>
auto && get(record<field_ids, field_types> && r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(std::move(r));
}

//!\overload
template <field f, typename field_ids, typename field_types>
auto const && get(record<field_ids, field_types> const && r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(std::move(r));
}
//!\}

// Implementation note: for some reason, the following is already "related" do bio::io::record
/*!\name Convenience functions for creating bio::io::record.
 * \{
 */
//-------------------------------------------------------------------------------
// make_record
//-------------------------------------------------------------------------------

/*!\brief Create a deep bio::io::record from the arguments (like std::make_tuple for std::tuple).
 * \param[in] tag    A tag that specifies the identifiers of the subsequent arguments.
 * \param[in] fields The arguments to put into the record.
 * \returns A bio::io::record with copies of the field arguments.
 * \details
 *
 * The record will contain copies of the arguments.
 *
 * For more information, see \ref record_type and \ref record_make_tie
 *
 * ### Example
 *
 * \snippet test/snippet/record.cpp make_and_tie_record
 */
template <auto... field_ids, typename... field_type_ts>
constexpr auto make_record(meta::vtag_t<field_ids...> BIOCPP_IO_DOXYGEN_ONLY(tag), field_type_ts &&... fields)
  -> record<meta::vtag_t<field_ids...>, bio::meta::type_list<std::decay_t<field_type_ts>...>>
{
    return {std::forward<field_type_ts>(fields)...};
}

//-------------------------------------------------------------------------------
// tie_record
//-------------------------------------------------------------------------------

/*!\brief Create a shallow bio::io::record from the arguments (like std::tie for std::tuple).
 * \param[in] tag    A tag that specifies the identifiers of the subsequent arguments.
 * \param[in] fields The arguments to represent in the record.
 * \returns A bio::io::record with references to the field arguments.
 * \details
 *
 * The record will contain references to the arguments.
 *
 * For more information, see \ref record_type and \ref record_make_tie
 *
 * ### Example
 *
 * \snippet test/snippet/record.cpp make_and_tie_record
 */
template <auto... field_ids, typename... field_type_ts>
constexpr auto tie_record(meta::vtag_t<field_ids...> BIOCPP_IO_DOXYGEN_ONLY(tag), field_type_ts &... fields)
  -> record<meta::vtag_t<field_ids...>, bio::meta::type_list<field_type_ts &...>>
{
    return {fields...};
}

//!\}

} // namespace bio::io
