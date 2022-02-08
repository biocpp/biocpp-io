// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::record template and the bio::field enum.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

#include <bio/detail/concept.hpp>
#include <bio/misc.hpp>

namespace bio
{

// ----------------------------------------------------------------------------
// enum field
// ----------------------------------------------------------------------------

/*!\brief An enumerator for the fields used in file formats.
 * \ingroup bio
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

} // namespace bio

namespace bio::detail
{

//!\brief Checks whether a type is a bio::vtag_t over bio::field.
template <typename t>
inline constexpr bool is_fields_tag = false;

//!\brief Checks whether a type is a bio::vtag_t over bio::field.
template <field... vs>
inline constexpr bool is_fields_tag<vtag_t<vs...>> = true;

} // namespace bio::detail

namespace bio
{

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

/*!\brief The class template that file records are based on; behaves like an std::tuple.
 * \implements seqan3::tuple_like
 * \ingroup bio
 * \tparam field_types The types of the fields in this record as a seqan3::type_list.
 * \tparam field_ids   A vtag_t type with bio::field IDs corresponding to field_types.
 * \details
 *
 * This class template behaves just like an std::tuple, with the exception that it provides an additional
 * get-interface that takes a bio::field identifier. The traditional get interfaces (via index and
 * via type) are also supported, but discouraged, because accessing via bio::field is unambiguous and
 * better readable.
 *
 * ### Example
 *
 * For input files this template is specialised automatically and provided by the file via its `record_type` member.
 * For output files you my define it locally and pass instances of this to the output file's `push_back()`.
 *
 * This is how it works:
 *
 * \todo include test/snippet/io/record_2.cpp
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

    static_assert(seqan3::list_traits::size<field_types> == field_ids::size,
                  "You must give as many IDs as types to bio::record.");

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

//!\brief A macro that defines all getter functions for fields contained in bio::record.
#define BIO_RECORD_MEMBER(F)                                                                                           \
    /*!\brief Return the bio::field F if available.*/                                                                  \
    decltype(auto) F() noexcept(noexcept(get<field::F>()))                                                             \
    {                                                                                                                  \
        return get<field::F>();                                                                                        \
    }                                                                                                                  \
    /*!\brief Return the bio::field F if available. [const-qualified version] */                                       \
    decltype(auto) F() const noexcept(noexcept(get<field::F>()))                                                       \
    {                                                                                                                  \
        return get<field::F>();                                                                                        \
    }

    /*!\name Member accessors
     * \brief This is the same as calling #get<field::X>(); functions are only defined if record has that element.
     * \{
     */
    BIO_RECORD_MEMBER(seq)
    BIO_RECORD_MEMBER(id)
    BIO_RECORD_MEMBER(qual)
    BIO_RECORD_MEMBER(seq_qual)
    BIO_RECORD_MEMBER(offset)
    BIO_RECORD_MEMBER(ref_id)
    BIO_RECORD_MEMBER(ref_seq)
    BIO_RECORD_MEMBER(pos)
    BIO_RECORD_MEMBER(qname)
    BIO_RECORD_MEMBER(flag)
    BIO_RECORD_MEMBER(mapq)
    BIO_RECORD_MEMBER(cigar)
    BIO_RECORD_MEMBER(next_ref_id)
    BIO_RECORD_MEMBER(next_pos)
    BIO_RECORD_MEMBER(tlen)
    BIO_RECORD_MEMBER(optionals)
    BIO_RECORD_MEMBER(chrom)
    BIO_RECORD_MEMBER(ref)
    BIO_RECORD_MEMBER(alt)
    BIO_RECORD_MEMBER(filter)
    BIO_RECORD_MEMBER(info)
    BIO_RECORD_MEMBER(genotypes)
    //!\}
#undef BIO_RECORD_MEMBER
};

} // namespace bio

//-------------------------------------------------------------------------------
// tuple traits
//-------------------------------------------------------------------------------

namespace std
{

/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::unary_type_trait
 * \relates bio::record
 * \see std::tuple_size_v
 */
template <typename field_ids, typename field_types>
struct tuple_size<bio::record<field_ids, field_types>>
{
    //!\brief The value member. Delegates to same value on base_type.
    static constexpr size_t value = tuple_size_v<typename bio::record<field_ids, field_types>::base_type>;
};

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::transformation_trait
 * \relates bio::record
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template <size_t elem_no, typename field_ids, typename field_types>
struct tuple_element<elem_no, bio::record<field_ids, field_types>>
{
    //!\brief The member type. Delegates to same type on base_type.
    using type = std::tuple_element_t<elem_no, typename bio::record<field_ids, field_types>::base_type>;
};

} // namespace std

namespace bio
{

//-------------------------------------------------------------------------------
// record_element
//-------------------------------------------------------------------------------

/*!\brief Like std::tuple_element but with bio::field on bio::record. [declaration]
 * \implements seqan3::transformation_trait
 * \relates bio::record
 */
template <field f, typename t>
struct record_element;

//!\brief Like std::tuple_element but with bio::field on bio::record. [implementation]
template <field f, typename field_ids, typename field_types>
    requires(field_ids::contains(f))
struct record_element<f, record<field_ids, field_types>> :
  public std::tuple_element<field_ids::index_of(f), record<field_ids, field_types>>
{};

//!\brief Like std::tuple_element but with bio::field on bio::record. [type trait shortcut]
template <field f, typename t>
    requires(requires { typename record_element<f, t>::type; })
using record_element_t = typename record_element<f, t>::type;

//-------------------------------------------------------------------------------
// bio::get
//-------------------------------------------------------------------------------

/*!\name Free function get() interface for bio::record based on bio::field.
 * \brief This is the tuple interface via bio::field, e.g. `get<field::seq>(tuple)`.
 * \relates bio::record
 * \{
 */
//!\brief Free function get() for bio::record based on bio::field.
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

// Implementation note: for some reason, the following is already "related" do bio::record
/*!\name Convenience functions for creating bio::record.
 * \{
 */
//-------------------------------------------------------------------------------
// make_record
//-------------------------------------------------------------------------------

/*!\brief Create a bio::record and deduce type from arguments (like std::make_tuple for std::tuple).
 * \details
 *
 * ### Example
 *
 * TODO
 */
template <auto... field_ids, typename... field_type_ts>
constexpr auto make_record(vtag_t<field_ids...>, field_type_ts &&... fields)
  -> record<vtag_t<field_ids...>, seqan3::type_list<std::decay_t<field_type_ts>...>>
{
    return {std::forward<field_type_ts>(fields)...};
}

//-------------------------------------------------------------------------------
// tie_record
//-------------------------------------------------------------------------------

/*!\brief Create a bio::record of references (like std::tie for std::tuple).
 * \details
 *
 * ### Example
 *
 * TODO
 */
template <auto... field_ids, typename... field_type_ts>
constexpr auto tie_record(vtag_t<field_ids...>, field_type_ts &... fields)
  -> record<vtag_t<field_ids...>, seqan3::type_list<field_type_ts &...>>
{
    return {fields...};
}

//!\}

} // namespace bio
