// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::format_input_handler and bio::format_input_handler_base.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <span>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/views/char_strictly_to.hpp>

#include <bio/detail/charconv.hpp>
#include <bio/detail/concept.hpp>
#include <bio/detail/range.hpp>
#include <bio/exception.hpp>
#include <bio/record.hpp>

namespace bio
{

/*!\brief A template that provides the reading and parsing functionality of the specified format.
 * \ingroup format
 * \details
 *
 * Formats need to specialise this template and provide a constructor that accepts an std::istream as well as a
 * public member function with the following signature:
 * ```cpp
 * void parse_next_record_into(auto & parsed_record);
 * ```
 * It must accept any bio::record and "fill" that record with the content found in that file.
 *
 * This template may be specialised with a user-provided type, however the process is non-trivial. Documentation
 * can be found here (TODO).
 */
template <typename format_t>
class format_input_handler;

/*!\brief A CRTP base class that helps implement bio::format_input_handler.
 * \ingroup format
 * \details
 *
 * The details of this base class are currently not exposed and shall thus not be relied upon by user code.
 */
template <typename derived_t>
class format_input_handler_base
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief Befriend the derived type so it can instantiate.
    friend derived_t;

    //!\brief Downcast self to derived type.
    derived_t * to_derived() { return static_cast<derived_t *>(this); }

    //!\brief Downcast self to derived type. [const-qualified version]
    derived_t const * to_derived() const { return static_cast<derived_t const *>(this); }
    //!\}

    /* members */
    //!\brief A pointer to the input stream.
    std::istream * stream = nullptr;

    /* stuff for turning raw record into parsed record */

    /*!\name Parsing individual fields - defaults (step 3)
     * \brief These are sane defaults for text based formats (not taken for binary formats).
     * \{
     */
    //!\brief Not parsing at all / *raw IO*; binary formats' raw record is bytes already.
    static inline void parse_field_aux(std::span<std::byte const> const in, std::span<std::byte const> & parsed_field)
    {
        parsed_field = in;
    }

    //!\brief Not parsing at all / *raw IO*.
    static inline void parse_field_aux(std::string_view const in, std::span<std::byte const> & parsed_field)
    {
        parsed_field = std::span<std::byte const>{reinterpret_cast<std::byte const *>(in.data()), in.size()};
    }

    //!\brief Parsing into string views. override this.
    static inline void parse_field_aux(std::string_view const in, std::string_view & parsed_field)
    {
        parsed_field = in;
    }

    //!\brief Parsing into transformed string views.
    template <typename fun_t>
    static void parse_field_aux(std::string_view const                                 in,
                                std::ranges::transform_view<std::string_view, fun_t> & parsed_field)
    {
        parsed_field = {in, fun_t{}};
    }

    //!\brief Parsing into transformed string views.
    template <typename fun1_t, typename fun2_t>
    static void parse_field_aux(
      std::string_view const                                                                       in,
      std::ranges::transform_view<std::ranges::transform_view<std::string_view, fun1_t>, fun2_t> & parsed_field)
    {
        parsed_field = {
          {in, fun1_t{}},
          fun2_t{          }
        };
    }

    //!\brief Parse into string-like types.
    template <detail::back_insertable parsed_field_t>
        requires detail::char_range<parsed_field_t>
    static void parse_field_aux(std::string_view const in, parsed_field_t & parsed_field)
    {
        detail::sized_range_copy(in, parsed_field);
    }

    //!\brief Parse into containers of alphabets.
    template <detail::back_insertable parsed_field_t>
        requires detail::deliberate_alphabet<std::ranges::range_reference_t<parsed_field_t>>
    static void parse_field_aux(std::string_view const in, parsed_field_t & parsed_field)
    {
        using target_alph_type = std::ranges::range_value_t<parsed_field_t>;
        //         detail::sized_range_copy(in | seqan3::views::char_strictly_to<target_alph_type>,
        detail::sized_range_copy(in | seqan3::views::char_strictly_to<target_alph_type>, parsed_field);
    }

    //!\brief Parse into a numerical type.
    static void parse_field_aux(std::string_view const in, seqan3::arithmetic auto & parsed_field)
    {
        detail::string_to_number(in, parsed_field);
    }
    //!\}

    /*!\name Parsing individual fields (step 2)
     * \brief The second step in parsing consists of field-specific parsing. This is typically provided by derived_t.
     * \{
     */
    //!\brief Various target types have sane default implementations.
    template <field field_id>
    void parse_field(vtag_t<field_id> const & /**/, auto & parsed_field)
    {
        if constexpr (requires { derived_t::parse_field_aux(get<field_id>(to_derived()->raw_record), parsed_field); })
        {
            to_derived()->parse_field_aux(get<field_id>(to_derived()->raw_record), parsed_field);
        }
        else
        {
            static_assert(seqan3::arithmetic<decltype(parsed_field)> /*always false*/,
                          "Format X does not know how to parse field Y into type Z. Provide different traits or a "
                          "custom format handler.");
        }
    }
    //!\}

    /*!\name Parsing record (step 1)
     * \brief The first step in parsing decomposes the record.
     * \{
     */
    //!\brief Only act on those fields that are present in the record and also provided by the format.
    template <field field_id, typename parsed_record_t>
    void parse_record_impl(vtag_t<field_id> const & /**/, parsed_record_t & parsed_record)
    {
        if constexpr (parsed_record_t::field_ids::contains(field_id))
        {
            auto & parsed_field = get<field_id>(parsed_record);
            to_derived()->parse_field(vtag<field_id>, parsed_field);
        }
        // fields that are not in format or not in target record are simply ignored
    }

    //!\brief Splits the record into individual fields.
    template <field... field_ids, typename parsed_record_t>
    void parse_record(vtag_t<field_ids...> const & /**/, parsed_record_t & parsed_record)
    {
        (to_derived()->parse_record_impl(vtag<field_ids>, parsed_record), ...);
    }
    //!\}

    /*!\name Constructors, destructor and assignment.
     * \brief These are all private to prevent wrong instantiation.
     * \{
     */
    format_input_handler_base()                                              = default; //!< Defaulted.
    format_input_handler_base(format_input_handler_base const &)             = delete;  //!< Deleted.
    format_input_handler_base(format_input_handler_base &&)                  = default; //!< Defaulted.
    ~format_input_handler_base()                                             = default; //!< Defaulted.
    format_input_handler_base & operator=(format_input_handler_base const &) = delete;  //!< Deleted.
    format_input_handler_base & operator=(format_input_handler_base &&)      = default; //!< Defaulted.

    //!\brief Construct from a std::istream.
    format_input_handler_base(std::istream & str) : stream{&str}
    {
        if (std::istreambuf_iterator<char>{*stream} == std::istreambuf_iterator<char>{})
            throw bio::file_open_error{"The file was empty."};
    }
    //!\}

public:
    /*!\brief Parse input into the record argument.
     * \param[out] parsed_record The out-parameter.
     * \details
     *
     * ### Attention
     *
     * Most users should not use this interface and should instead use formatted readers, e.g.
     * bio::map_io::reader, bio::seq_io::reader or bio::var_io::reader.
     */
    void parse_next_record_into(auto & parsed_record)
    {
        // create new raw record
        to_derived()->read_raw_record();

        // create new parsed record
        parsed_record.clear();
        to_derived()->parse_record(typename derived_t::format_fields{}, parsed_record);
    }
};

} // namespace bio
