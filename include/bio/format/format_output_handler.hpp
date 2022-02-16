// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the seqan3::format_output_handler_base.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <bio/detail/concept.hpp>
#include <bio/record.hpp>
#include <bio/stream/detail/fast_streambuf_iterator.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/views/to_char.hpp>

namespace bio
{

/*!\brief A template that provides the writing functionality of the specified format.
 * \ingroup format
 * \details
 *
 * Formats need to specialise this template and provide a constructor that accepts an std::ostream as well as a
 * public member function with the following signature:
 *
 * ```cpp
 * void write_record(bio::record<field_types, field_ids> && record)
 * ```
 * It must accept any bio::record and write that record's fields to the file.
 *
 * This template may be specialised with a user-provided type, however the process is non-trivial. Documentation
 * can be found here (TODO).
 */
template <typename format_t>
class format_output_handler;

/*!\brief A CRTP base class that helps implement bio::format_output_handler specialisation.
 * \ingroup format
 * \details
 *
 * The details of this base class are currently not exposed and shall thus not be relied upon by user code.
 */
template <typename derived_t>
class format_output_handler_base
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

    /*!\name State
     * \{
     */
    //!\brief A pointer to the output stream.
    std::ostream *                         stream = nullptr;
    //!\brief A more efficient stream iterator.
    detail::fast_ostreambuf_iterator<char> it;
    //!\}

    /*!\name Writing individual fields - defaults (step 3)
     * \{
     */
    //!\brief Write chracter ranges.
    template <std::ranges::input_range rng_t>
        requires(std::convertible_to<std::ranges::range_reference_t<rng_t>, char>)
    void write_field_aux(rng_t && range) { to_derived()->it->write_range(range); }

    //!\brief Write alphabet ranges.
    template <std::ranges::input_range rng_t>
        requires(detail::deliberate_alphabet<std::ranges::range_reference_t<rng_t>>)
    void write_field_aux(rng_t && range) { to_derived()->write_field_aux(range | seqan3::views::to_char); }

    //!\brief Write CStrings.
    void write_field_aux(char const * const cstr) { write_field_aux(std::string_view{cstr}); }

    //!\brief Write raw data.
    void write_field_aux(std::span<std::byte const> const range)
    {
        std::string_view const v{range.data(), range.size()};
        to_derived()->it->write_range(v);
    }

    //!\brief Write numbers.
    void write_field_aux(seqan3::arithmetic auto number) { to_derived()->it->write_number(number); }

    //!\brief Write bool.
    void write_field_aux(bool)
    {
        // TODO fix me
    }
    //!\}

    /*!\name Writing individual fields (step 2)
     * \{
     */
    //!\brief Various types have sane default implementations.
    template <field field_id>
    void write_field(vtag_t<field_id> /**/, auto & field)
    {
        if constexpr (requires(derived_t & d) { d.write_field_aux(field); })
        {
            to_derived()->write_field_aux(field);
        }
        else // no handler
        {
            static_assert(seqan3::arithmetic<decltype(field)> /*always false*/,
                          "Format X does not know how to write field Y of type Z. Provide different traits or a "
                          "custom format handler.");
            // TODO replace X Y and Z with actual strings generated from types.
        }
    }
    //!\}

    /*!\name Writing the record (step 1)
     * \{
     */
    //     template <typename field_types, typename field_ids>
    //     void write_record(seqan3::record<field_types, field_ids> && record)
    //     {
    //         // derived classes need to implement this as a public member
    //     }
    //!\}

    /*!\name Constructors, destructor and assignment.
     * \brief These are all private to prevent wrong instantiation.
     * \{
     */
    format_output_handler_base()                                               = delete;  //!< Deleted.
    format_output_handler_base(format_output_handler_base const &)             = delete;  //!< Deleted.
    format_output_handler_base(format_output_handler_base &&)                  = default; //!< Defaulted.
    ~format_output_handler_base()                                              = default; //!< Defaulted.
    format_output_handler_base & operator=(format_output_handler_base const &) = delete;  //!< Deleted.
    format_output_handler_base & operator=(format_output_handler_base &&)      = default; //!< Defaulted.

    //!\brief Construct from a std::istream.
    format_output_handler_base(std::ostream & str) : stream{&str}, it{str} {}
    //!\}
};

} // namespace bio
