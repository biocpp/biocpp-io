// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the seqan3::output_format_handler_base.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/io/utility.hpp>

namespace bio
{

template <typename format_t>
class output_format_handler;

template <typename derived_t>
class output_format_handler_base
{
private:
    /* CRTP STUFF */
    friend derived_t;

    derived_t * to_derived() { return static_cast<derived_t *>(this); }

    derived_t const * to_derived() const { return static_cast<derived_t const *>(this); }

    /* members */
    std::ostream *                         stream = nullptr;
    detail::fast_ostreambuf_iterator<char> it;

    /*!\name Writing individual fields - defaults (step 3)
     * \{
     */
    template <std::ranges::input_range rng_t>
        requires std::same_as < std::ranges::range_reference_t<rng_t>,
    char > void write_field_aux(rng_t && range) { it->write_range(range); }

    template <std::ranges::input_range rng_t>
        requires(!arithmetic<std::ranges::range_reference_t<rng_t>> && alphabet<std::ranges::range_reference_t<rng_t>>)
    void write_field_aux(rng_t && range) { to_derived()->write_field_aux(range | views::to_char); }

    void write_field_aux(char const * const cstr) { write_field_aux(std::string_view{cstr}); }

    void write_field_aux(std::span<std::byte> const range)
    {
        std::string_view const v{range.data(), range.size()};
        it->write_range(v);
    }

    void write_field_aux(seqan3::arithmetic auto number) { it->write_number(number); }

    void write_field_aux(bool)
    {
        // TODO fix me
    }
    //!\}

    /*!\name Writing individual fields (step 2)
     * \{
     */
    //!\brief Default is no handler.
    template <field field_id, typename field_t>
    void write_field(tag_t<field_id> /**/, field_t & field)
    {
        static_assert(arithmetic<field_t> /*always false*/,
                      "Format X does not know how to write field Y of type Z. Provide different traits or a "
                      "custom format handler.");
        // TODO replace X Y and Z with actual strings generated from types.
    }

    //!\brief Various types have sane default implementations.
    template <field field_id, typename field_t>
    void write_field(tag_t<field_id> /**/,
                     field_t & field) requires(requires(derived_t & d) { d.write_field_aux(field); })
    {
        to_derived()->write_field_aux(field);
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

    output_format_handler_base()                                   = default;
    output_format_handler_base(output_format_handler_base const &) = delete;
    output_format_handler_base(output_format_handler_base &&)      = default;
    output_format_handler_base & operator=(output_format_handler_base const &) = delete;
    output_format_handler_base & operator=(output_format_handler_base &&) = default;

    output_format_handler_base(std::ostream & str) : stream{&str}, it{str} {}
};

} // namespace bio
