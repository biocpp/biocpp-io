// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides utility types for bio::plain_io::reader and bio::plain_io::writer.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <stdexcept>
#include <string_view>
#include <vector>

#include <bio/platform.hpp>

namespace bio::plain_io
{

/*!\addtogroup plain_io
 * \{
 */

/*!\brief The value type of bio::plaintext_file_input if every line is split into fields.
 * \details
 *
 * Plain I/O records are always shallow.
 */
struct record
{
    //!\brief The entire line (exluding EOL characters but including delimiters).
    std::string_view              line;
    //!\brief A range of the individual fields (without delimiters or EOL characters).
    std::vector<std::string_view> fields;
};

//!\brief Option to switch between reading-by-line and splitting a line into fields.
enum class record_kind
{
    line,           //!< Only the line is provided.
    line_and_fields //!< The line is provided and also individual fields (bio::plaintext_record).
};

/*!\brief A helper for specifying the header of a bio::plaintext_file_input.
 * \tparam record_kind Whether to split lines on delimiter (e.g. TSV files) or not.
 *
 * \details
 *
 * This class can be used similar to an enum to specify the kind of header present in a plain IO file.
 * It can be constructed from the kinds of state:
 *
 *   * bio::plain_io::header_kind::none -> no header
 *   * bio::plain_io::header_kind::first_line -> first line is treated as header
 *   * bio::plain_io::header_kind::starts_with{c} -> all lines starting with character "c" are header
 *
 * ### Example
 *
 * ```cpp
 * bio::plain_io::reader r{"example.txt"};                                               // Implicitly no header
 * bio::plain_io::reader r{"example.txt", bio::plain_io::header_kind::none};             // Explicitly no header
 * bio::plain_io::reader r{"example.txt", bio::plain_io::header_kind::first_line};       // First line
 * bio::plain_io::reader r{"example.txt", bio::plain_io::header_kind::starts_with{'#'}}; // Lines starting with '#'
 * ```
 */
class header_kind
{
private:
    //!\brief Special value for "none".
    static constexpr int none_state       = -200;
    //!\brief Special value for "first_line".
    static constexpr int first_line_state = -300;

    //!\brief The state of the variable, encoded as the value of the char or a special value.
    int state = none_state;

public:
    //!\brief The state representing "no header".
    static constexpr struct
    {
    } none{};
    //!\brief The state representing "first line is header".
    static constexpr struct
    {
    } first_line{};

    //!\brief Type of the state representing "all lines that start with character X".
    struct starts_with
    {
        //!\privatesection
        //!\brief Store the chracter.
        char c;
    };

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr header_kind() noexcept                                = default; //!< Defaulted.
    constexpr header_kind(header_kind const &) noexcept             = default; //!< Defaulted.
    constexpr header_kind(header_kind &&) noexcept                  = default; //!< Defaulted.
    constexpr header_kind & operator=(header_kind const &) noexcept = default; //!< Defaulted.
    constexpr header_kind & operator=(header_kind &&) noexcept      = default; //!< Defaulted.

    //!\brief Initialise to the "none"-state.
    constexpr header_kind(decltype(none)) noexcept : state{none_state} {}
    //!\brief Initialise to the "first_line"-state.
    constexpr header_kind(decltype(first_line)) noexcept : state{first_line_state} {}
    //!\brief Initialise to the "starts_with"-state and a given character.
    constexpr header_kind(starts_with s) noexcept : state{s.c} {}
    //!\}

    /*!\name Functions for retrieving the state.
     * \{
     */
    //!\brief Whether this object is in the "none"-state.
    constexpr bool is_none() const noexcept { return state == none_state; }
    //!\brief Whether this object is in the "first_line"-state.
    constexpr bool is_first_line() const noexcept { return state == first_line_state; }
    //!\brief Whether this object is in any "starts_with"-state.
    constexpr bool is_starts_with() const noexcept { return state != none_state && state != first_line_state; }

    /*!\brief The character stored as the "starts_with"-state.
     * \throws std::logic_error If this object is not in a "starts_with"-state.
     */
    char get_starts_with() const
    {
        if (!is_starts_with())
        {
            throw std::logic_error{
              "Tried to read starts_with from header_kind but it was in a "
              "different state."};
        }

        return static_cast<char>(state);
    }
    //!\}
};

//!\}

} // namespace bio::plain_io
