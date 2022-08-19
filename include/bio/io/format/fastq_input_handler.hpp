// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::format_input_handler implementation for bio::io::fastq.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/core/range/type_traits.hpp>

#include <bio/ranges/views/to_char.hpp>
#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

#include <bio/io/detail/range.hpp>
#include <bio/io/format/fastq.hpp>
#include <bio/io/format/format_input_handler.hpp>
#include <bio/io/plain_io/reader.hpp>

namespace bio::io
{

/*!\brief Format input handler for the FastQ format (bio::io::fastq).
 * \ingroup format
 * \details
 *
 * ### Attention
 *
 * Most users should not perform I/O through input/output handlers but should instead use the respective
 * readers/writers. See the overview (TODO link) for more information.
 *
 * ### Options
 *
 * The following options are considered if the respective member variable is availabele in the object passed to
 * the constructor:
 *
 * | Member          | Type    | Default | Description                                         |
 * |-----------------|---------|---------|-----------------------------------------------------|
 * |`truncate_ids`   |`bool`   | `false` | Whether to truncate IDs on the first whitespace     |
 *
 * ### Performance
 *
 * The current implementation is not optimised for performance and performs an unnecessary copy operation––even
 * when shallow records are returned.
 */
template <>
class format_input_handler<fastq> : public format_input_handler_base<format_input_handler<fastq>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The type of the CRTP base class.
    using base_t = format_input_handler_base<format_input_handler<fastq>>;
    using base_t::parse_field;
    using base_t::parse_field_aux;
    using base_t::stream;

    //!\brief Befriend the base class to enable CRTP.
    friend base_t;
    //!\}

    //!\brief Print an error message with current line number in diagnostic.
    [[noreturn]] void error(auto const &... messages) const
    {
        std::string message = "[BioC++ FastQ format error in line " + detail::to_string(line) + "] ";
        ((message += detail::to_string(messages)), ...);

        throw parse_error{message};
    }

    /*!\name Options
     * \{
     */
    //!\brief Whether to truncate IDs on first whitespace.
    bool truncate_ids = false;
    //!\}

    /*!\name Raw record handling
     * \{
     */
    //!\brief The fields that this format supports [the base class accesses this type].
    using format_fields     = vtag_t<field::id, field::seq, field::qual>;
    //!\brief Type of the raw record.
    using raw_record_type   = record<format_fields, seqan3::list_traits::repeat<3, std::string_view>>;
    //!\brief Type of the low-level iterator.
    using lowlevel_iterator = detail::plaintext_input_iterator<plain_io::record_kind::line>;

    //!\brief The raw record.
    raw_record_type raw_record;
    //!\brief Buffer for the ID.
    std::string     id_buffer;
    //!\brief Buffer for the sequence.
    std::string     seq_buffer;
    //!\brief Buffer for the qualities.
    std::string     qual_buffer;

    //!\brief The low-level iterator.
    lowlevel_iterator it;
    //!\brief A line counter.
    size_t            line = -1;

    //!\brief Read the raw record [the base class invokes this function].
    void read_raw_record()
    {
        id_buffer.clear();
        seq_buffer.clear();
        qual_buffer.clear();

        /* READ ID */
        {
            ++it;
            ++line;

            if (it == std::default_sentinel)
                error("Reached end of file while trying to read ID.");

            std::string_view current_line = *it;

            if (current_line.empty())
                error("Expected to be on begin of record but line is empty.");

            if ((!seqan3::is_char<'@'>)(current_line[0]))
                error("ID-line does not begin with '@'.");

            if (truncate_ids)
            {
                size_t e = 1;
                for (; e < current_line.size() && (!seqan3::is_space)(current_line[e]); ++e)
                {}
                detail::string_copy(current_line.substr(1, e - 1), id_buffer);
            }
            else
            {
                detail::string_copy(current_line.substr(1), id_buffer);
            }

            get<field::id>(raw_record) = id_buffer;
        }

        /* READ SEQ */
        {
            ++it;
            ++line;

            if (it == std::default_sentinel)
                error("Reached end of file while trying to read SEQ.");

            detail::string_copy(*it, seq_buffer); // is allowed to be empty

            get<field::seq>(raw_record) = seq_buffer;
        }

        /* READ third line */
        {
            ++it;
            ++line;

            if (it == std::default_sentinel)
                error("Reached end of file while trying to read third FastQ record line.");

            if ((!seqan3::is_char<'+'>)((*it)[0]))
                error("Third FastQ record line does not begin with '+'.");

            // we don't process the rest of the line
        }

        /* READ QUAL */
        {
            ++it;
            ++line;

            if (it == std::default_sentinel)
                error("Reached end of file while trying to read QUALITIES.");

            detail::string_copy(*it, qual_buffer); // is allowed to be empty

            get<field::qual>(raw_record) = qual_buffer;
        }

        if (size_t ssize = seq_buffer.size(), qsize = qual_buffer.size(); ssize != qsize)
            error("Size mismatch between sequence (", ssize, ") and qualities (", qsize, ").");
    }
    //!\}

    /*!\name Parsed record handling
     * \brief This is mostly done via the defaults in the base class.
     * \{
     */
    //!\brief We can prevent another copy if the user wants a string.
    void parse_field(vtag_t<field::id> const & /**/, std::string & parsed_field) { std::swap(id_buffer, parsed_field); }

    //!\brief We can prevent another copy if the user wants a string.
    void parse_field(vtag_t<field::seq> const & /**/, std::string & parsed_field)
    {
        std::swap(seq_buffer, parsed_field);
    }

    //!\brief We can prevent another copy if the user wants a string.
    void parse_field(vtag_t<field::qual> const & /**/, std::string & parsed_field)
    {
        std::swap(qual_buffer, parsed_field);
    }
    //!\}

public:
    /*!\name Constructors, destructor and assignment.
     * \{
     */
    format_input_handler()                                         = default; //!< Defaulted.
    format_input_handler(format_input_handler const &)             = delete;  //!< Deleted.
    format_input_handler(format_input_handler &&)                  = default; //!< Defaulted.
    ~format_input_handler()                                        = default; //!< Defaulted.
    format_input_handler & operator=(format_input_handler const &) = delete;  //!< Deleted.
    format_input_handler & operator=(format_input_handler &&)      = default; //!< Defaulted.

    /*!\brief Construct with an options object.
     * \param[in,out] str The input stream.
     * \param[in] options An object with options for the input handler.
     * \details
     *
     * The options argument is typically bio::io::seq_io::reader_options, but any object with a subset of similarly
     * named members is also accepted. See bio::io::format_input_handler<bio::io::fastq> for the supported options and
     * defaults.
     */
    format_input_handler(std::istream & str, auto const & options) : base_t{str}, it{str, false}
    {
        if constexpr (requires { (bool)options.truncate_ids; })
        {
            truncate_ids = options.truncate_ids;
        }
    }

    //!\brief Construct with only an input stream.
    format_input_handler(std::istream & str) : format_input_handler{str, int{}} {}
    //!\}
};

} // namespace bio::io
