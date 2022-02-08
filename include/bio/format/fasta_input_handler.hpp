// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::format_input_handler implementation for bio::fasta.
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

#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

#include <bio/detail/range.hpp>
#include <bio/format/fasta.hpp>
#include <bio/format/format_input_handler.hpp>
#include <bio/plain_io/reader.hpp>

namespace bio
{

/*!\brief Format input handler for the FASTA format (bio::fasta).
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
 * | Member        | Type    | Default | Description                                      |
 * |---------------|---------|---------|--------------------------------------------------|
 * |`truncate_ids` |`bool`   | `false` | Whether to truncate IDs on the first whitespace  |
 *
 * ### Performance
 *
 * Since FASTA records can be exceptionally large, FASTA input always involves an additional buffering step.
 * Even when performing raw or lowlevel I/O (returning std::string_view or std::span<std::byte const>), all data
 * is typically copied once from the input stream buffer. This also means, that raw I/O does not preserve formatting
 * (e.g. linebreaks) as-is.
 *
 * If sequence and/or ID are requested as std::string, the record's element is swapped with the internal buffer to
 * prevent a second copy, but if output is requested as e.g. std::vector<seqan3::dna4>, a second copy needs to happen.
 * Requesting views never implies a second copy.
 */
template <>
class format_input_handler<fasta> : public format_input_handler_base<format_input_handler<fasta>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The type of the CRTP base class.
    using base_t = format_input_handler_base<format_input_handler<fasta>>;
    using base_t::parse_field;
    using base_t::parse_field_aux;
    using base_t::stream;

    //!\brief Befriend the base class to enable CRTP.
    friend base_t;
    //!\}

    //!\brief Print an error message with current line number in diagnostic.
    [[noreturn]] void error(auto const &... messages) const
    {
        std::string message = "[B.I.O. FASTA format error in line " + detail::to_string(line) + "] ";
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
    using format_fields     = vtag_t<field::id, field::seq>;
    //!\brief Type of the raw record.
    using raw_record_type   = record<format_fields, seqan3::type_list<std::string_view, std::string_view>>;
    //!\brief Type of the low-level iterator.
    using lowlevel_iterator = detail::plaintext_input_iterator<plain_io::record_kind::line>;

    //!\brief The raw record.
    raw_record_type raw_record;
    //!\brief Buffer for the ID.
    std::string     id_buffer;
    //!\brief Buffer for the sequence.
    std::string     seq_buffer;

    //!\brief The low-level iterator.
    lowlevel_iterator it;
    //!\brief A line counter.
    size_t            line = -1;

    //!\brief Read the raw record [the base class invokes this function].
    void read_raw_record()
    {
        constexpr auto not_id = !(seqan3::is_char<'>'> || seqan3::is_char<';'>);

        id_buffer.clear();
        seq_buffer.clear();
        raw_record.clear();

        /* READ ID */
        ++it;
        ++line;

        if (it == std::default_sentinel)
            error("Reached end of file while trying to read ID.");

        std::string_view current_line = *it;

        if (current_line.empty())
            error("Expected to be on begin of record but is on empty line.");

        if (not_id(current_line[0]))
            error("Record does not begin with '>' or ';'.");

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

        /* READ SEQ */
        /* Implementation NOTE: we create and compare a streambuf iterator here and do not check the line-iterator
         * (it != std::default_sentinel) && ...
         * because the line-iterator is "at end" one iteration later than the stream itself [it holds the line just
         * before the current stream-buffer pointer]. If we don't do this, the peak() might fail.
         */
        while ((std::istreambuf_iterator<char>{*stream} != std::istreambuf_iterator<char>{}) && not_id(it.peak()))
        {
            ++it;
            ++line;
            std::ranges::copy(*it | std::views::filter(!(seqan3::is_space || seqan3::is_digit)),
                              std::back_insert_iterator{seq_buffer});
        }

        if (seq_buffer.empty())
            error("No sequence or no valid sequence characters.");
        get<field::seq>(raw_record) = seq_buffer;
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
     * The options argument is typically bio::seq_io::reader_options, but any object with a subset of similarly named
     * members is also accepted. See bio::format_input_handler<bio::fasta> for the supported options and defaults.
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

} // namespace bio
