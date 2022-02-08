// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::reader_base.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <seqan3/utility/type_list/traits.hpp>

#include <bio/detail/in_file_iterator.hpp>
#include <bio/detail/misc.hpp>
#include <bio/exception.hpp>
#include <bio/format/format_input_handler.hpp>
#include <bio/record.hpp>
#include <bio/stream/transparent_istream.hpp>

namespace bio
{

// ----------------------------------------------------------------------------
// reader_base
// ----------------------------------------------------------------------------

/*!\brief This is a (non-CRTP) base-class for I/O readers.
 * \tparam options_t Type of the reader options.
 * \details
 *
 * Most I/O readers inherit from this class to reduce implementation overhead. It is not relevant for most users
 * of the library.
 */
template <typename options_t>
class reader_base : public std::ranges::view_base
{
protected:
    //!\privatesection
    /*!\name Format handling
     * \{
     */
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats = decltype(options_t::formats);
    //!\brief The seqan3::format_input_handler corresponding to the format.
    using format_handler_type =
      seqan3::detail::transfer_template_args_onto_t<seqan3::list_traits::transform<format_input_handler, valid_formats>,
                                                    std::variant>;
    //!\}

public:
    /*!\name Format handling
     * \{
     */
    //!\brief Type of the format, a std::variant over the `valid_formats`.
    using format_type = seqan3::detail::transfer_template_args_onto_t<valid_formats, std::variant>;
    //!\brief The seqan3::format_input_handler corresponding to the format.
    //!\}

    /*!\name Field types and record type
     * \brief The exact type of the record depends on the options!
     * \{
     */
    //!\brief The type of the record, a specialisation of bio::record; acts as a tuple of the selected field types.
    using record_type = record<decltype(options_t::field_ids), decltype(options_t::field_types)>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator    = detail::in_file_iterator<reader_base>;
    //!\brief The type returned by end().
    using sentinel    = std::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    reader_base()                                = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    reader_base(reader_base const &)             = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    reader_base & operator=(reader_base const &) = delete;
    //!\brief Move construction is defaulted.
    reader_base(reader_base &&)                  = default;
    //!\brief Move assignment is defaulted.
    reader_base & operator=(reader_base &&)      = default;
    //!\brief Destructor is defaulted.
    ~reader_base()                               = default;

    /*!\brief Construct from filename.
     * \param[in] filename  Path to the file you wish to open.
     * \param[in] fmt      The file format given as e.g. `fasta{}` [optional]
     * \param[in] opt       Reader options (exact type depends on specialisation). [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant, non-readable, unknown format.
     *
     * \details
     *
     * In addition to the file name, you may fix the format and/or provide options.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on compression and decompression (TODO) for more information.
     */
    reader_base(std::filesystem::path const & filename, format_type const & fmt, options_t const & opt = options_t{}) :
      options{opt}, stream{filename, opt.stream_options}, format{fmt}
    {}

    //!\overload
    explicit reader_base(std::filesystem::path const & filename, options_t const & opt = options_t{}) :
      options{opt}, stream{filename, opt.stream_options}
    {
        // initialise format handler or throw if format is not found
        detail::set_format(format, stream.truncated_filename());
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \param[in] str  The stream to operate on.
     * \param[in] fmt The file format given as e.g. `fasta{}`. [required]
     * \param[in] opt  Reader options (exact type depends on specialisation). [optional]
     *
     * \details
     *
     * In addition to the stream, you must fix the format and you may optionally provide options.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on compression and decompression (TODO) for more information.
     */
    reader_base(std::istream & str, format_type const & fmt, options_t const & opt = options_t{}) :
      options{opt}, stream{str, opt.stream_options}, format{fmt}
    {}

    //!\overload
    template <movable_istream temporary_stream_t>
        //!\cond REQ
        requires(!std::is_lvalue_reference_v<temporary_stream_t>)
    //!\endcond
    reader_base(temporary_stream_t && str, format_type const & fmt, options_t const & opt = options_t{}) :
      options{opt}, stream{std::move(str), opt.stream_options}, format{fmt}
    {}
    //!\}

    /*!\name Range interface
     * \brief Provides functions for record based reading of the file.
     * \{
     */
    /*!\brief Returns an iterator to current position in the file.
     * \returns An iterator pointing to the current position in the file.
     * \throws seqan3::format_error
     *
     * It is safe to call this function repeatedly, but it will always return an iterator pointing to the current
     * record in the file (and not seek back to the beginning).
     *
     * Equals end() if the file is at end.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Throws seqan3::format_error if the first record could not be read into the buffer.
     */
    iterator begin()
    {
        // buffer first record
        if (init_state)
        {
            // set format-handler
            std::visit([&](auto f) { format_handler = format_input_handler<decltype(f)>{stream, options}; }, format);

            // read first record
            read_next_record();
            init_state = false;
        }

        return {*this};
    }

    /*!\brief Returns a sentinel for comparison with iterator.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    sentinel end() noexcept { return {}; }

    /*!\brief Return the record we are currently at in the file.
     * \returns A reference to the currently buffered record.
     *
     * This function returns a reference to the currently buffered record, it is identical to calling `*begin()`.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    record_type & front() noexcept { return *begin(); }
    //!\}

protected:
    //!\privatesection

    /*!\name State
     * \{
     */
    //!\brief The object holding the options.
    options_t           options;
    //!\brief The input stream.
    transparent_istream stream;
    //!\brief Buffer for a single record.
    record_type         record_buffer;
    //!\brief Tracks whether the very first record is buffered when calling begin().
    bool                init_state = true;
    //!\brief File is at position 1 behind the last record.
    bool                at_end     = false;

    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type         format;
    //!\brief The respective input handler specialisation.
    format_handler_type format_handler;
    //!\}

    //!\brief Tell the format to move to the next record and update the buffer.
    void read_next_record()
    {
        if (at_end)
            return;

        // at end if we could not read further
        if (std::istreambuf_iterator<char>{stream} == std::istreambuf_iterator<char>{})
        {
            at_end = true;
            return;
        }

        assert(!format_handler.valueless_by_exception());
        std::visit([&](auto & f) { f.parse_next_record_into(record_buffer); }, format_handler);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

} // namespace bio
