// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::writer_base and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/detail/misc_output.hpp>
#include <seqan3/io/detail/out_file_iterator.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/format/output_format_handler_base.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/stream/transparent_ostream.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/utility/type_list/traits.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// writer_base
// ----------------------------------------------------------------------------

// TODO document format handling
template <typename options_t>
class writer_base
{
public:
    /*!\name Field types and record type
     * \brief These types size_Tare relevant for record/row-based reading; they may be manipulated via the \ref
     * traits_type to achieve different storage behaviour.
     * \{
     */
    //!\brief The field IDs encoded as a type.
    using field_ids = std::remove_cvref_t<decltype(options_t::field_ids)>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator  = detail::out_file_iterator<writer_base>;
    //!\brief The type returned by end().
    using sentinel  = std::default_sentinel_t;
    //!\}

    // protected:
    //!\privatesection
    /*!\name Format handling
     * \{
     */
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats = std::remove_cvref_t<decltype(options_t::formats)>;

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type = detail::transfer_template_args_onto_t<valid_formats, std::variant>;
    using format_handler_type =
      detail::transfer_template_args_onto_t<list_traits::transform<output_format_handler, valid_formats>, std::variant>;
    //!\}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    writer_base()                    = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    writer_base(writer_base const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    writer_base & operator=(writer_base const &) = delete;
    //!\brief Move construction is defaulted.
    writer_base(writer_base &&)                  = default;
    //!\brief Move assignment is defaulted.
    writer_base & operator=(writer_base &&) = default;
    //!\brief Destructor is defaulted.
    ~writer_base()                          = default;

    /*!\brief Construct from filename.
     * \param[in] filename  Path to the file you wish to open.
     * \param[in] frmt      The file format given as e.g. `format_fasta{}` [optional]
     * \param[in] opt       Writer options (exact type depends on specialisation). [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant, non-readable, unknown format.
     *
     * \details
     *
     * By default the format is determined from the extension of the provided filename. It may instead be specified
     * manually as specific type (e.g. `format_fasta{}`) or as a std::variant over all the allowed
     * formats.
     *
     *
     * ### Compression
     *
     * This constructor transparently applies a compression stream on top of the file stream in case
     * the extension indicates that you want compression.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    writer_base(std::filesystem::path const & filename, format_type const & fmt, options_t const & opt = options_t{}) :
      options{opt}, stream{filename, opt.stream_options}, format{fmt}
    {
        init();
    }

    //!\overload
    explicit writer_base(std::filesystem::path const & filename, options_t const & opt = options_t{}) :
      options{opt}, stream{filename, opt.stream_options}
    {
        // initialise format handler or throw if format is not found
        detail::set_format(format, stream.truncated_filename());

        init();
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \param[in] str  The stream to operate on.
     * \param[in] frmt The file format given as e.g. `format_fasta{}`.
     * \param[in] opt  Writer options (exact type depends on specialisation). [optional]
     *
     * \details
     *
     * ### Decompression
     *
     * This constructor transparently applies a compression stream on top of the file stream in case
     * the extension indicates that you want compression.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    writer_base(std::ostream & str, format_type const & frmt, options_t const & opt = options_t{}) :
      options{opt}, stream{str, opt.stream_options}, format{frmt}
    {
        init();
    }

    //!\overload
    template <typename temporary_stream_t>
        requires(movable_ostream<temporary_stream_t> && !std::is_lvalue_reference_v<temporary_stream_t>)
    writer_base(temporary_stream_t && str, format_type const & frmt, options_t const & opt = options_t{}) :
      options{opt}, stream{std::move(str), opt.stream_options}, format{frmt}
    {
        init();
    }
    //!\}

    /*!\name Range interface
     * \brief Provides functions for record based writing of the file.
     * \{
     */
    /*!\brief Returns an iterator to current position in the file.
     * \returns An iterator pointing to the current position in the file.
     *
     * You can write to the file by assigning to the iterator, but using push_back() is usually more intuitive.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_range_interface.cpp
     */
    iterator begin() noexcept { return {*this}; }

    /*!\brief Returns a sentinel for comparison with iterator.
     * \returns An end that is never reached.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour. It
     * always compares false against an iterator.
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

    /*!\brief Write a seqan3::record to the file.
     * \tparam field_types Types of the fields in the record.
     * \tparam field_ids   IDs of the fields in the record.
     * \param[in] r        The record to write.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant. TODO linear in the size of the written sequences?
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_push_back_record.cpp
     */
    template <typename field_types, typename field_ids>
    void push_back(record<field_types, field_ids> const & r)
    {
        write_record(r);
    }

    //!\overload
    template <typename field_types, typename field_ids>
    void push_back(record<field_types, field_ids> & r)
    {
        write_record(r); // pass as non-const to allow parsing views that are not const-iterable
    }

    //!\overload
    template <typename field_types, typename field_ids>
    void push_back(record<field_types, field_ids> && r)
    {
        write_record(r); // pass as non-const to allow parsing views that are not const-iterable
    }

    /*!\brief Write a record to the file by passing individual fields.
     * \tparam arg_types Types of the fields.
     * \param[in] args   The fields to be written.
     *
     * \details
     *
     * Uses the field IDs specified in the options to construct a record from the given arguments. Then writes that
     * record as if calling push_back().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_emplace_back.cpp
     */
    template <typename... arg_types>
    void emplace_back(arg_types &&... args)
    {
        static_assert(sizeof...(arg_types) == field_ids::size, "Wrong number of arguments provided to emplace_back()");
        push_back(tie_record<field_ids>(args...));
    }

    /*!\brief Write a range of records to the file.
     * \tparam rng_t     Type of the range, must satisfy std::ranges::output_range and have a reference type that
     *                   satisfies seqan3::tuple_like.
     * \param[in] range  The range to write.
     *
     * \details
     *
     * This function simply iterates over the argument and calls push_back() on each element.
     *
     * ### Complexity
     *
     * Linear in the number of records.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_batch_write.cpp
     */
    template <std::ranges::input_range rng_t>
    writer_base & operator=(rng_t && range)
    //!\cond
    //         requires tuple_like<std::ranges::range_reference_t<rng_t>>
    //!\endcond
    {
        for (auto && record : range)
            push_back(std::forward<decltype(record)>(record));
        return *this;
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy std::ranges::input_range and have a reference type that
     *                   satisfies seqan3::tuple_like.
     * \param[in] range  The range to write.
     * \param[in] f      The file being written to.
     *
     * \details
     *
     * This operator enables sequence_file_output to be at the end of a piping operation. It just calls
     * operator=() internally.
     *
     * ### Complexity
     *
     * Linear in the number of records.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_batch_write.cpp
     *
     * This is especially useful in combination with file-based filters:
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_view_pipeline.cpp
     */
    template <std::ranges::input_range rng_t>
    friend writer_base & operator|(rng_t && range, writer_base & f)
    //!\cond
    //         requires tuple_like<std::ranges::range_reference_t<rng_t>>
    //!\endcond
    {
        f = range;
        return f;
    }

    //!\overload
    template <std::ranges::input_range rng_t>
    friend writer_base operator|(rng_t && range, writer_base && f)
    //!\cond
    //         requires tuple_like<std::ranges::range_reference_t<rng_t>>
    //!\endcond
    {
#if defined(__GNUC__) && (__GNUC__ == 9) // an unreported build problem of GCC9
        for (auto && record : range)
            f.push_back(std::forward<decltype(record)>(record));
#else // ^^^ workaround | regular solution ↓↓↓
        f = range;
#endif
        return std::move(f);
    }
    //!\}

protected:
    //!\privatesection

    void init()
    {
        // set format-handler
        std::visit([&](auto f) { format_handler = output_format_handler<decltype(f)>{stream, options}; }, format);
    }

    void write_record(auto & r)
    {
        init_state = false;
        std::visit([&r](auto & handler) { handler.write_record(r); }, format_handler);
    }

    //!\brief The options.
    options_t                 options;
    //!\brief The stream.
    transparent_ostream<char> stream;

    //!\brief True as long as no IO has happened, yet.
    bool init_state = true;

    //!\brief The std::variant holding the detected/selected format.
    format_type         format;
    //!\brief The std::variant holding the respective handler.
    format_handler_type format_handler;

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::writer_base
 * \{
 */

//!\brief Deduces the sequence input file type from the stream and the format.
// template <typename options_type>
// writer_base(std::filesystem::path, options_type)
//     -> writer_base<options_type>;
//
// //!\brief Deduces the sequence input file type from the stream, the format and the field ids.
// template <typename options_type>
// writer_base(std::ostream && stream, options_type)
//     -> writer_base<options_type>;
//
// //!\overload
// template <typename options_type>
// writer_base(std::ostream & stream, options_type)
//     -> writer_base<options_type>;
//!\}

} // namespace seqan3
