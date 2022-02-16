// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::writer_base.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <bio/detail/misc.hpp>
#include <bio/detail/out_file_iterator.hpp>
#include <bio/exception.hpp>
#include <bio/format/format_output_handler.hpp>
#include <bio/record.hpp>
#include <bio/stream/concept.hpp>
#include <bio/stream/transparent_ostream.hpp>
#include <seqan3/utility/type_list/traits.hpp>

namespace bio
{

// ----------------------------------------------------------------------------
// writer_base
// ----------------------------------------------------------------------------

/*!\brief This is a (non-CRTP) base-class for file writers.
 * \tparam options_t Type of the reader options.
 * \details
 *
 * Most file writer inherit from this class to reduce implementation overhead. It is not relevant for most users
 * of the library.
 */
template <typename options_t>
class writer_base
{
protected:
    //!\privatesection
    /*!\name Format handling
     * \{
     */
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats            = decltype(options_t::formats);
    /*!\brief The seqan3::format_output_handler corresponding to the format(s).
     * \details
     * Metaprogramming shortcut to turn `type_list<vcf, bcf>` into
     * `std::variant<std::monostate, format_output_handler<vcf>, format_output_handler<bcf>>`.
     *
     * std::monostate is necessary, because the handlers are not default-constructible and the variant is
     * set later than construction.
     */
    using format_handler_variant_t = seqan3::detail::transfer_template_args_onto_t<
      seqan3::list_traits::concat<seqan3::type_list<std::monostate>,
                                  seqan3::list_traits::transform<format_output_handler, valid_formats>>,
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
     * \{
     */
    //!\brief The iterator type of this view (an input iterator).
    using iterator = detail::out_file_iterator<writer_base>;
    //!\brief The type returned by end().
    using sentinel = std::default_sentinel_t;
    //!\}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    writer_base()                    = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    writer_base(writer_base const &) = delete;
    //!\brief Move construction is defaulted.
    writer_base(writer_base &&)      = default;
    //!\brief Destructor which can potentially throw.
    ~writer_base() noexcept(false)
    {
        /* Implementation note:
         *
         * std::variant's destructor is always marked noexcept(true), so
         * when a format_handler throws on destruction, the program immediately
         * terminates. This is not desirable, because the exception cannot
         * be caught, and the functionality cannot be tested.
         *
         * This piece of code creates a temporary of the same type as the handler
         * that is then move-constructed from the handler, so it destructs (and
         * potentially throws) before the variant is itself destructed.
         * Thus the exception cannot be propagated outside of the destructor as usual.
         */
        auto del = []<typename handler_t>(handler_t & handler) { [[maybe_unused]] handler_t tmp = std::move(handler); };

        std::visit(del, format_handler);
    }
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    writer_base & operator=(writer_base const &) = delete;
    //!\brief Move assignment is defaulted.
    writer_base & operator=(writer_base &&)      = default;

    /*!\brief Construct from filename.
     * \param[in] filename  Path to the file you wish to open.
     * \param[in] fmt      The file format given as e.g. `fasta{}` [optional]
     * \param[in] opt       Writer options (exact type depends on specialisation). [optional]
     * \throws bio::file_open_error If the file could not be opened, e.g. non-existant, non-readable, unknown format.
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
     * See the section on compression and decompression for more information.
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
     * \param[in,out] str The stream to operate on.
     * \param[in] fmt     The file format given as e.g. `format_fasta{}`.
     * \param[in] opt     Writer options (exact type depends on specialisation). [optional]
     *
     * \details
     *
     * ### Decompression
     *
     * This constructor transparently applies a compression stream on top of the file stream in case
     * the extension indicates that you want compression.
     * See the section on compression and decompression for more information.
     */
    writer_base(std::ostream & str, format_type const & fmt, options_t const & opt = options_t{}) :
      options{opt}, stream{str, opt.stream_options}, format{fmt}
    {
        init();
    }

    //!\overload
    template <movable_ostream temporary_stream_t>
        //!\cond REQ
        requires(!std::is_lvalue_reference_v<temporary_stream_t>)
    //!\endcond
    writer_base(temporary_stream_t && str, format_type const & fmt, options_t const & opt = options_t{}) :
      options{opt}, stream{std::move(str), opt.stream_options}, format{fmt}
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

    /*!\brief Write a bio::record to the file.
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
     * \param[in] ids    The composition of fields; a bio::vtag over bio::field.
     * \param[in] args   The fields to be written.
     *
     * \details
     *
     * Uses the field IDs specified as the first arguments to construct a record from the following arguments.
     * Then writes that record as if calling push_back().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     */
    template <auto... field_ids>
    void emplace_back(vtag_t<field_ids...> ids, auto &&... args)
    {
        push_back(tie_record(ids, args...));
    }

    /*!\brief Write a range of records to the file.
     * \tparam rng_t     Type of the range, must be a std::ranges::input_range over bio::record.
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
     */
    template <std::ranges::input_range rng_t>
    requires seqan3::detail::template_specialisation_of<std::remove_cvref_t<std::ranges::range_reference_t<rng_t>>,
                                                        bio::record>
      writer_base & operator=(rng_t && range)
    {
        for (auto && record : range)
            push_back(std::forward<decltype(record)>(record));
        return *this;
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must be a std::ranges::input_range over bio::record.
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
     */
    template <std::ranges::input_range rng_t>
        requires seqan3::detail::template_specialisation_of<std::remove_cvref_t<std::ranges::range_reference_t<rng_t>>,
                                                            bio::record>
    friend writer_base & operator|(rng_t && range, writer_base & f)
    {
        f = range;
        return f;
    }

    //!\overload
    template <std::ranges::input_range rng_t>
        requires seqan3::detail::template_specialisation_of<std::remove_cvref_t<std::ranges::range_reference_t<rng_t>>,
                                                            bio::record>
    friend writer_base operator|(rng_t && range, writer_base && f)
    {
        f = range;
        return std::move(f);
    }
    //!\}

protected:
    //!\privatesection

    //!\brief Set the format handler.
    void init()
    {
        // set format-handler
        std::visit([&](auto f)
                   { format_handler.template emplace<format_output_handler<decltype(f)>>(stream, options); },
                   format);
    }

    //!\brief Implementation function for writing a record.
    void write_record(auto & r)
    {
        init_state = false;
        std::visit(detail::overloaded([](std::monostate) {}, [&r](auto & handler) { handler.write_record(r); }),
                   format_handler);
    }

    /*!\name State
     * \{
     */
    //!\brief True as long as no IO has happened, yet.
    bool init_state = true;

    //!\brief The options.
    options_t                options;
    //!\brief The stream.
    transparent_ostream      stream;
    //!\brief The std::variant holding the detected/selected format.
    format_type              format;
    //!\brief The std::variant holding the respective handler.
    format_handler_variant_t format_handler;

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
    //!\}
};

} // namespace bio
