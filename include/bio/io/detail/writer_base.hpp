// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::writer_base.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <bio/meta/tag/vtag.hpp>
#include <bio/meta/type_list/traits.hpp>

#include <bio/io/detail/misc.hpp>
#include <bio/io/detail/out_file_iterator.hpp>
#include <bio/io/exception.hpp>
#include <bio/io/format/format_output_handler.hpp>
#include <bio/io/record.hpp>
#include <bio/io/stream/concept.hpp>
#include <bio/io/stream/transparent_ostream.hpp>

namespace bio::io
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
template <typename derived_t, typename options_t>
class writer_base
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief Befriend the derived type so it can instantiate.
    friend derived_t;

    //!\brief Downcast self to derived type.
    derived_t & to_derived() { return *static_cast<derived_t *>(this); }

    //!\brief Downcast self to derived type. [const-qualified version]
    derived_t const & to_derived() const { return *static_cast<derived_t const *>(this); }
    //!\}

    /*!\name Format handling
     * \{
     */
    //!\brief A bio::meta::type_list with the possible formats.
    using valid_formats            = decltype(options_t::formats);
    /*!\brief The bio::io::format_output_handler corresponding to the format(s).
     * \details
     * Metaprogramming shortcut to turn `type_list<vcf, bcf>` into
     * `std::variant<std::monostate, format_output_handler<vcf>, format_output_handler<bcf>>`.
     *
     * std::monostate is necessary, because the handlers are not default-constructible and the variant is
     * set later than construction.
     */
    using format_handler_variant_t = bio::meta::transfer_template_args_onto_t<
      meta::list_traits::concat<meta::type_list<std::monostate>,
                                meta::list_traits::transform<format_output_handler, valid_formats>>,
      std::variant>;

    //!\}

public:
    /*!\name Format handling
     * \{
     */
    //!\brief Type of the format, a std::variant over the `valid_formats`.
    using format_type = bio::meta::transfer_template_args_onto_t<valid_formats, std::variant>;
    //!\brief The bio::io::format_input_handler corresponding to the format.
    //!\}

    /*!\name Field types and record type
     * \{
     */
    //!\brief The iterator type of this view (an input iterator).
    using iterator = detail::out_file_iterator<derived_t>;
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
     * \throws bio::io::file_open_error If the file could not be opened, e.g. non-existant, non-readable, unknown
     * format.
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
     */
    iterator begin() noexcept { return {to_derived()}; }

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

    /*!\brief Write a range of records to the file.
     * \tparam rng_t     Type of the range, must be a std::ranges::input_range over bio::io::record.
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
    derived_t & operator=(rng_t && range)
      //!\cond REQ
      requires(requires { to_derived().push_back(*std::ranges::begin(range)); })
    //!\endcond
    {
        for (auto && record : range)
            to_derived().push_back(std::forward<decltype(record)>(record));
        return to_derived();
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must be a std::ranges::input_range over bio::io::record.
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
    friend derived_t & operator|(rng_t && range, derived_t & f)
      //!\cond REQ
      requires(requires { f = range; })
    //!\endcond
    {
        f = range;
        return f;
    }

    //!\overload
    template <std::ranges::input_range rng_t>
    friend derived_t operator|(rng_t && range, derived_t && f)
      //!\cond REQ
      requires(requires { f.push_back(*std::ranges::begin(range)); })
    //!\endcond
    {
        //TODO(GCC11): replace with assignment once GCC10 is dropped
        for (auto && record : range)
            f.push_back(std::forward<decltype(record)>(record));
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

} // namespace bio::io
