// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::plain_io::writer.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <ranges>
#include <string_view>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/views/to_char.hpp>

#include <bio/plain_io/misc.hpp>
#include <bio/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/stream/transparent_ostream.hpp>

namespace bio::plain_io
{

/*!\brief A concept that describes types directly writable via plain_io
 * \ingroup plain_io
 *
 * \details
 *
 * The following types model this concept:
 *
 *  * `char`
 *  * Any number type (floating point and integral). Note that `signed char` und `unsigned char` are interpreted as
 * numbers.
 *  * String literals ("foobar").
 *  * Any std::ranges::input_range whose element type is convertible to `char`.
 *  * Any std::ranges::input_range whose element type is a seqan3::alphabet; seqan3::to_char is called per element.
 */
template <typename arg_t>
concept writable_to_output =
  (std::same_as<std::remove_cvref_t<arg_t>, char> || seqan3::arithmetic<std::remove_cvref_t<arg_t>> ||
   std::same_as<std::decay_t<arg_t>, char const *> ||
   (std::ranges::input_range<arg_t> && std::convertible_to<std::ranges::range_reference_t<arg_t>, char>) ||
   (std::ranges::input_range<arg_t> && seqan3::alphabet<std::ranges::range_reference_t<arg_t>>));

} // namespace bio::plain_io

namespace bio::detail
{

//!\brief An output iterator for writing into plaintext files via certain convenience functions.
//!\ingroup plain_io
template <plain_io::record_kind record_kind_ = plain_io::record_kind::line>
class plaintext_output_iterator
{
private:
    //!\brief The stream iterator.
    bio::detail::fast_ostreambuf_iterator<char> stream_it;

    //!\brief Delimiter between fields.
    char field_sep  = '\t';
    //!\brief Delimiter between records [not exposed to modification ATM].
    char record_sep = '\n';
    //!\brief Whether to add windows EOL.
    bool add_CR     = false;

    //!\brief Write a single argument to the stream.
    void write_single(plain_io::writable_to_output auto && arg)
    {
        using arg_t = decltype(arg);
        if constexpr (std::same_as<std::remove_cvref_t<arg_t>, char>)
        {
            *stream_it = arg;
        }
        else if constexpr (seqan3::arithmetic<std::remove_cvref_t<arg_t>>)
        {
            stream_it->write_number(arg);
        }
        else if constexpr (std::same_as<std::decay_t<arg_t>, char const *>)
        {
            stream_it->write_range(std::string_view{arg});
        }
        else if constexpr (std::ranges::input_range<arg_t> &&
                           std::convertible_to<std::ranges::range_reference_t<arg_t>, char>)
        {
            stream_it->write_range(arg);
        }
        else // if constexpr (std::ranges::input_range<arg_t> && alphabet<std::ranges::range_reference_t<arg_t>>)
        {
            stream_it->write_range(arg | seqan3::views::to_char);
        }
    }

    //!\brief Write the EOL character(s).
    void write_end_of_line() { stream_it->write_end_of_line(add_CR); }

public:
    /*!\name Associated types
     * \{
     */
    using difference_type   = ptrdiff_t;                //!< Defaults to ptrdiff_t.
    using value_type        = void;                     //!< The char type of the stream.
    using reference         = void;                     //!< The char type of the stream.
    using pointer           = void;                     //!< Has no pointer type.
    using iterator_category = std::output_iterator_tag; //!< Pure output iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    plaintext_output_iterator() noexcept                                              = default; //!< Defaulted.
    plaintext_output_iterator(plaintext_output_iterator const &) noexcept             = default; //!< Defaulted.
    plaintext_output_iterator(plaintext_output_iterator &&) noexcept                  = default; //!< Defaulted.
    ~plaintext_output_iterator() noexcept                                             = default; //!< Defaulted.
    plaintext_output_iterator & operator=(plaintext_output_iterator const &) noexcept = default; //!< Defaulted.
    plaintext_output_iterator & operator=(plaintext_output_iterator &&) noexcept      = default; //!< Defaulted.

    //!\brief Construct from a stream.
    explicit plaintext_output_iterator(std::ostream & ostr) : stream_it{ostr} {}

    //!\brief Construct from a stream and a seperator.
    plaintext_output_iterator(std::ostream & ostr,
                              char const     sep) requires(record_kind_ == plain_io::record_kind::line_and_fields) :
      stream_it{ostr}, field_sep{sep}
    {}
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief no op.
    plaintext_output_iterator & operator++() { return *this; }
    //!\overload
    plaintext_output_iterator & operator++(int) { return *this; }
    //!\}

    /*!\name Dereference operators
     * \{
     */
    //!\brief Dereference.
    plaintext_output_iterator & operator*() { return *this; }

    //!\brief Dereference.
    plaintext_output_iterator * operator->() { return this; }
    //!\}

    /*!\name Assignment operators
     * \{
     */
    //!\brief Writes a character to the associated output stream.
    plaintext_output_iterator & operator=(std::string_view const line) requires(record_kind_ ==
                                                                                plain_io::record_kind::line)
    {
        write_line(line);
        return *this;
    }

    //!\brief Writes a character to the associated output stream.
    template <typename range_of_fields_t>
        requires((record_kind_ == plain_io::record_kind::line_and_fields) &&
                 std::ranges::input_range<range_of_fields_t> &&
                 plain_io::writable_to_output<std::ranges::range_reference_t<range_of_fields_t>>)
    plaintext_output_iterator & operator=(range_of_fields_t && range_of_fields)
    {
        write_range_as_fields(std::forward<range_of_fields_t>(range_of_fields));
        return *this;
    }

    //!\brief Writes a character to the associated output stream.
    plaintext_output_iterator & operator=(plain_io::record const & record) requires(
      record_kind_ == plain_io::record_kind::line_and_fields)
    {
        write_range_as_fields(record.fields);
        return *this;
    }
    //!\}

    /*!\name Write functions.
     * \brief Write output to the stream.
     * \{
     */
    //!\brief Write one or more arguments to the file (no trailing newline).
    void write(plain_io::writable_to_output auto && head, plain_io::writable_to_output auto &&... tail)
    {
        write_single(std::forward<decltype(head)>(head));
        if constexpr (sizeof...(tail) > 0)
        {
            if constexpr (record_kind_ == plain_io::record_kind::line)
                (write_single(std::forward<decltype(tail)>(tail)), ...);
            else
                ((write_single(field_sep), write_single(std::forward<decltype(tail)>(tail))), ...);
        }
    }

    //!\brief Write zero or more arguments to the file plus a no trailing newline.
    void write_line(plain_io::writable_to_output auto &&... args)
    {
        if constexpr (sizeof...(args) > 0)
            write(std::forward<decltype(args)>(args)...);
        write_end_of_line();
    }

    //!\brief Write all elements from the given range to the file.
    template <std::ranges::input_range range_of_fields_t>
        requires((record_kind_ == plain_io::record_kind::line_and_fields) &&
                 plain_io::writable_to_output<std::ranges::range_reference_t<range_of_fields_t>>)
    void write_range_as_fields(range_of_fields_t && range_of_fields)
    {
        auto it = std::ranges::begin(range_of_fields);
        auto e  = std::ranges::end(range_of_fields);
        if (it != e)
        {
            write_single(*it);
            ++it;
        }

        while (it != e)
        {
            write_single(field_sep);
            write_single(*it);
            ++it;
        }

        write_end_of_line();
    }
    //!\}

    //!\brief Add carriage return characters before the linefeed (THIS IS NOT RECOMMENDED).
    void add_carriage_return(bool const add) { add_CR = add; }

    /*!\name Comparison operators
     * \brief We define comparison only against the sentinel.
     * \{
     */
    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend constexpr bool operator==(plaintext_output_iterator const &, std::default_sentinel_t const &) noexcept
    {
        return false;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend constexpr bool operator!=(plaintext_output_iterator const &, std::default_sentinel_t const &) noexcept
    {
        return false;
    }

    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend constexpr bool operator==(std::default_sentinel_t const &, plaintext_output_iterator const &) noexcept
    {
        return false;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend constexpr bool operator!=(std::default_sentinel_t const &, plaintext_output_iterator const &) noexcept
    {
        return false;
    }
    //!\}
};

} // namespace bio::detail

namespace bio::plain_io
{

/*!\brief Line-wise writer of plaintext files; supports transparent compression.
 * \tparam record_kind_ Whether to insert delimiters in lines (e.g. TSV files).
 * \ingroup plain_io
 *
 * \details
 *
 * ### Example
 *
 * Line-wise writing:
 * ```cpp
 * bio::plain_io::writer writer{"example.txt"};
 *
 * writer.push_back("line 1");
 * writer.emplace_back("line 2 contains the integer ", 3, " and the float ", 4.3, "!");
 * writer.push_back("line 3");
 * ```
 *
 * Column-wise writing:
 *```cpp
 * bio::plain_io::writer writer{"example.tsv.gz", '\t'}; // transparent compression supported
 *
 * writer.push_back("This is the header line");
 * writer.emplace_back("col1", 3.4, "col3", 4);          // will be printed tab-separated
 *```
 */
template <record_kind record_kind_>
class writer
{
private:
    //!\brief The element type of this range.
    using value_type = std::conditional_t<record_kind_ == record_kind::line, std::string_view, record>;

public:
    /*!\name Range associated types
     * \brief The types necessary to facilitate the behaviour of an input range (used in record-wise reading).
     * \{
     */
    //!\brief The iterator type of this view (an input iterator).
    using iterator       = detail::plaintext_output_iterator<record_kind_>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator = void;
    //!\brief The type returned by end().
    using sentinel       = std::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    writer()                           = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    writer(writer const &)             = delete;
    //!\brief Move construction is defaulted.
    writer(writer &&)                  = default;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    writer & operator=(writer const &) = delete;
    //!\brief Move assignment is defaulted.
    writer & operator=(writer &&)      = default;
    //!\brief Destructor is defaulted.
    ~writer()                          = default;

    /*!\brief Construct from filename.
     * \param[in] filename        Path that you wish to write to.
     * \param[in] field_separator Delimiter between fields in a line. [optional]
     * \param[in] ostream_options Options passed to the underlying stream; see seqan3::transparent_ostream_options.
     * [optional] \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant or non-readable.
     *
     * \details
     *
     * ### Compression
     *
     * This constructor transparently applies a compression stream on top of the file stream in case
     * the filename has an extension that indicates that. You can explicitly request a specific kind of compression
     * via the ostream_options.
     */
    writer(std::filesystem::path const &       filename,
           char const                          field_separator,
           transparent_ostream_options const & ostream_options =
             transparent_ostream_options{}) requires(record_kind_ == record_kind::line_and_fields) :
      stream{filename, ostream_options}, it{stream, field_separator}
    {}

    //!\overload
    explicit writer(std::filesystem::path const &       filename,
                    transparent_ostream_options const & ostream_options =
                      transparent_ostream_options{}) requires(record_kind_ == record_kind::line) :
      stream{filename, ostream_options}, it{stream}
    {}

    /*!\brief Construct from an existing stream and with specified format.
     * \param[in] str             The stream to open from; lvalues and rvalues are supported.
     * \param[in] field_separator Delimiter between fields in a line. [optional]
     * \param[in] ostream_options Options passed to the underlying stream; see seqan3::transparent_ostream_options.
     * [optional] \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant or non-readable.
     *
     * \details
     *
     * ### Compression
     *
     * This constructor transparently applies a compression stream on top of the stream in case you explicitly
     * request this via the ostream_options (default is no compression).
     */
    writer(std::ostream &                      str,
           char const                          field_separator,
           transparent_ostream_options const & ostream_options =
             transparent_ostream_options{}) requires(record_kind_ == record_kind::line_and_fields) :
      stream{str, ostream_options}, it{stream, field_separator}
    {}

    //!\overload
    explicit writer(std::ostream &                      str,
                    transparent_ostream_options const & ostream_options =
                      transparent_ostream_options{}) requires(record_kind_ == record_kind::line) :
      stream{str, ostream_options}, it{stream}
    {}

    //!\overload
    template <movable_ostream stream_t>
        //!\cond REQ
        requires(!std::is_lvalue_reference_v<stream_t> && record_kind_ == record_kind::line_and_fields)
    //!\endcond
    writer(stream_t &&                         str,
           char const                          field_separator,
           transparent_ostream_options const & ostream_options = transparent_ostream_options{}) :
      stream{std::move(str), ostream_options}, it{stream, field_separator}
    {}

    //!\overload
    template <movable_ostream stream_t>
        //!\cond REQ
        requires(!std::is_lvalue_reference_v<stream_t> && record_kind_ == record_kind::line)
    //!\endcond
    explicit writer(stream_t &&                         str,
                    transparent_ostream_options const & ostream_options = transparent_ostream_options{}) :
      stream{std::move(str), ostream_options}, it{stream}
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
    iterator begin() { return it; }

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
    //!\}

    /*!\name Write functions
     * \brief Functions for writing output.
     * \{
     */
    /*!\brief Write one or more values to the file without a terminating newline.
     * \tparam head_t Type of the first argument; see below for supported types.
     * \tparam tail_t Types of further arguments.
     * \param[in] head The first argument.
     * \param[in] tail Further arguments.
     *
     * \details
     *
     * See bio::plain_io::writable_to_output for a list of the types supported.
     *
     * If the file is a delimited file, the delimiter character is written between values -- otherwise they are simply
     * concatenated.
     *
     * No newline character(s) are written.
     */
    void write(writable_to_output auto && head, writable_to_output auto &&... tail)
    {
        it->write(std::forward<decltype(head)>(head), std::forward<decltype(tail)>(tail)...);
    }

    /*!\brief Write zero or more values to the file followed by a newline.
     * \tparam args_t Types of the arguments; see below for supported types.
     * \param[in] args The arguments.
     *
     * \details
     *
     * See bio::plain_io::writable_to_output for a list of the types supported.
     *
     * If the file is a delimited file, the delimiter character is written between values -- otherwise they are
     * simply concatenated.
     *
     * A newline is inserted after the arguments are written. This function can be called with zero arguments, in
     * which case only a newline is written.
     */
    void emplace_back(writable_to_output auto &&... args) { it->write_line(std::forward<decltype(args)>(args)...); }

    /*!\brief Write a line to the file.
     * \param[in] line The range of fields.
     *
     * \details
     *
     * This function is only available for non-delimited files.
     *
     * A newline is inserted after the argument is written.
     */
    void push_back(std::string_view const line) requires(record_kind_ == record_kind::line) { *it = line; }

    /*!\brief Write a range of fields separated by the delimiter the file.
     * \tparam range_of_fields_t A std::ranges::input_range.
     * \param[in] range_of_fields The range of fields.
     *
     * \details
     *
     * This function is only available for delimited files and the delimiter is inserted between each value.
     *
     * The element type of the range must be one of the types described in bio::plain_io::writable_to_output.
     *
     * A newline is inserted after the arguments are written.
     */
    template <std::ranges::input_range range_of_fields_t>
        requires((record_kind_ == record_kind::line_and_fields) &&
                 writable_to_output<std::ranges::range_reference_t<range_of_fields_t>>)
    void push_back(range_of_fields_t && range_of_fields) { *it = std::forward<range_of_fields_t>(range_of_fields); }

    /*!\brief Write a record to the file.
     * \param[in] record The record to be written.
     *
     * \details
     *
     * This function is only available for delimited files.
     *
     * The same as calling push_back(record.fields).
     */
    void push_back(record const & record) requires(record_kind_ == record_kind::line_and_fields)
    {
        *it = record.fields;
    }
    //!\}

    /*!\name Assignment functions
     * \brief Functions for writing output.
     * \{
     */
    /*!\brief            Write a range of records to the file.
     * \tparam rng_t     Type of the range, must satisfy std::ranges::input_range and have a reference type that
     *                   be push_back()-ed to this type.
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
        requires(std::convertible_to<std::ranges::range_reference_t<rng_t>, value_type>)
    writer & operator=(rng_t && range)
    {
        for (auto && record : range)
            push_back(std::forward<decltype(record)>(record));
        return *this;
    }

    /*!\brief            Write a range of records to the file.
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
        requires std::convertible_to<std::ranges::range_reference_t<rng_t>, value_type>
    friend writer & operator|(rng_t && range, writer & f)
    {
        f = range;
        return f;
    }

    //!\overload
    template <std::ranges::input_range rng_t>
        requires std::convertible_to<std::ranges::range_reference_t<rng_t>, value_type>
    friend writer operator|(rng_t && range, writer && f)
    {
        f = range;
        return std::move(f);
    }
    //!\}

    /*!\name Options
     * \brief Options that can be set after the file is created.
     * \{
     */
    //!\brief Add carriage return characters before the linefeed (THIS IS NOT RECOMMENDED).
    void add_carriage_return(bool add) { it->add_carriage_return(add); }
    /* DESIGN NOTE:
     * This is explicitly not part of constructor options, because it is used very, very rarely and we
     * shouldn't add config or options just because of this.
     */
    //!\}

protected:
    //!\brief The underlying stream object.
    transparent_ostream stream;
    //!\brief The stream iterator.
    iterator            it;
};

/*!\name Deduction guides
 * \relates bio::plain_io::writer
 * \{
 */
//!\brief Deduce to line_and_fields specialisation.
writer(auto &&, char const, transparent_ostream_options const & ostream_options = transparent_ostream_options{})
  -> writer<record_kind::line_and_fields>;

//!\brief Deduce to line-only specialisation.
writer(auto &&, transparent_ostream_options const & ostream_options = transparent_ostream_options{})
  -> writer<record_kind::line>;
//!\}

} // namespace bio::plain_io
