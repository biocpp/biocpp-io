// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::detail::fast_istreambuf_iterator and bio::io::detail::fast_ostreambuf_iterator.
 *        and bio::io::ostreambuf iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <charconv>
#include <iterator>
#include <ranges>
#include <span>

#include <bio/io/detail/charconv.hpp>
#include <bio/io/detail/concept.hpp>

namespace bio::io::detail
{
// ============================================================================
//  fast_istreambuf_iterator
// ============================================================================

/*!\brief Functionally the same as std::basic_streambuf<char_t, traits_t_>, but exposes protected members as public.
 * \ingroup stream
 * \tparam char_t   The stream's character type.
 * \tparam traits_t The stream's traits type.
 *
 * \details
 *
 * This wrapper adds no functionality to std::basic_streambuf and is only used to expose protected members to
 * access the get and put area of the std::basic_streambuf.
 */
template <typename char_t, typename traits_t = std::char_traits<char_t>>
struct stream_buffer_exposer : public std::basic_streambuf<char_t, traits_t>
{
    //!\brief The actual stream type.
    using base_t = std::basic_streambuf<char_t, traits_t>;

    //!\cond
    // Expose protected members:
    using base_t::eback;
    using base_t::egptr;
    using base_t::gbump;
    using base_t::gptr;
    using base_t::setg;
    using base_t::underflow;

    using base_t::epptr;
    using base_t::overflow;
    using base_t::pbase;
    using base_t::pbump;
    using base_t::pptr;
    using base_t::setp;
    //!\endcond
};

/*!\brief Functionally the same as std::istreambuf_iterator, but faster.
 * \ingroup stream
 * \tparam char_t       The stream's character type.
 * \tparam traits_t     The stream's traits type.
 *
 * \details
 *
 * Performs less virtual function calls than std::istreambuf_iterator.
 *
 * \todo Make this move-only after input iterators are allowed to be move-only.
 *
 */
template <typename char_t, typename traits_t = std::char_traits<char_t>>
class fast_istreambuf_iterator
{
private:
    //!\brief Down-cast pointer to the stream-buffer.
    stream_buffer_exposer<char_t, traits_t> * stream_buf = nullptr;

public:
    /*!\name Associated types
     * \{
     */
    using difference_type   = ptrdiff_t;               //!< Defaults to ptrdiff_t.
    using value_type        = char_t;                  //!< The char type of the stream.
    using reference         = char_t;                  //!< The char type of the stream.
    using pointer           = void;                    //!< Has no pointer type.
    using iterator_category = std::input_iterator_tag; //!< Pure input iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    fast_istreambuf_iterator() noexcept                                             = default; //!< Defaulted.
    fast_istreambuf_iterator(fast_istreambuf_iterator const &) noexcept             = default; //!< Defaulted.
    fast_istreambuf_iterator(fast_istreambuf_iterator &&) noexcept                  = default; //!< Defaulted.
    ~fast_istreambuf_iterator() noexcept                                            = default; //!< Defaulted.
    fast_istreambuf_iterator & operator=(fast_istreambuf_iterator const &) noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator & operator=(fast_istreambuf_iterator &&) noexcept      = default; //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit fast_istreambuf_iterator(std::basic_streambuf<char_t, traits_t> & ibuf) :
      stream_buf{reinterpret_cast<stream_buffer_exposer<char_t, traits_t> *>(&ibuf)}
    {
        assert(stream_buf != nullptr);
        stream_buf->underflow(); // ensure the stream buffer has content on construction
    }

    //!\brief Construct from a stream.
    explicit fast_istreambuf_iterator(std::basic_istream<char_t, traits_t> & istr) :
      fast_istreambuf_iterator{*istr.rdbuf()}
    {}
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advance by one and rebuffer if necessary (vtable lookup iff rebuffering).
    fast_istreambuf_iterator & operator++()
    {
        assert(stream_buf != nullptr);
        if ((stream_buf->gptr() + 1) == stream_buf->egptr())
            stream_buf->snextc(); // move right, then underflow()
        else
            stream_buf->gbump(1);
        return *this;
    }

    //!\overload
    void operator++(int) { ++(*this); }
    //!\}

    //!\brief Read current value from buffer (no vtable lookup, safe if not at end).
    reference operator*() const
    {
        assert(stream_buf != nullptr);
        return *stream_buf->gptr();
    }

    //!\brief Skip n characters in the input stream (works on non-seekable streams).
    void skip_n(size_t const n)
    {
        ptrdiff_t todo      = n;
        ptrdiff_t available = stream_buf->egptr() - stream_buf->gptr();
        while (todo > 0)
        {
            if (todo <= available)
            {
                stream_buf->gbump(todo);
                todo = 0;
            }
            else
            {
                todo -= available;

                stream_buf->gbump(available);
                stream_buf->underflow();
                available = stream_buf->egptr() - stream_buf->gptr();

                if (available == 0)
                {
                    throw unexpected_end_of_input{"Trying to read ",
                                                  n,
                                                  " characters, but only ",
                                                  n - todo,
                                                  " were available."};
                }
            }
        }
    }

    /*!\brief Read n characters from the stream.
     * \param[in] n Number of characters to read.
     * \param[out] out The place to write to [must hold at least n bytes capacity].
     */
    void read_n_chars_into(size_t const n, char_t * out)
    {
        /* TODO
         * This duplicates certain logic with the line reader.
         * At some point this should be unified.
         */
        ptrdiff_t todo      = n;
        ptrdiff_t available = stream_buf->egptr() - stream_buf->gptr();
        while (todo > 0)
        {
            if (todo <= available)
            {
                std::ranges::copy_n(stream_buf->gptr(), todo, out);
                stream_buf->gbump(todo);
                todo = 0;
            }
            else
            {
                std::ranges::copy_n(stream_buf->gptr(), available, out);
                todo -= available;

                stream_buf->gbump(available);
                stream_buf->underflow();
                available = stream_buf->egptr() - stream_buf->gptr();

                if (available == 0)
                {
                    throw unexpected_end_of_input{"Trying to read ",
                                                  n,
                                                  " characters, but only ",
                                                  n - todo,
                                                  " were available."};
                }
            }
        }
    }

    //!\brief Read a POD datastructure bytewise from the stream.
    template <typename t>
        requires(std::is_trivially_copyable_v<t> && !std::ranges::range<t>)
    void read_as_binary(t & pod)
    {
        // TODO enforce little endian on numbers
        read_n_chars_into(sizeof(t), reinterpret_cast<char *>(&pod));
    }

    /*!\name Comparison operators
     * \brief We define comparison only against the sentinel.
     * \{
     */
    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(fast_istreambuf_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        assert(lhs.stream_buf != nullptr);
        // compare size of remaining buffer; since ++ always resizes if possible, safe to compare pointers here
        return (lhs.stream_buf->gptr() == lhs.stream_buf->egptr());
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(fast_istreambuf_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::default_sentinel);
    }

    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(std::default_sentinel_t const &, fast_istreambuf_iterator const & rhs) noexcept
    {
        return rhs == std::default_sentinel;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(std::default_sentinel_t const &, fast_istreambuf_iterator const & rhs) noexcept
    {
        return !(rhs == std::default_sentinel);
    }
    //!\}
};

/*!\brief Functionally the same as std::ostreambuf_iterator, but offers writing a range more efficiently.
 * \ingroup stream
 * \tparam char_t       The stream's character type.
 * \tparam traits_t     The stream's traits type.
 *
 * \details
 *
 * The functions bio::io::fast_ostreambuf_iterator::write_range and bio::io::fast_ostreambuf_iterator::write_n allow
 * more efficient writing of ranges by writing in chunks that avoiding overflow checks.
 *
 * \include test/snippet/io/detail/iterator_write_range.cpp
 */
template <typename char_t, typename traits_t = std::char_traits<char_t>>
class fast_ostreambuf_iterator
{
private:
    //!\brief Down-cast pointer to the stream-buffer.
    stream_buffer_exposer<char_t, traits_t> * stream_buf = nullptr;
    //!\brief Buffer used for certain arithmetic conversions.
    std::array<char, 150>                     buffer{};

public:
    /*!\name Associated types
     * \{
     */
    using difference_type   = ptrdiff_t;                //!< Defaults to ptrdiff_t.
    using value_type        = char_t;                   //!< The char type of the stream.
    using reference         = char_t;                   //!< The char type of the stream.
    using pointer           = void;                     //!< Has no pointer type.
    using iterator_category = std::output_iterator_tag; //!< Pure output iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    fast_ostreambuf_iterator() noexcept                                             = default; //!< Defaulted.
    fast_ostreambuf_iterator(fast_ostreambuf_iterator const &) noexcept             = default; //!< Defaulted.
    fast_ostreambuf_iterator(fast_ostreambuf_iterator &&) noexcept                  = default; //!< Defaulted.
    ~fast_ostreambuf_iterator() noexcept                                            = default; //!< Defaulted.
    fast_ostreambuf_iterator & operator=(fast_ostreambuf_iterator const &) noexcept = default; //!< Defaulted.
    fast_ostreambuf_iterator & operator=(fast_ostreambuf_iterator &&) noexcept      = default; //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit fast_ostreambuf_iterator(std::basic_streambuf<char_t, traits_t> & ibuf) :
      stream_buf{reinterpret_cast<stream_buffer_exposer<char_t, traits_t> *>(&ibuf)}
    {
        assert(stream_buf != nullptr);
        if (stream_buf->pptr() == stream_buf->epptr())
            stream_buf->overflow(); // ensures that put area has space available
    }

    //!\brief Construct from a stream buffer.
    explicit fast_ostreambuf_iterator(std::basic_ostream<char_t, traits_t> & ostr) :
      fast_ostreambuf_iterator{*ostr.rdbuf()}
    {}
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief no op.
    fast_ostreambuf_iterator & operator++() { return *this; }
    //!\overload
    fast_ostreambuf_iterator & operator++(int) { return *this; }
    //!\}

    //!\brief no op.
    fast_ostreambuf_iterator & operator*() { return *this; }

    //!\brief no op.
    fast_ostreambuf_iterator * operator->() { return this; }

    //!\brief Writes a character to the associated output stream.
    fast_ostreambuf_iterator & operator=(char_t const c)
    {
        assert(stream_buf != nullptr);
        // TODO: evaluate whether this is actually faster than just calling: stream_buf->sputc(c);
        if (stream_buf->pptr() == stream_buf->epptr())
        {
            if (stream_buf->sputc(c) == traits_t::eof()) // sputc() [virtual], then write character
            {
                // LCOV_EXCL_START
                throw std::ios_base::failure{"Cannot write to output stream (reached traits::eof() condition)."};
                // LCOV_EXCL_STOP
            }
        }
        else
        {
            *stream_buf->pptr() = c;
            stream_buf->pbump(1); // advance pptr() in put area without any checks
        }

        return *this;
    }

    //!\brief Returns `true if this iterator has encountered the end-of-file condition on output, `false` otherwise.
    bool failed() const noexcept { return stream_buf->overflow() == traits_t::eof(); }

    /*!\brief Writes a range to the associated output.
     * \tparam rng_t The type of range to write; Must model std::ranges::forward_range.
     * \param[in] rng The range to write.
     * \returns If `rng_t` models `std::ranges::borrowed_range` returns an iterator pointing to end of the range
     *          (rng) else returns `void`.
     *
     * This function avoids the buffer-at-end check by writing the range in chunks, where a chunks has the size of
     * the remaining space in the put area of the buffer.
     * If the range type models `std::ranges::sized_range` the chunks are written using `std::ranges::copy_n`, which
     * may use memcpy if applicable. Otherwise, a simple for loop iterates over the chunk.
     *
     * \attention You can only use the return value (end iterator) if your range type models
     *            `std::ranges::borrowed_range`.
     *
     * Example:
     *
     * \include test/snippet/io/detail/iterator_write_range.cpp
     */
    template <std::ranges::forward_range rng_t>
    auto write_range(rng_t && rng)
    {
        using sen_t = std::ranges::sentinel_t<rng_t>;
        using it_t  = std::ranges::iterator_t<rng_t>;

        it_t  it  = std::ranges::begin(rng);
        sen_t end = std::ranges::end(rng);

        if (stream_buf->epptr() - stream_buf->pptr() == 0 && it != end)
        {
            stream_buf->sputc(*it);
            ++it;
        }

        while (it != end)
        {
            size_t const buffer_space = stream_buf->epptr() - stream_buf->pptr();
            assert(buffer_space > 0);

            if constexpr (std::sized_sentinel_for<sen_t, it_t>)
            {
                size_t const characters_to_write = std::min<size_t>(end - it, buffer_space);
                auto         copy_res            = std::ranges::copy_n(it, characters_to_write, stream_buf->pptr());
                it                               = copy_res.in;
                stream_buf->pbump(characters_to_write);
            }
            else
            {
                size_t i = 0;
                for (; it != end && i < buffer_space; ++it, ++i)
                    *stream_buf->pptr() = *it;
                stream_buf->pbump(i);
            }

            if (it == end) // no more characters to write
                break;

            // Push one more character and flush
            if (stream_buf->sputc(*it) == traits_t::eof())
            {
                // LCOV_EXCL_START
                throw std::ios_base::failure{"Cannot write to output stream (reached traits::eof() condition)."};
                // LCOV_EXCL_STOP
            }

            ++it; // drop 1 character that has been written in sputc() above
        }

        if constexpr (std::ranges::borrowed_range<rng_t>)
            return it;
        else
            return;
    }

    //!\overload
    template <std::ranges::contiguous_range rng_t>
        //!\cond
        requires(std::ranges::sized_range<rng_t> && std::same_as<std::ranges::range_value_t<rng_t> const, char const>)
    //!\endcond
    auto write_range(rng_t && rng)
    {
        stream_buf->sputn(std::ranges::data(rng), std::ranges::size(rng));

        if constexpr (std::ranges::borrowed_range<rng_t>)
            return std::ranges::begin(rng) + std::ranges::size(rng);
        else
            return;
    }

    //!\overload
    void write_range(char const * const cstring) { write_range(std::string_view{cstring}); }

    /*!\brief Writes a number to the underlying stream buffer using std::to_chars.
     * \param[in] num The number to write.
     */
    void write_number(meta::arithmetic auto num)
    {
        if (stream_buf->epptr() - stream_buf->pptr() > 150) // enough space for any number, should be likely
        {
            auto res = to_chars(stream_buf->pptr(), stream_buf->epptr(), num);
            stream_buf->pbump(res.ptr - stream_buf->pptr()); // advance pptr
        }
        else
        {
            auto res = to_chars(&buffer[0], &buffer[0] + sizeof(buffer), num);
            write_range(std::span<char>{&buffer[0], res.ptr});
        }
    }

    //!\brief Write the binary representation of an object byte-wise to the output stream.
    template <typename t>
        requires(std::is_trivially_copyable_v<t> && !std::ranges::range<t>)
    void write_as_binary(t const & num)
    {
        // TODO enforce little endian on numbers
        std::string_view v{reinterpret_cast<char const *>(&num), sizeof(num)};
        write_range(v);
    }

    //!\brief Write the binary representation of a contiguous range byte-wise to the output stream.
    template <std::ranges::contiguous_range rng_t>
        requires std::ranges::sized_range<rng_t>
    void write_as_binary(rng_t && rng)
    {
        // TODO enforce little endian on elements?
        std::string_view v{reinterpret_cast<char const *>(std::ranges::data(rng)),
                           reinterpret_cast<char const *>(std::ranges::data(rng) + std::ranges::size(rng))};
        write_range(v);
    }

    /*!\brief Write `"\n"` or `"\r\n"` to the stream buffer, depending on arguments.
     * \param add_cr Whether to add carriage return, too.
     * \ingroup io
     */
    void write_end_of_line(bool const add_cr)
    {
        if (add_cr)
            *this = '\r';
        *this = '\n';
    }
};

} // namespace bio::io::detail
