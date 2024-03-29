// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::detail::in_file_iterator class template.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <cassert>
#include <ranges>

#include <bio/io.hpp>

namespace bio::io::detail
{

/*!\brief Input iterator that provides a range-like interface for readers.
 * \tparam file_type The data structure on which the iterator operates.
 * \implements std::input_Iterator
 * \ingroup io
 *
 * This iterator is a single-pass input iterator for input files. All member types are resolved
 * via `file_type`'s member types, dereference is implemented via file's `front()` member
 * function, and increment calls the `buffer_next_record()` member of file.
 *
 * Note that since this is a single-pass input iterator, post-increment returns void because
 * previous iterators are always invalid (all iterators point to the current position in single-pass
 * ranges).
 *
 * This iterator may be compared against std::default_sentinel_t, this check delegates to
 * calling the `eof()` member function on the file's stream.
 */
template <typename file_type>
class in_file_iterator
{
    static_assert(!std::is_const_v<file_type>,
                  "You cannot iterate over const files, because the iterator changes the file.");

public:
    /*!\name Member types
     * \brief The associated types are derived from the `file_type`.
     * \{
     */

    //!\brief The value type.
    using value_type        = typename file_type::record_type;
    //!\brief The reference type.
    using reference         = typename file_type::record_type &;
    //!\brief The const reference type.
    using const_reference   = typename file_type::record_type &;
    //!\brief The size type.
    using size_type         = size_t;
    //!\brief The difference type. A signed integer type, usually std::ptrdiff_t.
    using difference_type   = ptrdiff_t;
    //!\brief The pointer type.
    using pointer           = typename file_type::record_type *;
    //!\brief Tag this class as an input iterator.
    using iterator_category = std::input_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment.
     * \{
     */
    in_file_iterator()                                     = default; //!< Defaulted.
    in_file_iterator(in_file_iterator const &)             = default; //!< Defaulted.
    in_file_iterator(in_file_iterator &&)                  = default; //!< Defaulted.
    ~in_file_iterator()                                    = default; //!< Defaulted.
    in_file_iterator & operator=(in_file_iterator const &) = default; //!< Defaulted.
    in_file_iterator & operator=(in_file_iterator &&)      = default; //!< Defaulted.

    //!\brief Construct with reference to host.
    in_file_iterator(file_type & _host) noexcept : host{&_host} {}
    //!\}

    /*!\name Iterator operations
     * \{
     */
    //!\brief Move to the next record in the file and return a reference to it.
    in_file_iterator & operator++()
    {
        assert(host != nullptr);
        host->read_next_record();
        return *this;
    }

    //!\brief Post-increment is the same as pre-increment, but returns void.
    void operator++(int)
    {
        assert(host != nullptr);
        ++(*this);
    }

    //!\brief Dereference returns the currently buffered record.
    reference operator*() noexcept
    {
        assert(host != nullptr);
        return host->record_buffer;
    }

    //!\brief Dereference returns the currently buffered record.
    reference operator*() const noexcept
    {
        assert(host != nullptr);
        return host->record_buffer;
    }

    //!\brief Dereference returns the currently buffered record.
    value_type * operator->() noexcept
    {
        assert(host != nullptr);
        return &host->record_buffer;
    }

    //!\brief Dereference returns the currently buffered record.
    value_type const * operator->() const noexcept
    {
        assert(host != nullptr);
        return &host->record_buffer;
    }

    //!\}

    /*!\name Comparison operators
     * \brief Only (in-)equality comparison of iterator with end() is supported.
     * \{
     */

    //!\brief Checks whether `*this` is equal to the sentinel.
    constexpr bool operator==(std::default_sentinel_t const &) const noexcept
    {
        assert(host != nullptr);
        return host->at_end;
    }

    //!\brief Checks whether `*this` is not equal to the sentinel.
    constexpr bool operator!=(std::default_sentinel_t const &) const noexcept
    {
        assert(host != nullptr);
        return !host->at_end;
    }

    //!\brief Checks whether `it` is equal to the sentinel.
    constexpr friend bool operator==(std::default_sentinel_t const &, in_file_iterator const & it) noexcept
    {
        return (it == std::default_sentinel);
    }

    //!\brief Checks whether `it` is not equal to the sentinel.
    constexpr friend bool operator!=(std::default_sentinel_t const &, in_file_iterator const & it) noexcept
    {
        return (it != std::default_sentinel);
    }
    //!\}

private:
    //!\brief Pointer to file host.
    file_type * host{};
};

} // namespace bio::io::detail
