// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the bio::io::format_input_handler<bio::io::bcf> .
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <ranges>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

#include <bio/meta/tag/vtag.hpp>
#include <bio/meta/type_list/traits.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/detail/misc.hpp>
#include <bio/io/detail/range.hpp>
#include <bio/io/detail/to_little_endian.hpp>
#include <bio/io/format/bcf.hpp>
#include <bio/io/format/format_input_handler.hpp>
#include <bio/io/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/io/var/header.hpp>
#include <bio/io/var/misc.hpp>
#include <bio/io/var/reader_options.hpp> //TODO for field_types_raw; move somewhere else

namespace bio::io::detail
{
//!\brief Low-level input iterator that points to the byte-region of a BCF record in the input file.
class bcf_input_iterator
{
private:
    //!\brief Down-cast pointer to the stream-buffer.
    bio::io::detail::stream_buffer_exposer<char> * stream_buf = nullptr;

    //!\brief Place to store lines that overlap buffer boundaries.
    std::basic_string<char> overflow_buffer;
    //!\brief Temporary storage for field delimiter positions.
    std::vector<size_t>     field_end_positions;

    //!\brief The record, consisting of byte-span and position of genotypes-begin.
    std::pair<std::span<std::byte const>, size_t> record;

    //!\brief Whether iterator is at end.
    bool at_end = false;

public:
    //!\brief The bcf format header.
    var::detail::bcf_header header;

    /*!\name Associated types
     * \{
     */
    using difference_type   = ptrdiff_t; //!< Defaults to ptrdiff_t.
    //!\brief The record type.
    using value_type        = std::pair<std::span<std::byte const>, size_t>;
    using reference         = value_type;              //!< Same as record type (passed by value).
    using pointer           = value_type const *;      //!< Has no pointer type.
    using iterator_category = std::input_iterator_tag; //!< Pure input iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bcf_input_iterator() noexcept                                       = default; //!< Defaulted.
    bcf_input_iterator(bcf_input_iterator const &) noexcept             = default; //!< Defaulted.
    bcf_input_iterator(bcf_input_iterator &&) noexcept                  = default; //!< Defaulted.
    bcf_input_iterator & operator=(bcf_input_iterator const &) noexcept = default; //!< Defaulted.
    bcf_input_iterator & operator=(bcf_input_iterator &&) noexcept      = default; //!< Defaulted.
    ~bcf_input_iterator() noexcept                                      = default; //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit bcf_input_iterator(std::basic_streambuf<char> & ibuf, bool const init = true) :
      stream_buf{reinterpret_cast<bio::io::detail::stream_buffer_exposer<char> *>(&ibuf)}
    {
        assert(stream_buf != nullptr);
        stream_buf->underflow(); // ensure the stream buffer has content on construction

        // READ HEADER
        if (stream_buf->egptr() - stream_buf->gptr() < 9)
            throw format_error{"File is too small to be a BCF file."};

        header.magic[0]      = *stream_buf->gptr();
        header.magic[1]      = *(stream_buf->gptr() + 1);
        header.magic[2]      = *(stream_buf->gptr() + 2);
        header.major_version = *(stream_buf->gptr() + 3);
        header.minor_version = *(stream_buf->gptr() + 4);
        header.l_text        = detail::to_little_endian(*(reinterpret_cast<uint32_t *>(stream_buf->gptr() + 5)));
        stream_buf->gbump(9);

        if (header.magic[0] != 'B' || header.magic[1] != 'C' || header.magic[2] != 'F')
            throw format_error{"File does not start with BCF magic header."};

        if (header.major_version != 2 || (header.minor_version != 1 && header.minor_version != 2))
            throw format_error{"Only BCF 2.1 and BCF 2.2 are supported."};

        header.text.resize(header.l_text);
        ptrdiff_t written = 0;
        while (written < header.l_text)
        {
            if (stream_buf->egptr() == stream_buf->gptr())
                throw format_error{"End-of-BCF-Stream encountered before the header could be read completely."};

            ptrdiff_t remaining         = header.l_text - written;
            ptrdiff_t to_be_written_now = std::min<ptrdiff_t>(remaining, stream_buf->egptr() - stream_buf->gptr());

            std::ranges::copy(stream_buf->gptr(),
                              stream_buf->gptr() + to_be_written_now,
                              header.text.begin() + written);
            written += to_be_written_now;

            stream_buf->gbump(to_be_written_now);
            if (written < header.l_text)
                stream_buf->underflow();
        }

        header.text.resize(header.l_text - 1); // remove trailing '\0'

        /* zero records is allowed */
        if (stream_buf->gptr() == stream_buf->egptr())
            at_end = true;

        if (init) // read first record
        {
            operator++();
        }
    }

    //!\brief Construct from a stream.
    explicit bcf_input_iterator(std::basic_istream<char> & istr, bool const init = true) :
      bcf_input_iterator{*istr.rdbuf(), init}
    {}
    //!\}

    //!\brief Resets the state of this iterator to cope with changes in underlying stream.
    void reset(std::basic_streambuf<char> & ibuf)
    {
        stream_buf = reinterpret_cast<bio::io::detail::stream_buffer_exposer<char> *>(&ibuf);
        overflow_buffer.clear();
        field_end_positions.clear();
        at_end = false;
    }

    //!\overload
    void reset(std::basic_istream<char> & istr) { reset(*istr.rdbuf()); }

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advance by one and rebuffer if necessary (vtable lookup iff rebuffering).
    bcf_input_iterator & operator++()
    {
        assert(stream_buf != nullptr);

        if (at_end)
            return *this;

        if (stream_buf->egptr() - stream_buf->gptr() == 0) // possible to be on empty buffer
        {
            stream_buf->underflow();
            if (stream_buf->egptr() - stream_buf->gptr() == 0)
            {
                at_end = true;
                return *this;
            }
        }

        auto & [span, genotype_offset] = record;

        overflow_buffer.clear();

        uint32_t l_shared = 0;
        uint32_t l_indiv  = 0;

        // we have buffer, but it is not large enough to contain l_shared and l_indiv
        if (ptrdiff_t remaining_buffer_size = stream_buf->egptr() - stream_buf->gptr(); remaining_buffer_size < 8)
        {
            overflow_buffer.resize(8);
            std::ranges::copy(stream_buf->gptr(), stream_buf->egptr(), overflow_buffer.begin());
            stream_buf->gbump(remaining_buffer_size);
            stream_buf->underflow();

            ptrdiff_t remaining_data_size = 8 - remaining_buffer_size;

            if (stream_buf->egptr() - stream_buf->gptr() < remaining_data_size)
                throw format_error{"End-of-BCF-Stream encountered before the record data could read."};

            std::ranges::copy(stream_buf->gptr(),
                              stream_buf->gptr() + remaining_data_size,
                              overflow_buffer.begin() + remaining_buffer_size);
            stream_buf->gbump(remaining_data_size);

            l_shared = detail::to_little_endian(*(reinterpret_cast<uint32_t *>(overflow_buffer.data())));
            l_indiv  = detail::to_little_endian(*(reinterpret_cast<uint32_t *>(overflow_buffer.data() + 4)));
            overflow_buffer.clear();
        }
        else // read l_shared and l_indiv directly
        {
            l_shared = detail::to_little_endian(*(reinterpret_cast<uint32_t *>(stream_buf->gptr())));
            l_indiv  = detail::to_little_endian(*(reinterpret_cast<uint32_t *>(stream_buf->gptr() + 4)));
            stream_buf->gbump(8); // skip l_shared and l_indiv
        }

        ptrdiff_t record_size = l_shared + l_indiv;
        genotype_offset       = l_shared;

        if (record_size == 0)
            throw format_error{"Record of size 0 found when reading BCF stream."};

        ptrdiff_t remaining_buffer_size = stream_buf->egptr() - stream_buf->gptr();

        if (record_size <= remaining_buffer_size)
        {
            span =
              std::span<std::byte const>{reinterpret_cast<std::byte const *>(stream_buf->gptr()), size_t(record_size)};
            stream_buf->gbump(record_size);
        }
        else // record is not completely in buffer
        {
            overflow_buffer.resize(record_size);
            ptrdiff_t written = 0;

            while (written < record_size) // record might span multiple buffer-lengths (unlikely)
            {
                if (stream_buf->egptr() == stream_buf->gptr())
                    throw format_error{"End-of-BCF-Stream encountered before the record ended."};

                ptrdiff_t remaining         = record_size - written;
                ptrdiff_t to_be_written_now = std::min<ptrdiff_t>(remaining, stream_buf->egptr() - stream_buf->gptr());

                std::ranges::copy(stream_buf->gptr(),
                                  stream_buf->gptr() + to_be_written_now,
                                  overflow_buffer.begin() + written);
                written += to_be_written_now;

                stream_buf->gbump(to_be_written_now);
                if (written < record_size)
                    stream_buf->underflow();
            }

            span = std::span<std::byte const>{reinterpret_cast<std::byte const *>(overflow_buffer.data()),
                                              overflow_buffer.size()};
        }

        return *this;
    }

    //!\overload
    void operator++(int) { ++(*this); }
    //!\}

    /*!\name Dereference operators
     * \{
     */
    //!\brief Read current value from buffer (no vtable lookup, safe even at end).
    reference operator*() const { return record; }

    //!\brief Read current value from buffer (no vtable lookup, safe even at end).
    pointer operator->() const { return &record; }
    //!\}

    /*!\name Comparison operators
     * \brief We define comparison only against the sentinel.
     * \{
     */
    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(bcf_input_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return lhs.at_end;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(bcf_input_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::default_sentinel);
    }

    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(std::default_sentinel_t const &, bcf_input_iterator const & rhs) noexcept
    {
        return rhs == std::default_sentinel;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(std::default_sentinel_t const &, bcf_input_iterator const & rhs) noexcept
    {
        return !(rhs == std::default_sentinel);
    }
    //!\}
};

} // namespace bio::io::detail

namespace bio::io
{
/*!\brief Format input handler for the BCF format (bio::io::bcf).
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
 * | Member          | Type    | Default | Description                                                      |
 * |-----------------|---------|---------|------------------------------------------------------------------|
 * |`print_warnings` |`bool`   | `false` | Whether to print non-critical warngings to std::cerr             |
 *
 * ### Performance
 *
 * TODO after genotype redesign
 */
template <>
class format_input_handler<bcf> :
  public format_input_handler_base<format_input_handler<bcf>>,
  public var::format_handler_mixin
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The type of the CRTP base class.
    using base_t = format_input_handler_base<format_input_handler<bcf>>;
    using base_t::parse_field;
    using base_t::parse_field_aux;
    using base_t::stream;

    //!\brief Befriend the base class to enable CRTP.
    friend base_t;
    //!\}

    //!\brief Print an error message with current record number in diagnostic.
    [[noreturn]] void error(auto const &... messages) const
    {
        throw format_error{"[BCF format error in record ", record_no, "] ", messages...};
    }

    /*!\name Options
     * \{
     */
    //!\brief Whether to print warnings or not.
    bool print_warnings = true;
    //!\}

    /*!\name Binary decoders
     * \brief These interpret binary data.
     * \{
     */
    //!\brief Decode the type descripter byte into a var::detail::bcf_type_descriptor and a number.
    std::pair<var::detail::bcf_type_descriptor, uint8_t> decode_type_descriptor_byte(std::byte const b) const
    {
        uint8_t desc = static_cast<uint8_t>(b) & 0b00001111;
        uint8_t size = static_cast<uint8_t>(b) >> 4;

        switch (static_cast<uint8_t>(desc))
        {
            case 4:
            case 6:
            case 8:
            case 12:
            case 14:
                error("Cannot decode BCF type descriptor byte. Value: ", desc, " is unknown.");
            default:
                break;
        }

        return {var::detail::bcf_type_descriptor{desc}, size};
    }

    //!\brief Decodes into a number and advances the cache_ptr.
    int32_t decode_integral(std::byte const *& cache_ptr) const
    {
        assert(cache_ptr < &*file_it->first.end());

        auto [desc, size] = decode_type_descriptor_byte(*cache_ptr);
        ++cache_ptr;

        if (size != 1)
            error("Trying to decode single integral value but found ", size, " values.");

        int32_t number = -1;

        switch (desc)
        {
            case var::detail::bcf_type_descriptor::int8:
                number = *reinterpret_cast<int8_t const *>(cache_ptr);
                cache_ptr += 1;
                break;
            case var::detail::bcf_type_descriptor::int16:
                number = *reinterpret_cast<int16_t const *>(cache_ptr);
                cache_ptr += 2;
                break;
            case var::detail::bcf_type_descriptor::int32:
                number = *reinterpret_cast<int32_t const *>(cache_ptr);
                cache_ptr += 4;
                break;
            default:
                error("Trying to decode an integral value but found something else.");
                break;
        }

        assert(cache_ptr <= &*file_it->first.end());
        return number;
    }

    //!\brief Decodes into a string and advances the cache_ptr.
    std::string_view decode_string(std::byte const *& cache_ptr) const
    {
        assert(cache_ptr < &*file_it->first.end());

        auto [desc, size] = decode_type_descriptor_byte(*cache_ptr);
        ++cache_ptr;

        if (desc != var::detail::bcf_type_descriptor::char8)
            error("Trying to decode a string but found something else.");

        size_t real_size = size < 15 ? size : decode_integral(cache_ptr);

        std::string_view ret{reinterpret_cast<char const *>(cache_ptr), real_size};

        cache_ptr += real_size;

        assert(cache_ptr <= &*file_it->first.end());

        return ret;
    }

    //!\brief Return current field as bytes (no actual decoding) and advance cache_ptr.
    std::span<std::byte const> decode_any_field(std::byte const *& cache_ptr) const
    {
        assert(cache_ptr < &*file_it->first.end());

        std::byte const * b = cache_ptr;

        auto [desc, size] = decode_type_descriptor_byte(*cache_ptr);
        ++cache_ptr;

        size_t real_size = size < 15 ? size : decode_integral(cache_ptr);

        // TODO shouldn't we take into account the element size?
        cache_ptr += real_size;

        std::span<std::byte const> ret{b, size_t(cache_ptr - b)};

        assert(cache_ptr <= &*file_it->first.end());

        return ret;
    }

    // IMPLEMENTATION NOTE: the following is the only decoder that is used in parsed record reading.
    //!\brief Decodes any range of integral and stores in any range of integral (no check of sizes).
    template <ranges::back_insertable out_t>
        requires detail::int_range<out_t>
    void decode_numbers_into(std::span<std::byte const> in, out_t & out)
    {
        std::byte const *                  cache_ptr     = in.data();
        [[maybe_unused]] std::byte const * cache_end_ptr = cache_ptr + in.size();

        auto [desc, size] = decode_type_descriptor_byte(*cache_ptr);
        ++cache_ptr;

        size_t real_size = size < 15 ? size : decode_integral(cache_ptr);

        switch (desc)
        {
            case var::detail::bcf_type_descriptor::missing:
                return;
            case var::detail::bcf_type_descriptor::int8:
                {
                    assert((cache_end_ptr - cache_ptr) == (ptrdiff_t)real_size);
                    std::span<int8_t const> data{reinterpret_cast<int8_t const *>(cache_ptr), real_size};
                    detail::sized_range_copy(data, out);
                    break;
                }
            case var::detail::bcf_type_descriptor::int16:
                {
                    assert((cache_end_ptr - cache_ptr) == (ptrdiff_t)real_size * 2);
                    std::span<int16_t const> data{reinterpret_cast<int16_t const *>(cache_ptr), real_size};
                    detail::sized_range_copy(data, out);
                    break;
                }
            case var::detail::bcf_type_descriptor::int32:
                {
                    assert((cache_end_ptr - cache_ptr) == (ptrdiff_t)real_size * 4);
                    std::span<int32_t const> data{reinterpret_cast<int32_t const *>(cache_ptr), real_size};
                    detail::sized_range_copy(data, out);
                    break;
                }
            default:
                error("Expected numbers but got something else.");
                break;
        }
    }
    //!\}

    /*!\name Raw record handling
     * \{
     */
    //!\brief The fields that this format supports [the base class accesses this type].
    using format_fields = std::remove_cvref_t<decltype(var::detail::field_ids)>;
    //!\brief Type of the raw record.
    using raw_record_type =
      io::detail::tuple_record<format_fields,
                               meta::list_traits::concat<meta::list_traits::repeat<9, std::span<std::byte const>>,
                                                         meta::type_list<var::record_private_data>>>;

    //!\brief The raw record.
    raw_record_type raw_record;

    //!\brief This is cached during raw-record reading already.
    std::string_view                     id_cache;
    //!\brief This is cached during raw-record reading already.
    std::string_view                     ref_cache;
    //!\brief This is cached during raw-record reading already.
    std::vector<std::string_view>        alts_cache;
    //!\brief This is cached during raw-record reading already.
    var::detail::bcf_record_core const * record_core = nullptr;
    //!\brief Storage for GT values which are not encoded as strings (but pretend to be).
    std::vector<std::string>             gt_cache; // TODO concatenated_sequences
    //!\brief Temporary storage for string/number conversions.
    std::string                          number_cache;

    //!\brief The header.
    var::header                header;
    //!\brief Lowlevel stream iterator.
    detail::bcf_input_iterator file_it;
    //!\brief Current record number in file.
    size_t                     record_no = 0;

    //!\brief Read the raw record [the base class invokes this function].
    void read_raw_record()
    {
        using raw_f_t = std::span<std::byte const>;
        ++record_no;
        ++file_it;

        auto [record_data, genotype_offset] = *file_it;

        alts_cache.clear();

        record_core = reinterpret_cast<var::detail::bcf_record_core const *>(record_data.data());
        gt_cache.resize(record_core->n_sample);

        std::byte const * cache_ptr           = nullptr;
        std::byte const * cache_ptr_bak       = nullptr;
        // offset, count
        get<detail::field::chrom>(raw_record) = record_data.subspan(0, 4);
        get<detail::field::pos>(raw_record)   = record_data.subspan(4, 4);
        get<detail::field::qual>(raw_record)  = record_data.subspan(12, 4);

        cache_ptr = record_data.data() + 24; // == sizeof(record_core)

        cache_ptr_bak                      = cache_ptr;
        id_cache                           = decode_string(cache_ptr); // store string_view already
        get<detail::field::id>(raw_record) = raw_f_t{cache_ptr_bak, size_t(cache_ptr - cache_ptr_bak)};

        cache_ptr_bak                       = cache_ptr;
        ref_cache                           = decode_string(cache_ptr); // store string_view already
        get<detail::field::ref>(raw_record) = raw_f_t{cache_ptr_bak, size_t(cache_ptr - cache_ptr_bak)};

        cache_ptr_bak = cache_ptr;
        for (size_t i = 1; i < record_core->n_allele; ++i)
            alts_cache.push_back(decode_string(cache_ptr));
        get<detail::field::alt>(raw_record) = raw_f_t{cache_ptr_bak, size_t(cache_ptr - cache_ptr_bak)};

        get<detail::field::filter>(raw_record) = decode_any_field(cache_ptr);

        size_t distance_to_genotypes         = (record_data.data() + genotype_offset) - cache_ptr;
        get<detail::field::info>(raw_record) = raw_f_t{cache_ptr, distance_to_genotypes};

        get<detail::field::genotypes>(raw_record) = record_data.subspan(genotype_offset /*, end*/);
    }
    //!\}

    /*!\name Dynamic type initialisation and parsing
     * \{
     */
    //!\brief Set to single value.
    template <var::type_enum id_, typename elem_t>
    static inline void element_value_type_init_single(std::byte const *& cache_ptr, auto & output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();
        output_                  = *reinterpret_cast<elem_t const *>(cache_ptr);
        cache_ptr += sizeof(elem_t);
    }

    //!\brief Set to string.
    template <var::type_enum id_>
    static inline void element_value_type_init_string(size_t const size, std::byte const *& cache_ptr, auto & output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();

        std::string_view tmp{reinterpret_cast<char const *>(cache_ptr), size};
        detail::string_copy(tmp, output_);

        cache_ptr += size;
    }

    //!\brief Set to vector.
    template <var::type_enum id_, typename elem_t>
    static inline void element_value_type_init_vector(size_t const size, std::byte const *& cache_ptr, auto & output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();

        std::span<elem_t const> tmp{reinterpret_cast<elem_t const *>(cache_ptr), size};
        detail::sized_range_copy(tmp, output_);

        cache_ptr += size * sizeof(elem_t);
    }

    //!\brief Set to vector-of-string (genotypes and info have different implementation here).
    template <var::type_enum id_>
    static inline void info_element_value_type_init_vector_of_string(size_t const       size,
                                                                     std::byte const *& cache_ptr,
                                                                     auto &             output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();

        std::string_view tmp{reinterpret_cast<char const *>(cache_ptr), size};
        for (std::string_view const s : tmp | detail::eager_split(','))
        {
            output_.emplace_back();
            detail::string_copy(s, output_.back());
        }

        cache_ptr += size;
    }

    //!\brief Set to vector-of-string (genotypes and info have different implementation here).
    template <var::type_enum id_>
    static inline void genotype_variant_init_vector_of_string(size_t const       outer_size,
                                                              size_t const       inner_size,
                                                              std::byte const *& cache_ptr,
                                                              auto &             output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();

        output_.reserve(outer_size);
        // TODO this is not yet concatenated_sequences
        //         output_.concat_reserve(outer_size * inner_size);

        for (size_t i = 0; i < outer_size; ++i)
        {
            std::string_view tmp{reinterpret_cast<char const *>(cache_ptr) + i * inner_size, inner_size};

            // vectors can be smaller by being padded with end-of-vector values
            size_t s = tmp.size();
            while (s > 0 && tmp[s - 1] == var::detail::end_of_vector<char>)
                --s;
            tmp = tmp.substr(0, s);

            output_.emplace_back();
            detail::string_copy(tmp, output_.back());
        }
        cache_ptr += outer_size * inner_size * sizeof(char);
    }

    //!\brief Set to vector-of-vector.
    template <var::type_enum id_, typename elem_t>
    static inline void element_value_type_init_vector_of_vector(size_t const       outer_size,
                                                                size_t const       inner_size,
                                                                std::byte const *& cache_ptr,
                                                                auto &             output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();

        output_.reserve(outer_size);
        // in cases where vectors are padded, we might be over-allocating here, but it's better than not doing it:
        output_.concat_reserve(outer_size * inner_size);

        for (size_t i = 0; i < outer_size; ++i)
        {
            std::span<elem_t const> tmp{reinterpret_cast<elem_t const *>(cache_ptr) + i * inner_size, inner_size};

            // vectors can be smaller by being padded with end-of-vector values
            size_t s = tmp.size();
            while (s > 0 && tmp[s - 1] == var::detail::end_of_vector<elem_t>)
                --s;
            tmp = tmp.subspan(0, s);

            output_.push_back(tmp);
        }
        cache_ptr += outer_size * inner_size * sizeof(elem_t);
    }

    template <var::detail::is_info_variant dyn_t>
    void parse_element_value_type(var::type_enum const                   id_from_header,
                                  var::detail::bcf_type_descriptor const desc,
                                  size_t const                           size,
                                  std::byte const *&                     cache_ptr,
                                  dyn_t &                                output); // implementation below class

    template <var::detail::is_genotype_variant dyn_t>
    void parse_element_value_type(var::type_enum const                   id_from_header,
                                  var::detail::bcf_type_descriptor const desc,
                                  size_t const                           outer_size,
                                  size_t const                           inner_size,
                                  std::byte const *&                     cache_ptr,
                                  dyn_t &                                output); // implementation below class
    //!\}

    /*!\name Parsed record handling
     * \{
     */
    //!\brief Reading of CHROM field.
    void parse_field(meta::vtag_t<detail::field::chrom> const & /**/, std::integral auto & parsed_field)
    {
        parsed_field = record_core->chrom;
    }

    //!\overload
    void parse_field(meta::vtag_t<detail::field::chrom> const & /**/, auto & parsed_field)
    {
        parse_field_aux(header.contig_idx_to_string_map().at(record_core->chrom), parsed_field);
    }

    //!\brief Reading of POS field.
    void parse_field(meta::vtag_t<detail::field::pos> const & /**/, std::integral auto & parsed_field)
    {
        parsed_field = record_core->pos + 1; // one-based positions
    }

    //!\brief Reading of ID field.
    void parse_field(meta::vtag_t<detail::field::id> const & /**/, auto & parsed_field)
    {
        if (std::ranges::empty(id_cache))
            parse_field_aux(std::string_view{"."}, parsed_field);
        else
            parse_field_aux(id_cache, parsed_field);
    }

    //!\brief Reading of REF field.
    void parse_field(meta::vtag_t<detail::field::ref> const & /**/, auto & parsed_field)
    {
        parse_field_aux(ref_cache, parsed_field);
    }

    //!\brief Reading of ALT field.
    void parse_field(meta::vtag_t<detail::field::alt> const & /**/, std::vector<std::string_view> & parsed_field)
    {
        parsed_field = alts_cache;
    }

    //!\overload
    template <ranges::back_insertable parsed_field_t>
        requires(std::ranges::range<std::ranges::range_reference_t<parsed_field_t>>)
    void parse_field(meta::vtag_t<detail::field::alt> const & /**/, parsed_field_t & parsed_field)
    {
        for (std::string_view const alt : alts_cache)
        {
            std::ranges::range_value_t<parsed_field_t> out;

            parse_field_aux(alt, out);

            parsed_field.push_back(std::move(out));
        }
    }

    //!\brief Reading of QUAL field.
    void parse_field(meta::vtag_t<detail::field::qual> const & /**/, meta::arithmetic auto & parsed_field)
    {
        parsed_field = record_core->qual;
    }

    //!\brief Reading of FILTER field.
    template <ranges::back_insertable parsed_field_t>
        requires detail::int_range<parsed_field_t>
    void parse_field(meta::vtag_t<detail::field::filter> const & /**/, parsed_field_t & parsed_field)
    {
        decode_numbers_into(get<detail::field::filter>(raw_record), parsed_field);
    }

    //!\overload
    template <ranges::back_insertable parsed_field_t>
        requires detail::out_string<std::ranges::range_reference_t<parsed_field_t>>
    void parse_field(meta::vtag_t<detail::field::filter> const & /**/, parsed_field_t & parsed_field)
    {
        std::vector<int32_t> tmp; // ATTENTION this allocates, TODO change

        decode_numbers_into(get<detail::field::filter>(raw_record), tmp);

        for (int32_t const idx : tmp)
        {
            std::ranges::range_value_t<parsed_field_t> out;
            parse_field_aux(header.idx_to_string_map().at(idx), out);
            parsed_field.push_back(std::move(out));
        }
    }

    //!\brief Reading of the INFO field.
    template <ranges::back_insertable parsed_field_t>
        requires var::detail::info_element_reader_concept<std::ranges::range_value_t<parsed_field_t>>
    void parse_field(meta::vtag_t<detail::field::info> const & /**/, parsed_field_t & parsed_field)
    {
        std::span<std::byte const> raw_field = get<detail::field::info>(raw_record);
        std::byte const *          cache_ptr = raw_field.data();

        for (size_t i = 0; i < record_core->n_info; ++i)
        {
            std::ranges::range_value_t<parsed_field_t> new_element;

            auto & [id, variant] = new_element;

            int32_t idx = decode_integral(cache_ptr);

            assert(cache_ptr < raw_field.data() + raw_field.size());

            auto [desc, size] = decode_type_descriptor_byte(*cache_ptr);
            ++cache_ptr;

            size_t real_size = size < 15 ? size : decode_integral(cache_ptr);

            auto info_it = header.infos.find(header.idx_to_string_map().at(idx));
            assert(info_it != header.infos.end());
            auto && [id_string, info] = *info_it;

            parse_element_value_type(info.type_id, desc, real_size, cache_ptr, variant);

            if constexpr (detail::out_string<decltype(id)>)
                detail::string_copy(id_string, id);
            else
                id = idx;

            parsed_field.push_back(std::move(new_element));
        }

        assert(cache_ptr == raw_field.data() + raw_field.size());
    }

    //!\brief Auxilliary function for reading GT which isn't actually encoded as a string.
    void parse_gt_field(auto && in, std::string & number_cache, std::string & output) const
    {
        using int_t = std::ranges::range_value_t<decltype(in)>;
        bool first  = true;

        for (int_t i : in)
        {
            if (i == var::detail::end_of_vector<int_t>)
                return;

            bool phased = i % 2;
            i >>= 1;

            if (first)
                first = false;
            else
                output.push_back(phased ? '|' : '/');

            if (i == 0)
            {
                output.push_back('.');
            }
            else
            {
                --i;
                number_cache.reserve(12); // large enough to hold all int32_t representations and inside SSO

                if (auto [ptr, errc] = std::to_chars(number_cache.data(), number_cache.data() + 12, i);
                    errc == std::errc{})
                {
                    std::ranges::copy(number_cache.data(), ptr, std::back_insert_iterator{output});
                }
                else
                {
                    error("Error parsing binary representation of number in GT field.");
                }
            }
        }
    }

    //!\brief Implementation for parsing into genotypes.
    void parse_genotypes_impl(auto & parsed_field)
    {
        std::span<std::byte const> raw_field = get<detail::field::genotypes>(raw_record);

        std::byte const * cache_ptr = raw_field.data();

        std::string number_cache;

        for (size_t i = 0; i < static_cast<size_t>(record_core->n_fmt); ++i)
        {
            std::ranges::range_value_t<decltype(parsed_field)> new_element;

            auto & [id, parsed_variant] = new_element;

            int32_t fmt_key = decode_integral(cache_ptr);
            auto    fmt_it  = header.formats.find(header.idx_to_string_map().at(fmt_key));
            assert(fmt_it != header.formats.end());
            auto && [id_string, format] = *fmt_it;

            if constexpr (detail::out_string<decltype(id)>)
                detail::string_copy(id_string, id);
            else
                id = fmt_key;

            auto [fmt_type, fmt_size] = decode_type_descriptor_byte(*cache_ptr);
            ++cache_ptr;

            // TODO fmt_size might overflow because it is only uint8_t
            fmt_size = fmt_size < 15 ? fmt_size : decode_integral(cache_ptr);

            if (id_string == "GT") // this needs custom decoding, it is not a string
            {
                /* we explicitly parse into integer format for now: */
                parse_element_value_type(var::type_enum::vector_of_int32,
                                         fmt_type,
                                         record_core->n_sample,
                                         fmt_size,
                                         cache_ptr,
                                         parsed_variant);

                /* we transform number to string and store in caches */
                std::visit(
                  [&]<typename rng_t>(rng_t && rng)
                  {
                      if constexpr (std::ranges::range<std::ranges::range_value_t<rng_t>> &&
                                    detail::int_range<std::ranges::range_value_t<rng_t>>)
                      {
                          if (std::ranges::size(rng) != record_core->n_sample)
                              error("Expected exactly one GT string per sample.");

                          for (size_t sample = 0; sample < record_core->n_sample; ++sample)
                          {
                              gt_cache[sample].clear();
                              number_cache.clear();
                              parse_gt_field(rng[sample], number_cache, gt_cache[sample]);
                          }
                      }
                      else
                      {
                          // THIS NEVER HAPPENS BUT WE SAVE INSTANTIATIONS with if-constexpr
                          error("Unreachable state reached at code line: ", __LINE__);
                      }
                  },
                  parsed_variant);

                /* now reset the variant to string-state and create copies/views */
                constexpr size_t string_id = static_cast<size_t>(var::type_enum::string);
                auto &           strings   = parsed_variant.template emplace<string_id>();
                strings.resize(record_core->n_sample);

                for (size_t sample = 0; sample < record_core->n_sample; ++sample)
                {
                    if constexpr (std::same_as<std::string &, decltype(strings[sample])>)
                        std::swap(strings[sample], gt_cache[sample]); // TODO better to copy one-way because of SSO?
                    else
                        strings[sample] = gt_cache[sample]; // string to string_view
                }
            }
            else
            {
                parse_element_value_type(format.type_id,
                                         fmt_type,
                                         record_core->n_sample,
                                         fmt_size,
                                         cache_ptr,
                                         parsed_variant);
            }

            parsed_field.push_back(std::move(new_element));
        }

        assert(cache_ptr == raw_field.data() + raw_field.size());
    }

    //!\brief Reading of the GENOTYPES field.
    template <ranges::back_insertable field_t>
        requires var::detail::genotype_reader_concept<std::ranges::range_value_t<field_t>>
    void parse_field(meta::vtag_t<detail::field::genotypes> const & /**/, field_t & parsed_field)
    {
        parse_genotypes_impl(parsed_field);
    }

    //!\brief Overload for parsing the private data.
    void parse_field(meta::vtag_t<detail::field::_private> const & /**/, var::record_private_data & parsed_field)
    {
        parsed_field.header_ptr  = &header;
        parsed_field.raw_record  = &raw_record;
        parsed_field.record_core = record_core; // already a pointer
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
     * The options argument is typically bio::io::var::reader_options, but any object with a subset of similarly
     * named members is also accepted. See bio::io::format_input_handler<bio::io::bcf> for the supported options and
     * defaults.
     */
    template <typename options_t>
    format_input_handler(std::istream & str, options_t const & options) : base_t{str}, file_it{str, false /*no_init!*/}

    {
        // extract options
        if constexpr (requires { (bool)options.print_warnings; })
        {
            print_warnings = options.print_warnings;
        }

        std::string text_tmp;
        std::swap(text_tmp, file_it.header.text);
        header = var::header{std::move(text_tmp)};

        /* checks on header */
        if (header.formats.contains("GT") && header.formats["GT"].type_id != var::type_enum::string)
        {
            error("The \"GT\" field must always be encoded as a string.");
        }
        // TODO more checks
    }

    //!\brief Construct with only an input stream.
    format_input_handler(std::istream & str) : format_input_handler{str, int{}} {}
    //!\}

    //!\brief Return a reference to the header contained in the input handler.
    var::header const & get_header() const { return header; }

    //!\brief This resets the stream iterator after region-seek.
    void reset_stream() { file_it.reset(*stream); }
};

/*!\brief Parse a "dynamically typed" field out of a BCF stream and store the content in a variant.
 * \tparam dyn_t             Type of the variant; specialisation of bio::io::var::info_element_value_type.
 * \param[in] id_from_header A value of bio::io::var::type_enum that denotes the expected type.
 * \param[in] desc           A value of bio::io::detail::bcf_type_descriptor that notes the detected type.
 * \param[in] size           The number of values belonging to this field.
 * \param[in,out] cache_ptr  Pointer into the BCF stream; will be updated to point past the end of read data.
 * \param[out] output        The variant to hold the parsed value.
 */
template <var::detail::is_info_variant dyn_t>
inline void format_input_handler<bcf>::parse_element_value_type(var::type_enum const                   id_from_header,
                                                                var::detail::bcf_type_descriptor const desc,
                                                                size_t const                           size,
                                                                std::byte const *&                     cache_ptr,
                                                                dyn_t &                                output)
{
    // TODO DRY out the boilerplate error messages
    if (static_cast<size_t>(id_from_header) < static_cast<size_t>(var::type_enum::string) && size != 1)
        error("BCF data field expected exactly one element, got:", size);

    switch (id_from_header)
    {
        case var::type_enum::char8:
            {
                if (desc != var::detail::bcf_type_descriptor::char8)
                    error("Attempting to create char but the byte descriptor does not indicate char type.");

                element_value_type_init_single<var::type_enum::char8, char>(cache_ptr, output);
                return;
            }
        case var::type_enum::int8:
        case var::type_enum::int16:
        case var::type_enum::int32:
            {
                switch (desc)
                {
                    case var::detail::bcf_type_descriptor::int8:
                        element_value_type_init_single<var::type_enum::int8, int8_t>(cache_ptr, output);
                        break;
                    case var::detail::bcf_type_descriptor::int16:
                        element_value_type_init_single<var::type_enum::int16, int16_t>(cache_ptr, output);
                        break;
                    case var::detail::bcf_type_descriptor::int32:
                        element_value_type_init_single<var::type_enum::int32, int32_t>(cache_ptr, output);
                        break;
                    default:
                        error("Attempting to create int but the byte descriptor does not indicate int type.");
                }
                return;
            }
        case var::type_enum::float32:
            {
                if (desc != var::detail::bcf_type_descriptor::float32)
                    error("Attempting to create float but the byte descriptor does not indicate float type.");

                element_value_type_init_single<var::type_enum::float32, float>(cache_ptr, output);
                return;
            }
        case var::type_enum::string:
            {
                if (desc != var::detail::bcf_type_descriptor::char8)
                    error("Attempting to creates string but the byte descriptor does not indicate string type.");

                element_value_type_init_string<var::type_enum::string>(size, cache_ptr, output);
                return;
            }
        case var::type_enum::vector_of_int8:
        case var::type_enum::vector_of_int16:
        case var::type_enum::vector_of_int32:
            {
                switch (desc)
                {
                    case var::detail::bcf_type_descriptor::int8:
                        {
                            element_value_type_init_vector<var::type_enum::vector_of_int8, int8_t>(size,
                                                                                                   cache_ptr,
                                                                                                   output);
                            break;
                        }
                    case var::detail::bcf_type_descriptor::int16:
                        {
                            element_value_type_init_vector<var::type_enum::vector_of_int16, int16_t>(size,
                                                                                                     cache_ptr,
                                                                                                     output);
                            break;
                        }
                    case var::detail::bcf_type_descriptor::int32:
                        {
                            element_value_type_init_vector<var::type_enum::vector_of_int32, int32_t>(size,
                                                                                                     cache_ptr,
                                                                                                     output);
                            break;
                        }
                    default:
                        error("Attempting to create vector of int but the byte descriptor does not indicate int type.");
                }
                return;
            }
        case var::type_enum::vector_of_float32:
            {
                if (desc != var::detail::bcf_type_descriptor::float32)
                    error("Attempting to create vector of float but the byte descriptor does not indicate float type.");

                element_value_type_init_vector<var::type_enum::vector_of_float32, float>(size, cache_ptr, output);
                return;
            }
        case var::type_enum::vector_of_string:
            {
                if (desc != var::detail::bcf_type_descriptor::char8)
                    error(
                      "Attempting to create vector of string but the byte descriptor does not indicate char alphabet");

                info_element_value_type_init_vector_of_string<var::type_enum::vector_of_string>(size,
                                                                                                cache_ptr,
                                                                                                output);
                return;
            }
        case var::type_enum::flag:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::flag);
                output.template emplace<id>(true);

                cache_ptr += size; // This should be 0, but is allowed to be something else
                return;
            }
    }
}

/*!\brief Parse a "dynamically typed" field out of a BCF stream and store the content in a vector-variant.
 * \tparam dyn_t             Type of the variant; specialisation of bio::io::var::info_element_value_type.
 * \param[in] id_from_header A value of bio::io::var::type_enum that denotes the expected type.
 * \param[in] desc           A value of bio::io::detail::bcf_type_descriptor that notes the detected type.
 * \param[in] outer_size     The number of values belonging to this field.
 * \param[in] inner_size     The number of values per inner vector in case of vector-of-vector.
 * \param[in,out] cache_ptr  Pointer into the BCF stream; will be updated to point past the end of read data.
 * \param[out] output        The variant to hold the parsed value.
 */
template <var::detail::is_genotype_variant dyn_t>
inline void format_input_handler<bcf>::parse_element_value_type(var::type_enum const                   id_from_header,
                                                                var::detail::bcf_type_descriptor const desc,
                                                                size_t const                           outer_size,
                                                                size_t const                           inner_size,
                                                                std::byte const *&                     cache_ptr,
                                                                dyn_t &                                output)
{
    // TODO DRY out the boilerplate error messages
    if (static_cast<size_t>(id_from_header) < static_cast<size_t>(var::type_enum::string) && inner_size != 1)
        error("BCF data field expected exactly one element, got:", inner_size);

    switch (id_from_header)
    {
        case var::type_enum::char8:
            {
                if (desc != var::detail::bcf_type_descriptor::char8)
                    error("Attempting to create char but the byte descriptor does not indicate char type.");

                element_value_type_init_vector<var::type_enum::char8, char>(outer_size, cache_ptr, output);
                return;
            }
        case var::type_enum::int8:
        case var::type_enum::int16:
        case var::type_enum::int32:
            {
                switch (desc)
                {
                    case var::detail::bcf_type_descriptor::int8:
                        {
                            element_value_type_init_vector<var::type_enum::int8, int8_t>(outer_size, cache_ptr, output);
                            break;
                        }
                    case var::detail::bcf_type_descriptor::int16:
                        {
                            element_value_type_init_vector<var::type_enum::int16, int16_t>(outer_size,
                                                                                           cache_ptr,
                                                                                           output);
                            break;
                        }
                    case var::detail::bcf_type_descriptor::int32:
                        {
                            element_value_type_init_vector<var::type_enum::int32, int32_t>(outer_size,
                                                                                           cache_ptr,
                                                                                           output);
                            break;
                        }
                    default:
                        error("Attempting to create int but the byte descriptor does not indicate int type.");
                }
                return;
            }
        case var::type_enum::float32:
            {
                if (desc == var::detail::bcf_type_descriptor::float32)
                    element_value_type_init_vector<var::type_enum::float32, float>(outer_size, cache_ptr, output);
                else
                    error("Attempting to create float but the byte descriptor does not indicate float type.");
                return;
            }
        case var::type_enum::string:
            {
                if (desc != var::detail::bcf_type_descriptor::char8)
                    error("Attempting to creates string but the byte descriptor does not indicate string type.");

                genotype_variant_init_vector_of_string<var::type_enum::string>(outer_size,
                                                                               inner_size,
                                                                               cache_ptr,
                                                                               output);
                return;
            }
        case var::type_enum::vector_of_int8:
        case var::type_enum::vector_of_int16:
        case var::type_enum::vector_of_int32:
            {
                switch (desc)
                {
                    case var::detail::bcf_type_descriptor::int8:
                        {
                            element_value_type_init_vector_of_vector<var::type_enum::vector_of_int8, int8_t>(outer_size,
                                                                                                             inner_size,
                                                                                                             cache_ptr,
                                                                                                             output);
                            break;
                        }
                    case var::detail::bcf_type_descriptor::int16:
                        {
                            element_value_type_init_vector_of_vector<var::type_enum::vector_of_int16, int16_t>(
                              outer_size,
                              inner_size,
                              cache_ptr,
                              output);
                            break;
                        }
                    case var::detail::bcf_type_descriptor::int32:
                        {
                            element_value_type_init_vector_of_vector<var::type_enum::vector_of_int32, int32_t>(
                              outer_size,
                              inner_size,
                              cache_ptr,
                              output);
                            break;
                        }
                    default:
                        error("Attempting to create vector of int but the byte descriptor does not indicate int type.");
                }
                return;
            }
        case var::type_enum::vector_of_float32:
            {
                if (desc == var::detail::bcf_type_descriptor::float32)
                {
                    element_value_type_init_vector_of_vector<var::type_enum::vector_of_float32, float>(outer_size,
                                                                                                       inner_size,
                                                                                                       cache_ptr,
                                                                                                       output);
                    break;
                }
                else
                {
                    error("Attempting to create vector of float but the byte descriptor does not indicate float type.");
                }
                return;
            }
        case var::type_enum::vector_of_string:
            {
                if (desc != var::detail::bcf_type_descriptor::char8)
                    error(
                      "Attempting to create vector of string but the byte descriptor does not indicate char alphabet");

                // TODO this definitely needs a test
                constexpr size_t id      = static_cast<size_t>(var::type_enum::vector_of_string);
                auto &           output_ = output.template emplace<id>();
                std::string_view tmp{reinterpret_cast<char const *>(cache_ptr), outer_size * inner_size};
                for (size_t sample = 0; sample < outer_size; ++sample)
                {
                    output_.emplace_back();

                    std::string_view tmp_inner = tmp.substr(sample * inner_size, inner_size);

                    // string might be padded with end_of_vector values
                    size_t s = tmp_inner.size();
                    while (s > 0 && tmp_inner[s - 1] == var::detail::end_of_vector<char>)
                        --s;
                    tmp_inner = tmp_inner.substr(0, s);

                    for (std::string_view const s : tmp_inner | detail::eager_split(','))
                    {
                        output_.back().push_back({});
                        detail::string_copy(s, output_.back().back());
                    }
                }
                cache_ptr += outer_size * inner_size;
                return;
            }
        case var::type_enum::flag:
            {
                error("bio::io::var::genotype_variant cannot be initialised to flag state.");
                return;
            }
    }
}

} // namespace bio::io
