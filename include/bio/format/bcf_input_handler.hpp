// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the seqan3::format_input_handler<seqan3::format_bcf> .
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

// #include <seqan3/alphabet/views/char_strictly_to.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/debug_stream.hpp> //TODO evaluate if there is a better solution
#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/detail/to_little_endian.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/views/to.hpp>

#include <bio/detail/misc.hpp>
#include <bio/detail/range.hpp>
#include <bio/format/bcf.hpp>
#include <bio/format/format_input_handler.hpp>
#include <bio/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/var_io/dynamic_type.hpp>
#include <bio/var_io/header.hpp>
#include <bio/var_io/misc.hpp>
#include <bio/var_io/reader_options.hpp> //TODO for field_types_raw; move somewhere else

namespace bio::detail
{
//!\brief Low-level input iterator that points to the byte-region of a BCF record in the input file.
class bcf_input_iterator
{
private:
    //!\brief Down-cast pointer to the stream-buffer.
    bio::detail::stream_buffer_exposer<char> * stream_buf = nullptr;

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
    bcf_header header;

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
    bcf_input_iterator() noexcept                           = default;             //!< Defaulted.
    bcf_input_iterator(bcf_input_iterator const &) noexcept = default;             //!< Defaulted.
    bcf_input_iterator(bcf_input_iterator &&) noexcept      = default;             //!< Defaulted.
    bcf_input_iterator & operator=(bcf_input_iterator const &) noexcept = default; //!< Defaulted.
    bcf_input_iterator & operator=(bcf_input_iterator &&) noexcept = default;      //!< Defaulted.
    ~bcf_input_iterator() noexcept                                 = default;      //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit bcf_input_iterator(std::basic_streambuf<char> & ibuf, bool const init = true) :
      stream_buf{reinterpret_cast<bio::detail::stream_buffer_exposer<char> *>(&ibuf)}
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
        header.l_text = seqan3::detail::to_little_endian(*(reinterpret_cast<uint32_t *>(stream_buf->gptr() + 5)));
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

            l_shared = seqan3::detail::to_little_endian(*(reinterpret_cast<uint32_t *>(overflow_buffer.data())));
            l_indiv  = seqan3::detail::to_little_endian(*(reinterpret_cast<uint32_t *>(overflow_buffer.data() + 4)));
            overflow_buffer.clear();
        }
        else // read l_shared and l_indiv directly
        {
            l_shared = seqan3::detail::to_little_endian(*(reinterpret_cast<uint32_t *>(stream_buf->gptr())));
            l_indiv  = seqan3::detail::to_little_endian(*(reinterpret_cast<uint32_t *>(stream_buf->gptr() + 4)));
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

} // namespace bio::detail

namespace bio
{
/*!\brief Format input handler for the BCF format (bio::bcf).
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
 * |`print_warnings` |`bool`   | `false` | Whether to print non-critical warngings to seqan3::debug_stream  |
 *
 * ### Performance
 *
 * TODO after genotype redesign
 */
template <>
class format_input_handler<bcf> : public format_input_handler_base<format_input_handler<bcf>>
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
        std::string message = "[SeqAn3 BCF format error in record " + detail::to_string(record_no) + "] ";
        ((message += detail::to_string(messages)), ...);

        throw format_error{message};
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
    //!\brief Decode the type descripter byte into a detail::bcf_type_descriptor and a number.
    std::pair<detail::bcf_type_descriptor, uint8_t> decode_type_descriptor_byte(std::byte const b) const
    {
        uint8_t desc = static_cast<uint8_t>(b) & 0b00001111;
        uint8_t size = static_cast<uint8_t>(b) >> 4;

        if (desc == 4 || desc == 6 || desc >= 8)
            error("Cannot decode BCF type descriptor byte. Value: ", desc, " is unknown.");

        return {detail::bcf_type_descriptor{desc}, size};
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
            case detail::bcf_type_descriptor::int8:
                number = *reinterpret_cast<int8_t const *>(cache_ptr);
                cache_ptr += 1;
                break;
            case detail::bcf_type_descriptor::int16:
                number = *reinterpret_cast<int16_t const *>(cache_ptr);
                cache_ptr += 2;
                break;
            case detail::bcf_type_descriptor::int32:
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

        if (desc != detail::bcf_type_descriptor::char8)
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

        cache_ptr += real_size;

        std::span<std::byte const> ret{b, size_t(cache_ptr - b)};

        assert(cache_ptr <= &*file_it->first.end());

        return ret;
    }

    // IMPLEMENTATION NOTE: the following is the only decoder that is used in parsed record reading.
    //!\brief Decodes any range of integral and stores in any range of integral (no check of sizes).
    template <detail::back_insertable out_t>
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
            case detail::bcf_type_descriptor::missing:
                return;
            case detail::bcf_type_descriptor::int8:
                {
                    assert((cache_end_ptr - cache_ptr) == (ptrdiff_t)real_size);
                    std::span<int8_t const> data{reinterpret_cast<int8_t const *>(cache_ptr), real_size};
                    detail::sized_range_copy(data, out);
                    break;
                }
            case detail::bcf_type_descriptor::int16:
                {
                    assert((cache_end_ptr - cache_ptr) == (ptrdiff_t)real_size * 2);
                    std::span<int16_t const> data{reinterpret_cast<int16_t const *>(cache_ptr), real_size};
                    detail::sized_range_copy(data, out);
                    break;
                }
            case detail::bcf_type_descriptor::int32:
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
    using format_fields   = std::remove_cvref_t<decltype(var_io::default_field_ids)>;
    //!\brief Type of the raw record.
    using raw_record_type = record<format_fields, decltype(var_io::field_types_raw)>;

    //!\brief The raw record.
    raw_record_type raw_record;

    //!\brief This is cached during raw-record reading already.
    std::string_view                id_cache;
    //!\brief This is cached during raw-record reading already.
    std::string_view                ref_cache;
    //!\brief This is cached during raw-record reading already.
    std::vector<std::string_view>   alts_cache;
    //!\brief This is cached during raw-record reading already.
    detail::bcf_record_core const * record_core = nullptr;
    //!\brief Storage for GT values which are not encoded as strings (but pretend to be).
    std::vector<std::string>        gt_cache; // TODO concatenated_sequences
    //!\brief Temporary storage for string/number conversions.
    std::string                     number_cache;

    //!\brief This is only used as temporary storage when reading VCF style genotypes.
    std::vector<var_io::genotype_element_bcf<ownership::deep>>    deep_genotypes_cache;
    //!\brief This is only used as temporary storage when reading VCF style genotypes.
    std::vector<var_io::genotype_element_bcf<ownership::shallow>> shallow_genotypes_cache;

    //!\brief The header.
    var_io::header             header;
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

        record_core = reinterpret_cast<detail::bcf_record_core const *>(record_data.data());
        gt_cache.resize(record_core->n_sample);

        std::byte const * cache_ptr     = nullptr;
        std::byte const * cache_ptr_bak = nullptr;
        // offset, count
        get<field::chrom>(raw_record)   = record_data.subspan(0, 4);
        get<field::pos>(raw_record)     = record_data.subspan(4, 4);
        get<field::qual>(raw_record)    = record_data.subspan(12, 4);

        cache_ptr = record_data.data() + 24; // == sizeof(record_core)

        cache_ptr_bak              = cache_ptr;
        id_cache                   = decode_string(cache_ptr); // store string_view already
        get<field::id>(raw_record) = raw_f_t{cache_ptr_bak, size_t(cache_ptr - cache_ptr_bak)};

        cache_ptr_bak               = cache_ptr;
        ref_cache                   = decode_string(cache_ptr); // store string_view already
        get<field::ref>(raw_record) = raw_f_t{cache_ptr_bak, size_t(cache_ptr - cache_ptr_bak)};

        cache_ptr_bak = cache_ptr;
        for (size_t i = 1; i < record_core->n_allele; ++i)
            alts_cache.push_back(decode_string(cache_ptr));
        get<field::alt>(raw_record) = raw_f_t{cache_ptr_bak, size_t(cache_ptr - cache_ptr_bak)};

        get<field::filter>(raw_record) = decode_any_field(cache_ptr);

        size_t distance_to_genotypes = (record_data.data() + genotype_offset) - cache_ptr;
        get<field::info>(raw_record) = raw_f_t{cache_ptr, distance_to_genotypes};

        get<field::genotypes>(raw_record) = record_data.subspan(genotype_offset /*, end*/);
    }
    //!\}

    /*!\name Dynamic type initialisation and parsing
     * \{
     */
    //!\brief Set to single value.
    template <var_io::dynamic_type_id id_, typename elem_t>
    static inline void dynamic_type_init_single(std::byte const *& cache_ptr, auto & output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();
        output_                  = *reinterpret_cast<elem_t const *>(cache_ptr);
        cache_ptr += sizeof(elem_t);
    }

    //!\brief Set to string.
    template <var_io::dynamic_type_id id_>
    static inline void dynamic_type_init_string(size_t const size, std::byte const *& cache_ptr, auto & output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();

        std::string_view tmp{reinterpret_cast<char const *>(cache_ptr), size};
        detail::string_copy(tmp, output_);

        cache_ptr += size;
    }

    //!\brief Set to vector.
    template <var_io::dynamic_type_id id_, typename elem_t>
    static inline void dynamic_type_init_vector(size_t const size, std::byte const *& cache_ptr, auto & output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();

        std::span<elem_t const> tmp{reinterpret_cast<elem_t const *>(cache_ptr), size};
        detail::sized_range_copy(tmp, output_);

        cache_ptr += size * sizeof(elem_t);
    }

    //!\brief Set to vector-of-string.
    template <var_io::dynamic_type_id id_>
    static inline void dynamic_type_init_vector_of_string(size_t const       size,
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

    //!\brief Set to vector-of-vector.
    template <var_io::dynamic_type_id id_, typename elem_t>
    static inline void dynamic_type_init_vector_of_vector(size_t const       outer_size,
                                                          size_t const       inner_size,
                                                          std::byte const *& cache_ptr,
                                                          auto &             output)
    {
        constexpr size_t id      = static_cast<size_t>(id_);
        auto &           output_ = output.template emplace<id>();
        for (size_t i = 0; i < outer_size; ++i)
        {
            std::span<elem_t const> tmp{reinterpret_cast<elem_t const *>(cache_ptr) + i * inner_size, inner_size};

            // vectors can be smaller by being padded with end-of-vector values
            size_t s = tmp.size();
            while (s > 0 && tmp[s - 1] == detail::end_of_vector<elem_t>)
                --s;
            tmp = tmp.subspan(0, s);

            output_.emplace_back();
            detail::sized_range_copy(tmp, output_.back());
        }
        cache_ptr += outer_size * inner_size * sizeof(elem_t);
    }

    template <detail::is_dynamic_type dyn_t>
    void parse_dynamic_type(var_io::dynamic_type_id const     id_from_header,
                            detail::bcf_type_descriptor const desc,
                            size_t const                      size,
                            std::byte const *&                cache_ptr,
                            dyn_t &                           output); // implementation below class

    template <detail::is_dynamic_vector_type dyn_t>
    void parse_dynamic_type(var_io::dynamic_type_id const     id_from_header,
                            detail::bcf_type_descriptor const desc,
                            size_t const                      outer_size,
                            size_t const                      inner_size,
                            std::byte const *&                cache_ptr,
                            dyn_t &                           output); // implementation below class
    //!\}

    /*!\name Parsed record handling
     * \{
     */
    //!\brief Reading of CHROM field.
    void parse_field(vtag_t<field::chrom> const & /**/, std::integral auto & parsed_field)
    {
        parsed_field = record_core->chrom;
    }

    //!\overload
    void parse_field(vtag_t<field::chrom> const & /**/, auto & parsed_field)
    {
        parse_field_aux(header.contigs[header.idx_to_contig_pos().at(record_core->chrom)].id, parsed_field);
    }

    //!\brief Reading of POS field.
    void parse_field(vtag_t<field::pos> const & /**/, std::integral auto & parsed_field)
    {
        parsed_field = record_core->pos + 1; // one-based positions
    }

    //!\brief Reading of ID field.
    void parse_field(vtag_t<field::id> const & /**/, auto & parsed_field)
    {
        if (std::ranges::empty(id_cache))
            parse_field_aux(std::string_view{"."}, parsed_field);
        else
            parse_field_aux(id_cache, parsed_field);
    }

    //!\brief Reading of REF field.
    void parse_field(vtag_t<field::ref> const & /**/, auto & parsed_field) { parse_field_aux(ref_cache, parsed_field); }

    //!\brief Reading of ALT field.
    void parse_field(vtag_t<field::alt> const & /**/, std::vector<std::string_view> & parsed_field)
    {
        parsed_field = alts_cache;
    }

    //!\overload
    template <detail::back_insertable parsed_field_t>
        requires(std::ranges::range<std::ranges::range_reference_t<parsed_field_t>>)
    void parse_field(vtag_t<field::alt> const & /**/, parsed_field_t & parsed_field)
    {
        for (std::string_view const alt : alts_cache)
        {
            std::ranges::range_value_t<parsed_field_t> out;

            parse_field_aux(alt, out);

            parsed_field.push_back(std::move(out));
        }
    }

    //!\brief Reading of CHROM field.
    void parse_field(vtag_t<field::qual> const & /**/, seqan3::arithmetic auto & parsed_field)
    {
        parsed_field = record_core->qual;
    }

    //!\brief Reading of FILTER field.
    template <detail::back_insertable parsed_field_t>
        requires detail::int_range<parsed_field_t>
    void parse_field(vtag_t<field::filter> const & /**/, parsed_field_t & parsed_field)
    {
        decode_numbers_into(get<field::filter>(raw_record), parsed_field);
    }

    //!\overload
    template <detail::back_insertable parsed_field_t>
        requires detail::char_range<std::ranges::range_reference_t<parsed_field_t>>
    void parse_field(vtag_t<field::filter> const & /**/, parsed_field_t & parsed_field)
    {
        std::vector<int32_t> tmp;

        decode_numbers_into(get<field::filter>(raw_record), tmp);

        for (int32_t const idx : tmp)
        {
            std::ranges::range_value_t<parsed_field_t> out;
            parse_field_aux(header.filters[header.idx_to_filter_pos().at(idx)].id, out);
            parsed_field.push_back(std::move(out));
        }
    }

    //!\brief Reading of the INFO field.
    template <detail::back_insertable parsed_field_t>
        requires detail::info_element_concept<std::ranges::range_reference_t<parsed_field_t>>
    void parse_field(vtag_t<field::info> const & /**/, parsed_field_t & parsed_field)
    {
        std::span<std::byte const> raw_field = get<field::info>(raw_record);
        std::byte const *          cache_ptr = raw_field.data();

        for (size_t i = 0; i < record_core->n_info; ++i)
        {
            parsed_field.push_back({});

            auto & [id, variant] = parsed_field.back();
            ;
            int32_t idx = decode_integral(cache_ptr);

            assert(cache_ptr < raw_field.data() + raw_field.size());

            auto [desc, size] = decode_type_descriptor_byte(*cache_ptr);
            ++cache_ptr;

            size_t real_size = size < 15 ? size : decode_integral(cache_ptr);

            parse_dynamic_type(header.infos[header.idx_to_info_pos().at(idx)].type,
                               desc,
                               real_size,
                               cache_ptr,
                               variant);

            if constexpr (detail::char_range<decltype(id)>)
                detail::string_copy(header.infos[header.idx_to_info_pos().at(idx)].id, id);
            else
                id = idx;
        }

        assert(cache_ptr == raw_field.data() + raw_field.size());
    }

    //!\brief Auxilliary function for reading GT which isn't actually encoded as a string.
    template <typename int_t>
    void parse_gt_field(std::span<int_t const> const in, std::string & number_cache, std::string & output) const
    {
        bool first = true;

        for (int_t i : in)
        {
            if (i == detail::end_of_vector<int_t>)
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

    //!\brief Implementation for parsing into bcf-style genotypes.
    void parse_genotypes_impl(auto & parsed_field)
    {
        std::span<std::byte const> raw_field = get<field::genotypes>(raw_record);

        std::byte const * cache_ptr = raw_field.data();

        std::string number_cache;

        for (size_t i = 0; i < static_cast<size_t>(record_core->n_fmt); ++i)
        {
            parsed_field.emplace_back();

            auto & [parsed_idx, parsed_variant] = parsed_field.back();

            int32_t fmt_key                 = decode_integral(cache_ptr);
            parsed_idx                      = fmt_key;
            var_io::header::format_t format = header.formats[header.idx_to_format_pos().at(fmt_key)];

            auto [fmt_type, fmt_size] = decode_type_descriptor_byte(*cache_ptr);
            ++cache_ptr;

            fmt_size = fmt_size < 15 ? fmt_size : decode_integral(cache_ptr);

            if (format.id == "GT") // this needs custom decoding, it is not a string
            {
                for (size_t sample = 0; sample < record_core->n_sample; ++sample)
                {
                    gt_cache[sample].clear();
                    number_cache.clear();

                    switch (fmt_type)
                    {
                        case detail::bcf_type_descriptor::int8:
                            {
                                std::span<int8_t const> data{reinterpret_cast<int8_t const *>(cache_ptr), fmt_size};
                                parse_gt_field(data, number_cache, gt_cache[sample]);
                                cache_ptr += fmt_size;
                                break;
                            }
                        case detail::bcf_type_descriptor::int16:
                            {
                                std::span<int16_t const> data{reinterpret_cast<int16_t const *>(cache_ptr), fmt_size};
                                parse_gt_field(data, number_cache, gt_cache[sample]);
                                fmt_size += fmt_size * 2;
                                break;
                            }
                        case detail::bcf_type_descriptor::int32:
                            {
                                std::span<int32_t const> data{reinterpret_cast<int32_t const *>(cache_ptr), fmt_size};
                                parse_gt_field(data, number_cache, gt_cache[sample]);
                                fmt_size += fmt_size * 4;
                                break;
                            }
                        default:
                            error("GT field must always be encoded as number but was something else.");
                            break;
                    }
                }

                constexpr size_t string_id = static_cast<size_t>(var_io::dynamic_type_id::string);
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
                parse_dynamic_type(format.type, fmt_type, record_core->n_sample, fmt_size, cache_ptr, parsed_variant);
            }
        }

        assert(cache_ptr == raw_field.data() + raw_field.size());
    }

    //!\brief Reading of the GENOTYPES field (BCF-style).
    template <detail::back_insertable field_t>
        requires detail::genotype_bcf_style_concept<std::ranges::range_reference_t<field_t>>
    void parse_field(vtag_t<field::genotypes> const & /**/, field_t & parsed_field)
    {
        parse_genotypes_impl(parsed_field);
    }

    //!\brief Implementation for parsing into vcf-style genotypes.
    void parse_genotypes_impl_vcf_style(auto & genotypes_cache, auto & parsed_field)
    {
        genotypes_cache.clear();
        parse_genotypes_impl(genotypes_cache); // we parse into BCF-style caches and convert

        auto & [parsed_formats, parsed_samples] = parsed_field;

        /* formats */
        for (auto & [idx, dyn_type] : genotypes_cache)
        {
            var_io::header::format_t & format = header.formats[header.idx_to_format_pos().at(idx)];

            parsed_formats.push_back({});
            detail::string_copy(format.id, parsed_formats.back());
        }

        /* samples */
        parsed_samples.resize(record_core->n_sample);
        for (size_t s = 0; s < record_core->n_sample; ++s)
        {
            parsed_samples[s].resize(parsed_formats.size());
            for (size_t f = 0; f < parsed_formats.size(); ++f)
            {
                auto & output  = parsed_samples[s][f];
                auto   visitor = [s, &output](auto & in) { output = std::move(in[s]); };

                std::visit(visitor, detail::get_second(genotypes_cache[f]));
            }
        }
    }

    //!\brief Reading of the GENOTYPES field (VCF-style).
    void parse_field(vtag_t<field::genotypes> const & /**/, detail::genotypes_vcf_style_concept auto & parsed_field)
    {
        using parsed_field_t = decltype(parsed_field);
        using dyn_t = std::ranges::range_value_t<std::ranges::range_reference_t<detail::second_elem_t<parsed_field_t>>>;

        if constexpr (std::same_as<var_io::dynamic_type<ownership::shallow>, dyn_t>)
            parse_genotypes_impl_vcf_style(shallow_genotypes_cache, parsed_field);
        else
            parse_genotypes_impl_vcf_style(deep_genotypes_cache, parsed_field);
    }

    //!\brief Overload for parsing the private data.
    void parse_field(vtag_t<field::_private> const & /**/, var_io::record_private_data & parsed_field)
    {
        parsed_field = var_io::record_private_data{.header_ptr = &header};
    }
    //!\}

public:
    /*!\name Constructors, destructor and assignment.
     * \{
     */
    format_input_handler()                             = default;            //!< Defaulted.
    format_input_handler(format_input_handler const &) = delete;             //!< Deleted.
    format_input_handler(format_input_handler &&)      = default;            //!< Defaulted.
    ~format_input_handler()                            = default;            //!< Defaulted.
    format_input_handler & operator=(format_input_handler const &) = delete; //!< Deleted.
    format_input_handler & operator=(format_input_handler &&) = default;     //!< Defaulted.

    /*!\brief Construct with an options object.
     * \param[in,out] str The input stream.
     * \param[in] options An object with options for the input handler.
     * \details
     *
     * The options argument is typically bio::var_io::reader_options, but any object with a subset of similarly named
     * members is also accepted. See bio::format_input_handler<bio::bcf> for the supported options and defaults.
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
        header = var_io::header{std::move(text_tmp)};

        /* checks on header */
        if (header.string_to_format_pos().contains("GT") &&
            header.formats[header.string_to_format_pos().at("GT")].type != var_io::dynamic_type_id::string)
        {
            error("The \"GT\" field must always be encoded as a string.");
        }
        // TODO more checks
    }

    //!\brief Construct with only an input stream.
    format_input_handler(std::istream & str) : format_input_handler
    {
        str, int {}
    }
    {}
    //!\}

    //!\brief Return a reference to the header contained in the input handler.
    var_io::header const & get_header() const { return header; }
};

/*!\brief Parse a "dynamically typed" field out of a BCF stream and store the content in a variant.
 * \tparam dyn_t             Type of the variant; specialisation of seqan3::var_io::dynamic_type.
 * \param[in] id_from_header A value of seqan3::var_io::dynamic_type_id that denotes the expected type.
 * \param[in] desc           A value of seqan3::detail::bcf_type_descriptor that notes the detected type.
 * \param[in] size           The number of values belonging to this field.
 * \param[in,out] cache_ptr  Pointer into the BCF stream; will be updated to point past the end of read data.
 * \param[out] output        The variant to hold the parsed value.
 */
template <detail::is_dynamic_type dyn_t>
inline void format_input_handler<bcf>::parse_dynamic_type(var_io::dynamic_type_id const     id_from_header,
                                                          detail::bcf_type_descriptor const desc,
                                                          size_t const                      size,
                                                          std::byte const *&                cache_ptr,
                                                          dyn_t &                           output)
{
    // TODO DRY out the boilerplate error messages
    if (static_cast<size_t>(id_from_header) < 3 /* char8, int32, float32 */ && size != 1)
        error("BCF data field expected exactly one element, got:", size);

    switch (id_from_header)
    {
        case var_io::dynamic_type_id::char8:
            {
                if (desc != detail::bcf_type_descriptor::char8)
                    error("Attempting to create char but the byte descriptor does not indicate char type.");

                dynamic_type_init_single<var_io::dynamic_type_id::char8, char>(cache_ptr, output);
                return;
            }
        case var_io::dynamic_type_id::int32:
            {
                switch (desc)
                {
                    case detail::bcf_type_descriptor::int8:
                        dynamic_type_init_single<var_io::dynamic_type_id::int32, int8_t>(cache_ptr, output);
                        break;
                    case detail::bcf_type_descriptor::int16:
                        dynamic_type_init_single<var_io::dynamic_type_id::int32, int16_t>(cache_ptr, output);
                        break;
                    case detail::bcf_type_descriptor::int32:
                        dynamic_type_init_single<var_io::dynamic_type_id::int32, int32_t>(cache_ptr, output);
                        break;
                    default:
                        error("Attempting to create int but the byte descriptor does not indicate int type.");
                }
                return;
            }
        case var_io::dynamic_type_id::float32:
            {
                if (desc != detail::bcf_type_descriptor::float32)
                    error("Attempting to create float but the byte descriptor does not indicate float type.");

                dynamic_type_init_single<var_io::dynamic_type_id::float32, float>(cache_ptr, output);
                return;
            }
        case var_io::dynamic_type_id::string:
            {
                if (desc != detail::bcf_type_descriptor::char8)
                    error("Attempting to creates string but the byte descriptor does not indicate string type.");

                dynamic_type_init_string<var_io::dynamic_type_id::string>(size, cache_ptr, output);
                return;
            }
        case var_io::dynamic_type_id::vector_of_char8:
            {
                if (desc != detail::bcf_type_descriptor::char8)
                    error("Attempting to create vector of char but the byte descriptor does not indicate char type.");

                dynamic_type_init_vector<var_io::dynamic_type_id::vector_of_char8, char>(size, cache_ptr, output);
                return;
            }
        case var_io::dynamic_type_id::vector_of_int8:
        case var_io::dynamic_type_id::vector_of_int16:
        case var_io::dynamic_type_id::vector_of_int32:
            {
                switch (desc)
                {
                    case detail::bcf_type_descriptor::int8:
                        {
                            dynamic_type_init_vector<var_io::dynamic_type_id::vector_of_int8, int8_t>(size,
                                                                                                      cache_ptr,
                                                                                                      output);
                            break;
                        }
                    case detail::bcf_type_descriptor::int16:
                        {
                            dynamic_type_init_vector<var_io::dynamic_type_id::vector_of_int16, int16_t>(size,
                                                                                                        cache_ptr,
                                                                                                        output);
                            break;
                        }
                    case detail::bcf_type_descriptor::int32:
                        {
                            dynamic_type_init_vector<var_io::dynamic_type_id::vector_of_int32, int32_t>(size,
                                                                                                        cache_ptr,
                                                                                                        output);
                            break;
                        }
                    default:
                        error("Attempting to create vector of int but the byte descriptor does not indicate int type.");
                }
                return;
            }
        case var_io::dynamic_type_id::vector_of_float32:
            {
                if (desc != detail::bcf_type_descriptor::float32)
                    error("Attempting to create vector of float but the byte descriptor does not indicate float type.");

                dynamic_type_init_vector<var_io::dynamic_type_id::vector_of_float32, float>(size, cache_ptr, output);
                return;
            }
        case var_io::dynamic_type_id::vector_of_string:
            {
                if (desc != detail::bcf_type_descriptor::char8)
                    error(
                      "Attempting to create vector of string but the byte descriptor does not indicate char alphabet");

                dynamic_type_init_vector_of_string<var_io::dynamic_type_id::vector_of_string>(size, cache_ptr, output);
                return;
            }
        case var_io::dynamic_type_id::flag:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::flag);
                output.template emplace<id>(true);

                cache_ptr += size; // This should be 0, but is allowed to be something else
                return;
            }
    }
}

/*!\brief Parse a "dynamically typed" field out of a BCF stream and store the content in a vector-variant.
 * \tparam dyn_t             Type of the variant; specialisation of seqan3::var_io::dynamic_type.
 * \param[in] id_from_header A value of seqan3::var_io::dynamic_type_id that denotes the expected type.
 * \param[in] desc           A value of seqan3::detail::bcf_type_descriptor that notes the detected type.
 * \param[in] outer_size     The number of values belonging to this field.
 * \param[in] inner_size     The number of values per inner vector in case of vector-of-vector.
 * \param[in,out] cache_ptr  Pointer into the BCF stream; will be updated to point past the end of read data.
 * \param[out] output        The variant to hold the parsed value.
 */
template <detail::is_dynamic_vector_type dyn_t>
inline void format_input_handler<bcf>::parse_dynamic_type(var_io::dynamic_type_id const     id_from_header,
                                                          detail::bcf_type_descriptor const desc,
                                                          size_t const                      outer_size,
                                                          size_t const                      inner_size,
                                                          std::byte const *&                cache_ptr,
                                                          dyn_t &                           output)
{
    if (static_cast<size_t>(id_from_header) < 3 /* char8, int32, float32 */ && inner_size != 1)
        error("BCF data field expected exactly one element, got:", inner_size);

    switch (id_from_header)
    {
        case var_io::dynamic_type_id::char8:
            {
                if (desc != detail::bcf_type_descriptor::char8)
                    error("Attempting to create char but the byte descriptor does not indicate char type.");

                dynamic_type_init_vector<var_io::dynamic_type_id::char8, char>(outer_size, cache_ptr, output);
                return;
            }
        case var_io::dynamic_type_id::int32:
            {
                switch (desc)
                {
                    case detail::bcf_type_descriptor::int8:
                        {
                            dynamic_type_init_vector<var_io::dynamic_type_id::int32, int8_t>(outer_size,
                                                                                             cache_ptr,
                                                                                             output);
                            break;
                        }
                    case detail::bcf_type_descriptor::int16:
                        {
                            dynamic_type_init_vector<var_io::dynamic_type_id::int32, int16_t>(outer_size,
                                                                                              cache_ptr,
                                                                                              output);
                            break;
                        }
                    case detail::bcf_type_descriptor::int32:
                        {
                            dynamic_type_init_vector<var_io::dynamic_type_id::int32, int32_t>(outer_size,
                                                                                              cache_ptr,
                                                                                              output);
                            break;
                        }
                    default:
                        error("Attempting to create int but the byte descriptor does not indicate int type.");
                }
                return;
            }
        case var_io::dynamic_type_id::float32:
            {
                if (desc != detail::bcf_type_descriptor::float32)
                    error("Attempting to create float but the byte descriptor does not indicate float type.");

                dynamic_type_init_vector<var_io::dynamic_type_id::float32, float>(outer_size, cache_ptr, output);
                return;
            }
        case var_io::dynamic_type_id::string:
            {
                if (desc != detail::bcf_type_descriptor::char8)
                    error("Attempting to creates string but the byte descriptor does not indicate string type.");

                dynamic_type_init_vector_of_string<var_io::dynamic_type_id::string>(outer_size, cache_ptr, output);
                return;
            }
        case var_io::dynamic_type_id::vector_of_char8:
            {
                if (desc != detail::bcf_type_descriptor::char8)
                    error("Attempting to create vector of char but the byte descriptor does not indicate char type.");

                dynamic_type_init_vector_of_vector<var_io::dynamic_type_id::vector_of_char8, char>(outer_size,
                                                                                                   inner_size,
                                                                                                   cache_ptr,
                                                                                                   output);
                return;
            }
        case var_io::dynamic_type_id::vector_of_int8:
        case var_io::dynamic_type_id::vector_of_int16:
        case var_io::dynamic_type_id::vector_of_int32:
            {
                switch (desc)
                {
                    case detail::bcf_type_descriptor::int8:
                        {
                            dynamic_type_init_vector_of_vector<var_io::dynamic_type_id::vector_of_int8, int8_t>(
                              outer_size,
                              inner_size,
                              cache_ptr,
                              output);
                            break;
                        }
                    case detail::bcf_type_descriptor::int16:
                        {
                            dynamic_type_init_vector_of_vector<var_io::dynamic_type_id::vector_of_int16, int16_t>(
                              outer_size,
                              inner_size,
                              cache_ptr,
                              output);
                            break;
                        }
                    case detail::bcf_type_descriptor::int32:
                        {
                            dynamic_type_init_vector_of_vector<var_io::dynamic_type_id::vector_of_int32, int32_t>(
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
        case var_io::dynamic_type_id::vector_of_float32:
            {
                if (desc != detail::bcf_type_descriptor::float32)
                    error("Attempting to create vector of float but the byte descriptor does not indicate float type.");

                dynamic_type_init_vector_of_vector<var_io::dynamic_type_id::vector_of_float32, float>(outer_size,
                                                                                                      inner_size,
                                                                                                      cache_ptr,
                                                                                                      output);
                return;
            }
        case var_io::dynamic_type_id::vector_of_string:
            {
                if (desc != detail::bcf_type_descriptor::char8)
                    error(
                      "Attempting to create vector of string but the byte descriptor does not indicate char alphabet");

                // TODO this definitely needs a test
                constexpr size_t id      = static_cast<size_t>(var_io::dynamic_type_id::vector_of_string);
                auto &           output_ = output.template emplace<id>();
                std::string_view tmp{reinterpret_cast<char const *>(cache_ptr), outer_size * inner_size};
                for (size_t i = 0; i < outer_size; ++i)
                {
                    output_.emplace_back();

                    std::string_view tmp_inner = tmp.substr(i * inner_size, (i + 1) * inner_size);
                    for (std::string_view const s : tmp_inner | detail::eager_split(','))
                    {
                        output_.back().push_back({});
                        detail::string_copy(s, output_.back().back());
                    }
                }
                cache_ptr += outer_size * inner_size;
                return;
            }
        case var_io::dynamic_type_id::flag:
            {
                error("seqan3::var_io::dynamic_vector_type cannot be initialised to flag state.");
                return;
            }
    }
}

} // namespace bio
