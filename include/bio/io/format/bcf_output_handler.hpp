// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::format_output_handler<bcf>.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bit>
#include <numeric>

#include <bio/meta/tag/vtag.hpp>

#include <bio/io/detail/magic_get.hpp>
#include <bio/io/detail/misc.hpp>
#include <bio/io/format/bcf.hpp>
#include <bio/io/format/format_output_handler.hpp>
#include <bio/io/misc/char_predicate.hpp>
#include <bio/io/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/io/var_io/header.hpp>
#include <bio/io/var_io/misc.hpp>
#include <bio/io/var_io/writer_options.hpp>

namespace bio::io
{

/*!\brief Format output handler for the VCF format (bio::io::bcf).
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
 * ### Performance
 *
 * TODO after genotype redesign
 */
template <>
class format_output_handler<bcf>
{
private:
    //!\brief Print an error message with current record number in diagnostic.
    [[noreturn]] void error(auto const &... messages) const
    {
        std::string message = "[SeqAn3 BCF format error writing record " + detail::to_string(record_no) + "] ";
        ((message += detail::to_string(messages)), ...);

        throw format_error{message};
    }

    /*!\name State
     * \{
     */
    //!\brief Whether the header has been written or not.
    bool                                                                    header_has_been_written = false;
    //!\brief Pointer to header that can be owning or non-owning.
    std::unique_ptr<var_io::header const, void (*)(var_io::header const *)> header                  = {nullptr,
                                                                                                       [](var_io::header const *) {}};
    //!\brief This is the descriptor used for IDX values [stored here, so it doesn't need to be recomputed].
    detail::bcf_type_descriptor                                             idx_desc;

    //!\brief Index of the current record.
    size_t record_no = 0;

    //!\brief Part of the record that can be written en-bloc.
    detail::bcf_record_core record_core;

    //!\brief A pointer to the output stream.
    std::ostream *                        stream = nullptr;
    //!\brief An intermediate stream that is guaranteed to always hold the complete current record in memory.
    std::ostringstream                    buffer_stream;
    //!\brief Accessor for the the buffer_stream.
    detail::stream_buffer_exposer<char> * streambuf_exposer =
      reinterpret_cast<detail::stream_buffer_exposer<char> *>(buffer_stream.rdbuf());
    //!\brief The iterator on that stream.
    detail::fast_ostreambuf_iterator<char> it{buffer_stream};

    //!\brief Distance from pbase()
    size_t this_record_offset = 0;

    //!\brief Used when writing vector-of-string.
    std::vector<size_t> string_size_buffer;

    //!\brief Set if this object was moved from.
    bool moved_from = false;
    //!\}

    /*!\name Options
     * \{
     */
    //!\brief Try to use the smallest possible integer type (creates smaller files but is potentially slower).
    bool compress_integers   = true;
    //!\brief Throw exceptions if the header type does not match the descriptor from the header.
    bool verify_header_types = false;
    //!\}

    /*!\name Arbitrary helpers
     * \{
     */
    //!\brief Write the type-descriptor byte(s).
    void write_type_descriptor(detail::bcf_type_descriptor const desc, size_t const size)
    {
        if (size < 15)
        {
            uint8_t byte = static_cast<uint8_t>(size) << 4;
            byte |= static_cast<uint8_t>(desc);
            it->write_as_binary(byte);
        }
        else
        {
            uint8_t byte = static_cast<uint8_t>(15) << 4;
            byte |= static_cast<uint8_t>(desc);
            it->write_as_binary(byte);
            detail::bcf_type_descriptor int_desc = detail::smallest_int_desc(size);
            write_type_descriptor1(int_desc);
            write_single_impl(size, int_desc);
        }
    }

    //!\brief Overload that is hardcoded for size==1
    void write_type_descriptor1(detail::bcf_type_descriptor const desc)
    {
        uint8_t byte = static_cast<uint8_t>(1) << 4;
        byte |= static_cast<uint8_t>(desc);
        it->write_as_binary(byte);
    }

    //!\brief Write a single value.
    void write_single_impl(meta::arithmetic auto num, detail::bcf_type_descriptor const desc)
    {
        switch (desc)
        {
            case detail::bcf_type_descriptor::char8:
                it->write_as_binary(static_cast<char>(num));
                break;
            case detail::bcf_type_descriptor::int8:
                it->write_as_binary(static_cast<int8_t>(num));
                break;
            case detail::bcf_type_descriptor::int16:
                it->write_as_binary(static_cast<int16_t>(num));
                break;
            case detail::bcf_type_descriptor::int32:
                it->write_as_binary(static_cast<int32_t>(num));
                break;
            case detail::bcf_type_descriptor::float32:
                it->write_as_binary(static_cast<float>(num));
                break;
            default:
                error("Trying to write an unknown type");
        }
    }

    /*!\brief Transform the given range into a run-length encoded one and write.
     * \param[in] range The data.
     * \param[in] desc The type descriptor.
     * \param[in,out] fun A functor for writing the value part of the (run_length, value)-pair.
     *
     * The functor takes the value and the descriptor as arguments.
     */
    void run_length_encode(auto & range, detail::bcf_type_descriptor const desc, auto && fun)
    {
        size_t run_length = 1;

        auto eq = meta::overloaded{[](std::ranges::range auto && lhs, std::ranges::range auto && rhs)
                                   { return std::ranges::equal(lhs, rhs); },
                                   [](auto && lhs, auto && rhs) { return lhs == rhs; }};

        for (auto rit = std::ranges::begin(range), next = std::ranges::next(rit); rit != std::ranges::end(range);
             ++rit, ++next)
        {
            if (next != std::ranges::end(range) && eq(*rit, *next))
            {
                ++run_length;
            }
            else
            {
                // write size
                detail::bcf_type_descriptor int_desc = detail::smallest_int_desc(run_length);
                write_type_descriptor1(int_desc);
                write_single_impl(run_length, int_desc);

                // write value
                fun(*rit, desc);

                // finish
                run_length = 1;
            }
        }
    }

    //!\brief Write a range of values.
    void write_range_impl(std::ranges::forward_range auto && range, detail::bcf_type_descriptor const desc)
    {
        using data_t = decltype(range);
        using elem_t = std::ranges::range_value_t<data_t>;

        auto dump_or_convert = [&]<typename target_t>(target_t)
        {
            if constexpr (std::same_as<target_t, elem_t> && std::ranges::contiguous_range<data_t> &&
                          std::ranges::sized_range<data_t>)
            {
                it->write_as_binary(range);
            }
            else if constexpr (std::same_as<target_t, elem_t> && std::same_as<target_t, char>)
            {
                it->write_range(range);
            }
            else if constexpr (detail::deliberate_alphabet<elem_t>)
            {
                it->write_range(range | bio::views::to_char);
            }
            else
            {
                // TODO this is probably not the most efficient implementation
                for (elem_t const elem : range)
                {
                    target_t buf = (elem == var_io::missing_value<elem_t>) ? var_io::missing_value<target_t>
                                                                           : static_cast<target_t>(elem);
                    it->write_as_binary(buf);
                }
            }
        };

        switch (desc)
        {
            case detail::bcf_type_descriptor::char8:
                dump_or_convert(char{});
                break;
            case detail::bcf_type_descriptor::int8:
                dump_or_convert(int8_t{});
                break;
            case detail::bcf_type_descriptor::int16:
                dump_or_convert(int16_t{});
                break;
            case detail::bcf_type_descriptor::int32:
                dump_or_convert(int32_t{});
                break;
            case detail::bcf_type_descriptor::float32:
                dump_or_convert(float{});
                break;
            default:
                error("Trying to write an unknown type");
        }
    }

    //!\brief Write a range of ranges.
    template <std::ranges::forward_range rng_t>
        requires std::ranges::forward_range<std::ranges::range_reference_t<rng_t>>
    void write_range_impl(rng_t & range, [[maybe_unused]] detail::bcf_type_descriptor const desc)
    {
        assert(desc == detail::bcf_type_descriptor::char8);
        using alph_t = std::ranges::range_value_t<std::ranges::range_reference_t<rng_t>>;

        bool first = true;

        // TODO it would be good to have a views::join_with and even better one with size
        for (auto && r : range)
        {
            if (first)
                first = false;
            else
                it = ',';

            if constexpr (std::same_as<alph_t, char>)
                it->write_range(r);
            else if constexpr (detail::deliberate_alphabet<alph_t>)
                it->write_range(r | bio::views::to_char);
            else
                static_assert(std::same_as<alph_t, char>, "Can't handle this alphabet type here.");
        }
    }

    //!\brief Write the GT field (has custom writer, because it pretends to be a string but is encoded differently).
    size_t write_GT_impl(std::string_view const gt_string, detail::bcf_type_descriptor const desc)
    {
        bool phased = false;

        size_t n_alleles = 0;

        // TODO use views::eager_split once that can handle predicates
        for (size_t i = 0, j = 0; i <= gt_string.size(); ++i)
        {
            if (i == gt_string.size() || (is_char<'/'> || is_char<'|'>)(gt_string[i]))
            {
                std::string_view substr = gt_string.substr(j, i - j);

                size_t encoded_buf = 0;
                if (substr != ".")
                {
                    size_t buf = 0;
                    detail::string_to_number(substr, buf);
                    // encode according to BCF spec
                    encoded_buf = (buf + 1) << 1 | (size_t)phased;
                }
                // else it remains 0

                switch (desc)
                {
                    case detail::bcf_type_descriptor::int8:
                        assert(encoded_buf < std::numeric_limits<int8_t>::max());
                        it->write_as_binary(static_cast<int8_t>(encoded_buf));
                        break;
                    case detail::bcf_type_descriptor::int16:
                        assert(encoded_buf < std::numeric_limits<int16_t>::max());
                        it->write_as_binary(static_cast<int16_t>(encoded_buf));
                        break;
                    case detail::bcf_type_descriptor::int32:
                        assert(encoded_buf < std::numeric_limits<int32_t>::max());
                        it->write_as_binary(static_cast<int32_t>(encoded_buf));
                        break;
                    default:
                        error("Trying to write an unknown type.");
                }

                ++n_alleles;
                j      = i + 1;
                phased = is_char<'|'>(gt_string[i]);
            }
        }

        return n_alleles;
    }

    /*!\brief Write padding values to make vector have fixed size.
     * \param[in] num  Number of padding values to write.
     * \param[in] desc The type of values to write.
     * \param[in] front_missing Whether to add an initial missing value [see below].
     * \details
     *
     * This function adds padding values. For an existing vector, num times the detail::end_of_vector
     * value ("EOV") is written.
     *
     * In the case where a vector in a vector-of-vectors is completely absent, you need to set
     * front_missing. In that case the first value written will be bio::io::var_io::missing_value instead of EOV.
     */
    void write_range_padding(size_t num, detail::bcf_type_descriptor const desc, bool front_missing)
    {
        if (num == 0)
            return;

        auto func = [&]<typename target_t>(target_t)
        {
            if (front_missing)
            {
                it->write_as_binary(var_io::missing_value<target_t>);
                --num;
            }

            // TODO this could be optimised with some manner of views::repeat_n | views::convert | views::join
            for (size_t i = 0; i < num; ++i)
                it->write_as_binary(detail::end_of_vector<target_t>);
        };

        switch (desc)
        {
            case detail::bcf_type_descriptor::char8:
                func(char{});
                break;
            case detail::bcf_type_descriptor::int8:
                func(int8_t{});
                break;
            case detail::bcf_type_descriptor::int16:
                func(int16_t{});
                break;
            case detail::bcf_type_descriptor::int32:
                func(int32_t{});
                break;
            case detail::bcf_type_descriptor::float32:
                func(float{});
                break;
            default:
                error("Trying to write an unknown type");
        }
    }

    //!\brief Generic writing function that dispatches depending on the type of the argument.
    template <typename elem_t>
        //!\cond REQ
        requires(alphabet::alphabet<elem_t> || meta::arithmetic<elem_t>)
    //!\endcond
    void write_typed_data(elem_t const num, detail::bcf_type_descriptor const desc)
    {
        write_type_descriptor1(desc);
        write_single_impl(num, desc);
    }

    //!\overload
    void write_typed_data(detail::char_range auto && data, [[maybe_unused]] detail::bcf_type_descriptor const desc)
    {
        assert(desc == detail::bcf_type_descriptor::char8);

        if (std::ranges::empty(data) || std::ranges::equal(data, std::string_view{"."}))
        {
            it->write_as_binary(uint8_t{0x07});
        }
        else
        {
            size_t size = std::ranges::distance(data);
            write_type_descriptor(detail::bcf_type_descriptor::char8, size);
            it->write_range(data);
        }
    }

    //!\overload
    void write_typed_data(char const * const cstring, detail::bcf_type_descriptor const desc)
    {
        assert(desc == detail::bcf_type_descriptor::char8);
        return write_typed_data(std::string_view{cstring}, desc);
    }

    //!\overload
    void write_typed_data(std::ranges::forward_range auto && data, detail::bcf_type_descriptor const desc)
    {
        static_assert(!std::ranges::range<std::ranges::range_reference_t<decltype(data)>>);
        write_type_descriptor(desc, std::ranges::distance(data));
        write_range_impl(data, desc);
    }

    //!\overload
    template <std::ranges::forward_range data_t>
        //!\cond REQ
        requires detail::char_range_or_cstring<std::ranges::range_value_t<data_t>>
    //!\endcond
    void write_typed_data(data_t && vector_of_string, [[maybe_unused]] detail::bcf_type_descriptor const desc)
    {
        assert(desc == detail::bcf_type_descriptor::char8);

        size_t num_strings = std::ranges::distance(vector_of_string);
        if (num_strings == 0)
        {
            it->write_as_binary(uint8_t{0x07});
            return;
        }

        auto to_size = meta::overloaded{[](char const * const cstr) { return std::string_view{cstr}.size(); },
                                        [](std::ranges::range auto & rng) { return std::ranges::distance(rng); }};

        size_t size_sum = std::transform_reduce(std::ranges::begin(vector_of_string),
                                                std::ranges::end(vector_of_string),
                                                0ull,
                                                std::plus<>{},
                                                to_size);
        // add additional size for the "," that will be interspersed
        size_sum += num_strings - 1;

        write_type_descriptor(desc, size_sum);
        write_range_impl(vector_of_string, desc);
    }

    //!\overload
    void write_typed_data(bool, [[maybe_unused]] detail::bcf_type_descriptor const desc)
    {
        assert(desc == detail::bcf_type_descriptor::int8);
        // This is the behaviour according to spec:
        // write_typed_data(int8_t{1}, false);
        // but bcftools and htslib expect this:
        it = '\0';
    }

    //!\brief This overload adds an automatically deduced descriptor.
    void write_typed_data(auto && num)
    {
        return write_typed_data(num, detail::type_2_bcf_type_descriptor<std::remove_cvref_t<decltype(num)>>);
    }
    //!\}

    /*!\name Core record setters
     * \brief These set up the core record.
     * \{
     */
    //!\brief Overload for CHROM and numeric IDs.
    void set_core_chrom(std::integral auto & field) { record_core.chrom = static_cast<int32_t>(field); }

    //!\brief Overload for CHROM and text IDs.
    void set_core_chrom(detail::char_range_or_cstring auto && field)
    {
        std::string_view const buf = detail::to_string_view(field);
        if (auto it = header->string_to_contig_pos().find(buf); it == header->string_to_contig_pos().end())
            error("The contig '", field, "' is not present in the header.");
        else
            record_core.chrom = header->contigs[it->second].idx;
    }

    //!\overload
    void set_core_chrom(char const * const field) { set_core_chrom(std::string_view{field}); }

    //!\brief Overload for POS.
    void set_core_pos(std::integral auto & field)
    {
        assert(field >= 1);
        record_core.pos = static_cast<int32_t>(field - 1); // BCF is 0-based, we are 1-based
    }

    //!\brief Overload for rlen.
    void set_core_rlen(auto && field) { record_core.rlen = static_cast<int32_t>(std::ranges::distance(field)); }

    //!\overload
    void set_core_rlen(char const * const field) { set_core_rlen(std::string_view{field}); }

    //!\brief Overload for QUAL.
    void set_core_qual(meta::arithmetic auto & field) { record_core.qual = static_cast<float>(field); }

    //!\brief Overload for n_info.
    void set_core_n_info(auto & field)
    {
        using field_t = decltype(field);

        if constexpr (detail::info_element_writer_concept<field_t>)
            record_core.n_info = 1;
        else
            record_core.n_info = detail::range_or_tuple_size(field);
    }

    //!\brief Overload for n_allele.
    void set_core_n_allele(auto & field)
    {
        record_core.n_allele = 1; // for REF

        using field_t = decltype(field);

        // assuming single string
        if constexpr (std::ranges::forward_range<field_t> &&
                      !std::ranges::forward_range<std::ranges::range_reference_t<field_t>>)
        {
            record_core.n_allele++;
        }
        else // range or tuple of sequence
        {
            record_core.n_allele += detail::range_or_tuple_size(field);
        }
    }

    //!\brief Overload for n_fmt.
    void set_core_n_fmt(auto & field) { record_core.n_fmt = detail::range_or_tuple_size(field); }
    //!\}

    /*!\name Field writers
     * \{
     */

    //!\brief Overload for ID.
    void write_field(meta::vtag_t<detail::field::id> /**/, auto && field)
    {
        static_assert(detail::type_2_bcf_type_descriptor<std::remove_cvref_t<decltype(field)>> ==
                        detail::bcf_type_descriptor::char8,
                      "ID field must be provided as string.");
        write_typed_data(field);
    }

    //!\brief Overload for REF.
    void write_field(meta::vtag_t<detail::field::ref> /**/, auto && field)
    {
        static_assert(detail::type_2_bcf_type_descriptor<std::remove_cvref_t<decltype(field)>> ==
                        detail::bcf_type_descriptor::char8,
                      "REF field must be provided as string.");
        write_typed_data(field);
    }

    //!\brief Overload for ALT (single argument).
    void write_field(meta::vtag_t<detail::field::alt> /**/, auto && field)
    {
        static_assert(detail::type_2_bcf_type_descriptor<std::remove_cvref_t<decltype(field)>> ==
                        detail::bcf_type_descriptor::char8,
                      "ALT field must be provided as (range of) string(s).");
        write_typed_data(field);
    }

    //!\brief Overload for ALT that is range-of-range.
    template <std::ranges::input_range rng_t>
        requires std::ranges::forward_range<std::ranges::range_reference_t<rng_t>>
    void write_field(meta::vtag_t<detail::field::alt> /**/, rng_t && range)
    {
        for (auto && elem : range)
            write_field(meta::vtag<detail::field::alt>, elem);
    }

    //!\brief Overload for FILTER; handles vector of numeric IDs and vector of IDX
    template <std::ranges::input_range rng_t>
        requires(!std::same_as<std::ranges::range_value_t<rng_t>, char>)
    void write_field(meta::vtag_t<detail::field::filter> /**/, rng_t && range)
    {
        write_type_descriptor(idx_desc, std::ranges::distance(range));

        if constexpr (std::integral<std::ranges::range_value_t<rng_t>>)
        {
            write_range_impl(range, idx_desc);
        }
        else
        {
            auto text_id_to_idx = [this](std::string_view const text_id)
            {
                auto it = header->string_to_filter_pos().find(text_id);

                if (it == header->string_to_filter_pos().end())
                    error("The filter '", text_id, "' is not present in the header.");

                return header->filters[it->second].idx;
            };

            write_range_impl(range | std::views::transform(text_id_to_idx), idx_desc);
        }
    }

    //!\brief Overload for FILTER; single string or single IDX
    void write_field(meta::vtag_t<detail::field::filter> /**/, auto && field)
    {
        std::span<std::remove_reference_t<decltype(field)>> s{&field, 1};
        write_field(meta::vtag<detail::field::filter>, s); // delegate to previous overload
    }

    //!\brief Deduce descriptor from parameter type and optionally compress (integers) and verify with header.
    detail::bcf_type_descriptor get_desc(auto & param, var_io::header::info_t const & hdr_entry)
    {
        using param_t                                = std::remove_cvref_t<decltype(param)>;
        constexpr detail::bcf_type_descriptor c_desc = detail::type_2_bcf_type_descriptor<param_t>;
        detail::bcf_type_descriptor           desc   = c_desc;

        if constexpr (c_desc == detail::bcf_type_descriptor::int8 || c_desc == detail::bcf_type_descriptor::int16 ||
                      c_desc == detail::bcf_type_descriptor::int32)
        {
            // explicit integer width given in header
            if (hdr_entry.other_fields.find("IntegerBits") != hdr_entry.other_fields.end())
            {
                desc = detail::value_type_id_2_type_descriptor(hdr_entry.type_id);
                if (!detail::type_descriptor_is_int(desc)) // ignore header value if it isn't intX
                    desc = c_desc;
            }
            else if constexpr (c_desc == detail::bcf_type_descriptor::int16 ||
                               c_desc == detail::bcf_type_descriptor::int32)
            {
                if (compress_integers)
                    desc = detail::smallest_int_desc(param);
            }
        }

        if (verify_header_types)
        {
            detail::bcf_type_descriptor header_desc = detail::value_type_id_2_type_descriptor(hdr_entry.type_id);
            if (desc != header_desc || !detail::type_descriptor_is_int(desc) ||
                !detail::type_descriptor_is_int(header_desc))
            {
                error("The type of field ",
                      hdr_entry.id,
                      " set in the header is different from the current record's "
                      "data.");
            }
        }

        return desc;
    }

    //!\brief Write a single info element.
    void write_info_element(auto & info_element)
    {
        auto & [id, value] = info_element;

        using id_t    = decltype(id);
        using value_t = decltype(value);

        /* ID */
        int32_t idx = 0;
        if constexpr (std::integral<id_t>)
        {
            assert((int64_t)id <= (int64_t)header->max_idx());
            idx = id;
        }
        else
        {
            auto it = header->string_to_info_pos().find(id);
            if (it == header->string_to_info_pos().end())
                error("The info '", id, "' is not present in the header.");

            idx = header->infos[it->second].idx;
        }
        write_typed_data(idx, idx_desc);

        var_io::header::info_t const & info = header->infos.at(header->idx_to_info_pos().at(idx));

        /* VALUE */
        if constexpr (detail::is_info_element_value_type<value_t>)
        {
            auto func = [&](auto & param) { write_typed_data(param, get_desc(param, info)); };
            std::visit(func, value);
        }
        else
        {
            write_typed_data(value, get_desc(value, info));
        }
    }

    //!\brief Overload for INFO; range of pairs.
    template <std::ranges::input_range rng_t>
        requires(detail::info_element_writer_concept<std::ranges::range_reference_t<rng_t>>)
    void write_field(meta::vtag_t<detail::field::info> /**/, rng_t && range)
    {
        for (auto & info_element : range)
            write_info_element(info_element);
    }

    //!\brief Overload for INFO; tuple of pairs.
    template <typename... elem_ts>
        requires(detail::info_element_writer_concept<elem_ts> &&...)
    void write_field(meta::vtag_t<detail::field::info> /**/, std::tuple<elem_ts...> & tup) // TODO add const version
    {
        auto func = [&](auto &... field) { (write_info_element(field), ...); };
        std::apply(func, tup);
    }

    /*!\brief Get the maximum allele count and maximum allele value in the set of GT strings.
     * \param range_of_strings The GT strings, e.g. "0/0", "1/1" ...
     * \returns A pair of "maximum allele count" and "maximum allele value".
     */
    std::pair<size_t, int32_t> get_GT_maxs(auto & range_of_strings)
    {
        size_t  max_alleles    = 0;
        int32_t max_allele_val = 0;

        for (auto && s_tmp : range_of_strings)
        {
            std::string_view s = detail::to_string_view(s_tmp);

            int32_t buf       = 0;
            size_t  n_alleles = 0;

            // TODO perform sanity check on string?

            // TODO use views::eager_split once that can handle predicates
            for (size_t i = 0, j = 0; i <= s.size(); ++i)
            {
                if (i == s.size() || (is_char<'/'> || is_char<'|'>)(s[i]))
                {
                    std::string_view substr = s.substr(j, i - j);

                    if (substr != ".")
                    {
                        detail::string_to_number(substr, buf);
                        max_allele_val = std::max(max_allele_val, buf);
                    }

                    j = i + 1;
                    ++n_alleles;
                }
            }

            max_alleles = std::max(max_alleles, n_alleles);
        }

        return {max_alleles, max_allele_val};
    }

    //!\brief Write a single Genotypes element.
    void write_genotypes_element(auto & genotype)
    {
        auto & [id, value] = genotype;

        using id_t    = decltype(id);
        using value_t = decltype(value);

        /* ID */
        int32_t idx = 0;
        if constexpr (std::integral<id_t>)
        {
            assert((int64_t)id <= (int64_t)header->max_idx());
            idx = id;
        }
        else
        {
            if (auto it = header->string_to_format_pos().find(id); it == header->string_to_format_pos().end())
                error("The genotype '", id, "' is not present in the header.");
            else
                idx = header->formats[it->second].idx;
        }
        write_typed_data(idx, idx_desc);

        var_io::header::format_t const & format = header->formats.at(header->idx_to_format_pos().at(idx));

        /* value */
        auto func = [&](auto & values)
        {
            using values_t = decltype(values);
            static_assert(std::ranges::range<values_t>, "This function handles only ranges; write_genotypes_element.");
            using value_t = std::ranges::range_value_t<values_t>;

            if (std::ranges::size(values) > record_core.n_sample)
            {
                error("There are ",
                      std::ranges::size(values),
                      " values in the genotype vector "
                      "for field ",
                      format.id,
                      " which is more than the number of samples.");
            }

            detail::bcf_type_descriptor desc = get_desc(values, format);

            // if there are no values, we can write the missing field and need no padding values at all
            if (std::ranges::size(values) == 0)
            {
                write_type_descriptor(desc, 0);
                return;
            }

            if constexpr (!std::ranges::range<value_t>) // case: one value per sample
            {
                write_type_descriptor1(desc); // size refers to size per-sample, which is one
                write_range_impl(values, desc);

                // per-sample padding
                write_range_padding(record_core.n_sample - std::ranges::size(values), desc, false);
            }
            // case: multiple values per sample or one string-per-sample
            else if constexpr (!std::ranges::range<std::ranges::range_reference_t<value_t>>)
            {
                // this looks like "one-string" but is actually a case with custom encoding
                if (format.id == "GT")
                {
                    if constexpr (detail::char_range_or_cstring<value_t>)
                    {
                        assert(desc == detail::bcf_type_descriptor::char8);

                        if (record_core.n_sample != std::ranges::size(values))
                            error("If the GT field is present, a value must be given for every sample.");

                        auto [max_alleles, max_allele_val] = get_GT_maxs(values);

                        // computation of binary value is "(allele_value + 1) << 1 | phased"
                        // we loose one bit because of +1 and one bit because of the shift and one because of signed:
                        if (max_allele_val <= 5)
                            desc = detail::bcf_type_descriptor::int8;
                        else if (max_allele_val <= 13)
                            desc = detail::bcf_type_descriptor::int16;
                        else
                            desc = detail::bcf_type_descriptor::int32;

                        auto fun = [&](auto && val, detail::bcf_type_descriptor const d)
                        {
                            size_t n_alleles = write_GT_impl(detail::to_string_view(val), d);

                            // per-value padding
                            write_range_padding(max_alleles - n_alleles, d, false);
                        };

                        // size refers to size per-sample
                        write_type_descriptor(desc, max_alleles);
                        for (auto && rng : values)
                            fun(rng, desc);

                        // no padding (GT required for all samples)
                    }
                    else
                    {
                        error("GT field must be provided as range of strings.");
                    }
                }
                else // the regular case of multiple values or one string (per sample)
                {
                    size_t max_length = std::ranges::size(*std::ranges::max_element(values, {}, std::ranges::size));

                    write_type_descriptor(desc, max_length); // size refers to size per-sample

                    for (auto && rng : values)
                    {
                        write_range_impl(rng, desc);
                        size_t const s = std::ranges::size(rng);
                        // per-value padding; we only insert missing value if this range is empty
                        write_range_padding(max_length - s, desc, s == 0);
                    }

                    // per-sample padding
                    for (size_t i = 0; i < record_core.n_sample - std::ranges::size(values); ++i)
                        write_range_padding(max_length, desc, true);
                }
            }
            else // this is the case of multiple strings per sample
            {
                if constexpr (std::same_as<char, bio::ranges::range_innermost_value_t<value_t>>)
                {
                    assert(desc == detail::bcf_type_descriptor::char8);

                    string_size_buffer.clear();
                    size_t max_length = 0;
                    for (auto && range_of_strings : values)
                    {
                        size_t sum = 0;

                        for (auto && string : range_of_strings)
                            sum += std::ranges::size(string);

                        // add one comma per string after the first
                        sum += std::ranges::size(range_of_strings) ? std::ranges::size(range_of_strings) - 1 : 0;
                        string_size_buffer.push_back(sum);

                        max_length = std::max(max_length, sum);
                    }

                    write_type_descriptor(desc, max_length); // size refers to size per-sample

                    size_t i = 0;
                    for (auto && range_of_strings : values)
                    {
                        // this choses overload that intersperses ','
                        write_range_impl(range_of_strings, desc);
                        // per-value padding
                        write_range_padding(max_length - string_size_buffer[i], desc, string_size_buffer[i] == 0);
                        ++i;
                    }

                    // TODO instead of storing the cumulative sizes, we could also track the bytes written

                    // per-sample padding
                    for (size_t i = 0; i < record_core.n_sample - std::ranges::size(values); ++i)
                        write_range_padding(max_length, desc, true);
                }
                else
                {
                    // this case is never reached. The constexpr prevents needless instantiations though.
                    static_assert(std::same_as<char, bio::ranges::range_innermost_value_t<value_t>>,
                                  "Cannot handle range-of-range-of-range unless the alphabet is char.");
                }
            }
        };

        if constexpr (detail::is_genotype_element_value_type<value_t>)
            std::visit(func, value);
        else
            func(value);
    }

    //!\brief Overload for GENOTYPES.
    template <std::ranges::forward_range range_t>
        requires(detail::genotype_writer_concept<std::ranges::range_reference_t<range_t>>)
    void write_field(meta::vtag_t<detail::field::genotypes> /**/, range_t && range)
    {
        for (auto && genotype : range)
            write_genotypes_element(genotype);
    }

    //!\brief Overload for GENOTYPES; tuple of pairs.
    template <typename... elem_ts>
        requires(detail::genotype_writer_concept<elem_ts> &&...)
    void write_field(meta::vtag_t<detail::field::genotypes> /**/,
                     std::tuple<elem_ts...> & tup) // TODO add const version
    {
        auto func = [&](auto &... field) { (write_genotypes_element(field), ...); };
        std::apply(func, tup);
    }
    //!\}

    //!\brief Write the header.
    void write_header()
    {
        if (header == nullptr)
            throw missing_header_error{"You need to call set_header() on the writer/format."};

        std::string vcf_header = header->to_plaintext();

        // BINARY HEADER
        it->write_range("BCF");
        it->write_as_binary(uint8_t{2});
        it->write_as_binary(uint8_t{2});
        uint32_t l_text = static_cast<uint32_t>(vcf_header.size()) + 1;
        it->write_as_binary(l_text);

        it->write_range(vcf_header);
        it                      = '\0';
        header_has_been_written = true;

        /* compute the smallest int type for the dictionaries */
        idx_desc = detail::smallest_int_desc(header->max_idx());
    }

    //!\brief Write the record (supports const and non-const lvalue ref).
    void write_record_impl(auto & record)
    {
        using record_t = std::remove_cvref_t<decltype(record)>;

        if (!header_has_been_written)
        {
            if (header == nullptr)
            {
                if (var_io::header const * ptr = record._private.header_ptr; ptr != nullptr)
                    set_header(*ptr);
            }

            write_header();
        }

        static_assert(meta::different_from<typename record_t::chrom_t, ignore_t>,
                      "The record must contain the CHROM field.");
        static_assert(meta::different_from<typename record_t::pos_t, ignore_t>,
                      "The record must contain the POS field.");
        static_assert(meta::different_from<typename record_t::ref_t, ignore_t>,
                      "The record must contain the REF field.");

        /* IMPLEMENTATION NOTE:
         * A problem when writing BCF is that the first fields of the record hold the size of the record
         * (or parts thereof). These are however not known to us, yet, because we haven't transformed the input
         * data to the respective binary representation, yet.
         * We therefor use an intermediate buffer stream implemented as an ostringstream and write placeholder
         * values initially. Later, after the full record is written to the buffer stream, we come back and
         * replace the placeholders with the actual values. Since we are using a stringstream, we know that
         * the full record is inside the buffer (in contrast to using the regular stream which might have been
         * flushed to disk already).
         * The stringstream's pbase() always points to the beginning of the string-buffer and we store the
         * current record's begin in a separate variable.
         * The contents of the stringstream are flushed to the actual stream periodically, but only after
         * it contains at least one full record. Flushing to the actual stream happens via sputn() so the
         * hope is that for sufficiently large buffers the write is not buffered again.
         *
         * Alternative implementations like SEEKing back in the fstream where considered, but seeking is not
         * guaranteed to work on arbitrary streams and is particularly inefficient on gzipped streams.
         */

        // begin position of this record in the output stream's buffer
        this_record_offset = streambuf_exposer->pptr() - streambuf_exposer->pbase();

        it->write_as_binary(uint32_t{}); // l_shared
        uint32_t l_shared_tmp = 0;
        it->write_as_binary(uint32_t{}); // l_indiv
        uint32_t l_indiv_tmp = 0;

        /* this prepares the record_core and doesn't write anything til the end */
        record_core = {}; // reset this
        set_core_chrom(record.chrom);
        set_core_pos(record.pos);

        set_core_rlen(record.ref);

        if constexpr (meta::different_from<typename record_t::qual_t, ignore_t>)
            set_core_qual(record.qual);

        if constexpr (meta::different_from<typename record_t::info_t, ignore_t>)
            set_core_n_info(record.info);

        if constexpr (meta::different_from<typename record_t::alt_t, ignore_t>)
            set_core_n_allele(record.alt);
        else
            record_core.n_allele = 1; // the REF allele

        record_core.n_sample = header->column_labels.size() > 9 ? header->column_labels.size() - 9 : 0;

        if constexpr (meta::different_from<typename record_t::genotypes_t, ignore_t>)
            set_core_n_fmt(record.genotypes);

        // write record core
        it->write_as_binary(record_core);

        /* After this point, the order of writers is important! */

        if constexpr (meta::different_from<typename record_t::id_t, ignore_t>)
            write_field(meta::vtag<detail::field::id>, record.id);
        else
            write_field(meta::vtag<detail::field::id>, std::string_view{});

        write_field(meta::vtag<detail::field::ref>, record.ref);

        if constexpr (meta::different_from<typename record_t::alt_t, ignore_t>)
            write_field(meta::vtag<detail::field::alt>, record.alt);
        else
            write_field(meta::vtag<detail::field::alt>, std::span<std::string_view>{});

        if constexpr (meta::different_from<typename record_t::filter_t, ignore_t>)
            write_field(meta::vtag<detail::field::filter>, record.filter);
        else
            write_field(meta::vtag<detail::field::filter>, std::span<std::string_view>{});

        if constexpr (meta::different_from<typename record_t::info_t, ignore_t>)
            write_field(meta::vtag<detail::field::info>, record.info);
        else
            write_field(meta::vtag<detail::field::info>, std::span<var_io::info_element<>>{});

        assert(streambuf_exposer->pptr() - streambuf_exposer->pbase() - this_record_offset > 0);
        //              (position in the buffer                              ) - (where this record starts)
        l_shared_tmp = streambuf_exposer->pptr() - streambuf_exposer->pbase() - this_record_offset;

        if (header->column_labels.size() > 8)
        {
            if constexpr (meta::different_from<typename record_t::genotypes_t, ignore_t>)
                write_field(meta::vtag<detail::field::genotypes>, record.genotypes);
            else
                write_field(meta::vtag<detail::field::genotypes>, std::span<var_io::genotype_element<>>{});
        }

        l_indiv_tmp = streambuf_exposer->pptr() - streambuf_exposer->pbase() - this_record_offset - l_shared_tmp;

        // now we get a pointer to where l_shared and l_indiv are stored in the buffer and write them
        uint32_t * iptr = reinterpret_cast<uint32_t *>(streambuf_exposer->pbase() + this_record_offset);
        *iptr           = l_shared_tmp - 8; // size of l_shared and l_indiv not counted
        ++iptr;
        *iptr = l_indiv_tmp;

        // possibly turn this into an option or take from the transparent_istream options
        ptrdiff_t min_write_size = 10 * 1024 * 1024; // 10MB
        if (streambuf_exposer->pptr() - streambuf_exposer->pbase() > min_write_size)
        {
            // write data from the buffer_stream into the actual stream
            // this *should* be unbuffered for large min_write_sizes
            stream->rdbuf()->sputn(streambuf_exposer->pbase(), streambuf_exposer->pptr() - streambuf_exposer->pbase());

            // reset the buffer_stream(first argument sets pbase AND pptr):
            streambuf_exposer->setp(streambuf_exposer->pbase(), streambuf_exposer->epptr());
        }

        ++record_no;
    }

public:
    /*!\name Constructors, destructor and assignment.
     * \{
     */
    format_output_handler()                                          = delete; //!< Deleted.
    format_output_handler(format_output_handler const &)             = delete; //!< Deleted.
    format_output_handler & operator=(format_output_handler const &) = delete; //!< Deleted.

    //!\brief Move construction.
    format_output_handler(format_output_handler && rhs) { move_impl(std::move(rhs)); }

    //!\brief Move assignment.
    format_output_handler & operator=(format_output_handler && rhs)
    {
        move_impl(std::move(rhs));
        return *this;
    }

    //!\brief Helper function for move-constructor and move-assignment.
    void move_impl(format_output_handler && rhs)
    {
        std::swap(stream, rhs.stream);
        std::swap(header_has_been_written, rhs.header_has_been_written);
        std::swap(header, rhs.header);
        std::swap(idx_desc, rhs.idx_desc);
        std::swap(record_no, rhs.record_no);
        std::swap(record_core, rhs.record_core);
        std::swap(this_record_offset, rhs.this_record_offset);
        std::swap(buffer_stream, rhs.buffer_stream);
        std::swap(compress_integers, rhs.compress_integers);
        std::swap(verify_header_types, rhs.verify_header_types);
        std::swap(string_size_buffer, rhs.string_size_buffer);

        streambuf_exposer     = reinterpret_cast<detail::stream_buffer_exposer<char> *>(buffer_stream.rdbuf());
        rhs.streambuf_exposer = nullptr;

        it = detail::fast_ostreambuf_iterator<char>{buffer_stream};
        // rhs.it does not have to be handled

        moved_from     = false;
        rhs.moved_from = true;
    }

    /*!\brief Construct with an options object.
     * \param[in,out] str The output stream.
     * \param[in] options An object with options for the output handler.
     * \details
     *
     * The options argument is typically bio::io::var_io::writer_options, but any object with a subset of similarly
     * named members is also accepted. See bio::io::format_output_handler<bcf> for the supported options and defaults.
     */
    format_output_handler(std::ostream & str, auto const & options) : stream{&str}
    {
        // extract options
        if constexpr (requires { (bool)options.compress_integers; })
            compress_integers = options.compress_integers;
        if constexpr (requires { (bool)options.verify_header_types; })
            verify_header_types = options.verify_header_types;
    }

    //!\brief Construct with only an output stream.
    format_output_handler(std::ostream & str) : format_output_handler(str, 1) {}

    //!\brief Destruction with flushing and cleanup.
    ~format_output_handler() noexcept(false)
    {
        // never throw if the stack is unwinding
        if (std::uncaught_exceptions() > 0)
            return;

        // no cleanup is needed if this object is in moved-from-state
        if (moved_from)
            return;

        // if no records were written, the header also wasn't written
        if (!header_has_been_written)
            write_header();

        // flush stringstream to regular stream so buffer content is written to disk
        if (streambuf_exposer->pptr() - streambuf_exposer->pbase() > 0)
            stream->rdbuf()->sputn(streambuf_exposer->pbase(), streambuf_exposer->pptr() - streambuf_exposer->pbase());
    }
    //!\}

    //!\brief Get the header.
    var_io::header const & get_header() const
    {
        if (header == nullptr)
            throw missing_header_error{"Attempting to read header, but no header was set."};

        return *header;
    }

    //!\brief Set the header.
    void set_header(var_io::header const & hdr)
    {
        // TODO verify that all header entries have IDX
        header = {&hdr, [](var_io::header const *) {}};
    }
    //!\overload
    void set_header(var_io::header const && hdr)
    {
        // TODO verify that all header entries have IDX
        header = {new var_io::header(std::move(hdr)), [](var_io::header const * ptr) { delete ptr; }};
    }
    //!\overload
    void set_header(var_io::header & hdr)
    {
        hdr.add_missing();
        set_header(std::as_const(hdr));
    }
    //!\overload
    void set_header(var_io::header && hdr)
    {
        hdr.add_missing();
        set_header(std::move(std::as_const(hdr)));
    }

    //!\brief Write the record.
    template <typename... field_types>
    void write_record(var_io::record<field_types...> const & record)
    {
        write_record_impl(record);
    }

    //!\overload
    template <typename... field_types>
    void write_record(var_io::record<field_types...> & record)
    {
        write_record_impl(record);
    }
};

} // namespace bio::io
