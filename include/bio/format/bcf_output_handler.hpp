// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::format_output_handler<bcf>.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bit>

#include <bio/detail/magic_get.hpp>
#include <bio/detail/misc.hpp>
#include <bio/format/format_output_handler.hpp>
#include <bio/format/bcf.hpp>
#include <bio/record.hpp>
#include <bio/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/var_io/header.hpp>
#include <bio/var_io/misc.hpp>
#include <bio/var_io/writer_options.hpp>

namespace bio
{

/*!\brief Format output handler for the VCF format (bio::bcf).
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
class format_output_handler<bcf> : public format_output_handler_base<format_output_handler<bcf>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The base class.
    using base_t = format_output_handler_base<format_output_handler<bcf>>;
    //!\brief Befriend the base class so we can instantiate.
    friend base_t;

    //DO NOT INHERIT it from base!
//     using base_t::it;
    using base_t::stream;
    using base_t::write_field;
//     using base_t::write_field_aux;
    //!\}

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
    detail::bcf_type_descriptor header_dictionary_desc;

    size_t record_no = 0;

    //!\brief This is cached during raw-record reading already.
    detail::bcf_record_core record_core;

    //!\brief An intermediate stream that is guaranteed to always hold the complete current record in memory.
    std::ostringstream buffer_stream;

    detail::stream_buffer_exposer<char> * streambuf_exposer =
        reinterpret_cast<detail::stream_buffer_exposer<char> *>(buffer_stream.rdbuf());
    //!\brief The iterator on that stream.
    detail::fast_ostreambuf_iterator<char> it{buffer_stream};

    //!\brief Distance from pbase()
    size_t this_record_offset = 0;
    //!\}

    /*!\name Options
     * \{
     */
    //!\brief Try to use the smallest possible integer type (creates smaller files but is potentially slower).
    bool compress_integers = false;
    //!\}
    //!\}
    /*!\name Arbitrary helpers
     * \{
     */

    detail::bcf_type_descriptor smallest_int_desc(std::unsigned_integral auto const num)
    {
        // bits required to represent number (the +1 because signed integral has smaller range)
        switch (std::bit_ceil(std::bit_width(num) + 1))
        {
            case 128:
            case 64:
                error("Could not write number '", num, "'. Value out of range (only int32 supported).");
                return {};
            case 32:
                return detail::bcf_type_descriptor::int32;
            case 16:
                return detail::bcf_type_descriptor::int16;
            default:
                return detail::bcf_type_descriptor::int8;
        }
    }

    detail::bcf_type_descriptor smallest_int_desc(std::signed_integral auto const num)
    {
        return smallest_int_desc(static_cast<uint64_t>(std::abs(num)));
    }

    detail::bcf_type_descriptor smallest_int_desc(std::ranges::forward_range auto && range)
    {
        auto max = std::ranges::max(range);
//         //TODO check if this is faster:
//         val_t max = std::numeric_limits<val_t>::lowest();
//         for (val_t elem : data)
//         {
//             if (elem > max)
//             {
//                 max = elem;
//                 desc = smallest_int_desc(elem);
//                 if (desc == detail::bcf_type_descriptor::int32) // this is max(type_descriptor)
//                     break;
//             }
//         }
        return smallest_int_desc(max);
    }

    template <std::ranges::forward_range rng_t>
        requires std::ranges::range<std::ranges::range_reference_t<rng_t>>
    detail::bcf_type_descriptor smallest_int_desc(rng_t & range)
    {
        return smallest_int_desc(range | std::views::join);
    }

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
            write_typed_data(size, true /*always shrink this number*/);
        }
    }

    //!\brief Overload that is hardcoded for size==1
    void write_type_descriptor1(detail::bcf_type_descriptor const desc)
    {
        uint8_t byte = static_cast<uint8_t>(1) << 4;
        byte |= static_cast<uint8_t>(desc);
        it->write_as_binary(byte);
    }

    void write_single_impl(seqan3::arithmetic auto num, detail::bcf_type_descriptor const desc)
    {
        write_type_descriptor1(desc);
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

    void write_range_impl(auto & range, detail::bcf_type_descriptor const desc)
    {
        using data_t = decltype(range);
        using elem_t = std::ranges::range_value_t<data_t>;

        auto dump_or_convert = [&] <typename target_t> (target_t)
        {
            if constexpr (std::same_as<target_t, elem_t> &&
                          std::ranges::contiguous_range<data_t> &&
                          std::ranges::sized_range<data_t>)
            {
                it->write_as_binary(range);
            }
            else if constexpr (detail::deliberate_alphabet<elem_t>)
            {
                it->write_range(range | seqan3::views::to_char);
            }
            else
            {
                //TODO this is probably not the most efficient implementation
                for (elem_t elem : range)
                {
                    target_t buf = static_cast<target_t>(elem);
                    std::string_view v{reinterpret_cast<char const *>(&buf), sizeof(buf)};
                    it->write_range(v);
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

    //!\brief Write single element.
    template <typename elem_t>
        requires (seqan3::alphabet<elem_t> || seqan3::arithmetic<elem_t>)
    void write_typed_data(elem_t const num, bool shrink_int)
    {
        detail::bcf_type_descriptor desc = detail::type_2_bcf_type_descriptor<elem_t>;
        if constexpr (std::integral<elem_t>)
            if (shrink_int)
                desc = smallest_int_desc(num);

        write_single_impl(num, desc);
    }

    void write_typed_data(detail::char_range auto && data, bool)
    {
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

    void write_typed_data(char const * const cstring, bool)
    {
        write_typed_data(std::string_view{cstring}, false);
    }

    void write_typed_data(std::ranges::forward_range auto && data, bool shrink_int)
    {
        using data_t = decltype(data);
        using elem_t = std::ranges::range_value_t<data_t>;

        detail::bcf_type_descriptor desc =detail::type_2_bcf_type_descriptor<elem_t>;
        if constexpr (std::integral<elem_t>)
            if (shrink_int)
                desc = smallest_int_desc(data);

        write_type_descriptor(desc, std::ranges::distance(data));

        write_range_impl(data, desc);
    }

    template <std::ranges::forward_range data_t>
        requires detail::char_range<std::ranges::range_value_t<data_t>>
    void write_typed_data(data_t && data, bool)
    {
        error("implement me");
    }

    void write_typed_data(bool)
    {
        error("You should never end up here.");
    }

    //!\brief This overload adds a default argument for integer compression based on the class member.
    void write_typed_data(auto && num)
    {
        write_typed_data(num, compress_integers);
    }

    //!\brief Write chracter ranges.
    void write_field_aux(detail::char_range auto & range) { write_typed_data(range); }

    //!\brief Write alphabet ranges.
    template <std::ranges::forward_range rng_t>
        requires(detail::deliberate_alphabet<std::ranges::range_reference_t<rng_t>>)
    void write_field_aux(rng_t & range) { write_typed_data(range | seqan3::views::to_char); }

    //!\brief Write CStrings.
    void write_field_aux(char const * const cstr) { write_typed_data(std::string_view{cstr}); }

    //!\brief Write numbers.
    void write_field_aux(seqan3::arithmetic auto number) { write_types_single_data(number); }

    /*!\name Core record setters
     * \brief These set up the core record.
     * \{
     */
    //!\brief Overload for CHROM and numeric IDs.
    void set_core_chrom(std::integral auto & field)
    {
        record_core.chrom = static_cast<int32_t>(field);
    }

    //!\brief Overload for CHROM and text IDs.
    void set_core_chrom(detail::char_range auto && field)
    {
        std::string_view const buf = field; // TODO make this more generic
        if (auto it = header->string_to_contig_pos().find(buf); it == header->string_to_contig_pos().end())
            error("The contig '", field, "' is not present in the header.");
        else
            record_core.chrom = header->contigs[it->second].idx;
    }

    //!\overload
    void set_core_chrom(char const * const field)
    {
        set_core_chrom(std::string_view{field});
    }

    //!\brief Overload for POS.
    void set_core_pos(std::integral auto & field)
    {
        record_core.pos = static_cast<int32_t>(field);
    }

    //!\brief Overload for rlen.
    void set_core_rlen(auto && field)
    {
        record_core.rlen = static_cast<int32_t>(std::ranges::distance(field));
    }

    //!\overload
    void set_core_rlen(char const * const field)
    {
        set_core_rlen(std::string_view{field});
    }

    //!\brief Overload for QUAL.
    void set_core_qual(seqan3::arithmetic auto & field)
    {
        record_core.qual = static_cast<float>(field);
    }

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
    void set_core_n_fmt(auto & field)
    {
        if constexpr (detail::genotypes_vcf_style_writer_concept<decltype(field)>)
            record_core.n_fmt = std::ranges::distance(detail::get_first(field));
        else
            record_core.n_fmt = detail::range_or_tuple_size(field);
    }
    //!\}

    /*!\name Field writers
     * \{
     */

    // ID, REF, (one-element) ALT handled by defaults

    //!\brief Overload for ALT that is range-of-range (single-range ALT is handled by defaults).
    template <std::ranges::input_range rng_t>
        requires std::ranges::forward_range<std::ranges::range_reference_t<rng_t>> // TOOD and requires
                                                                                   // write_field_aux(value)
    void write_field(vtag_t<field::alt> /**/, rng_t & range)
    {
        for (auto && elem : range)
            write_field_aux(elem);
    }

    //!\brief Overload for FILTER; handles vector of numeric IDs and vector of IDX
    template <std::ranges::input_range rng_t>
        requires(!std::same_as<std::ranges::range_value_t<rng_t>, char>)
    void write_field(vtag_t<field::filter> /**/, rng_t & range)
    {
        if constexpr (std::integral<std::ranges::range_value_t<rng_t>>)
        {
            write_typed_data(range);
        }
        else
        {
            auto text_id_to_idx = [this] (std::string_view const text_id)
            {
                auto it = header->string_to_filter_pos().find(text_id);

                if (it == header->string_to_filter_pos().end())
                    error("The filter '", text_id, "' is not present in the header.");

                return header->filters[it->second].idx;
            };

            // TODO we are not using header_dictionary_desc here
            write_typed_data(range | std::views::transform(text_id_to_idx));
        }
    }

    //!\brief Overload for FILTER; single string or single IDX
    void write_field(vtag_t<field::filter> /**/, auto & field)
    {
//         std::span<decltype(field)> s{&field, 1};
        write_field(vtag<field::filter>, std::initializer_list{field}); // delegate to previous overload
    }


    void write_info_element(auto & info_element)
    {
        auto & [ id, value ] = info_element;

        using id_t = decltype(id);
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
        write_single_impl(idx, header_dictionary_desc);

        /* VALUE */
        if constexpr (detail::is_dynamic_type<value_t>)
        {
            auto func = [&] (auto & param) { write_typed_data(param); };
            std::visit(func, value);
        }
        else
        {
            write_typed_data(value);
        }
    }

    //!\brief Overload for INFO; range of pairs.
    template <std::ranges::input_range rng_t>
        requires(detail::info_element_writer_concept<std::ranges::range_reference_t<rng_t>>)
    void write_field(vtag_t<field::info> /**/, rng_t & range)
    {
        for (auto & info_element : range)
            write_info_element(info_element);
    }

    //!\brief Overload for INFO; tuple of pairs.
    template <typename... elem_ts>
        requires(detail::info_element_writer_concept<elem_ts> &&...)
    void write_field(vtag_t<field::info> /**/, std::tuple<elem_ts...> & tup) // TODO add const version
    {
        auto func = [&](auto & ... field) { (write_info_element(field), ...); };
        std::apply(func, tup);
    }

    void write_genotypes_element(auto & genotype)
    {
        auto & [ id, value ] = genotype;

        using id_t = decltype(id);
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
            auto it = header->string_to_format_pos().find(id);
            if (it == header->string_to_format_pos().end())
                error("The genotype '", id, "' is not present in the header.");

            idx = header->formats[it->second].idx;
        }
        write_single_impl(idx, header_dictionary_desc);

        /* value */
        auto func = [&](auto & values)
        {
            using values_t = decltype(values);
            static_assert(std::ranges::range<values_t>, "This function handles only ranges; write_genotypes_element.");
            using value_t = std::ranges::range_value_t<values_t>;

            using alph_t = seqan3::range_innermost_value_t<values_t>;


            if (std::ranges::size(values) != 0 && std::ranges::size(values) != record_core.n_sample)
            {
                error("There are ", std::ranges::size(values), " values in the genotype vector, "
                      "but there must be exactly one per sample or none at all.");
            }

            detail::bcf_type_descriptor desc = detail::type_2_bcf_type_descriptor<alph_t>;
            if constexpr (std::integral<alph_t>)
                if (compress_integers)
                    desc = smallest_int_desc(values);

            if (std::ranges::size(values) == 0)
            {
                write_type_descriptor(desc, 0);
                return;
            }

            if constexpr (std::ranges::range<value_t>)
            {
                size_t max_length      = 1;
                bool   all_same_length = true;
                size_t last_length       = std::ranges::size(*std::ranges::begin(values));

                for (auto & rng : values)
                {
                    size_t s = std::ranges::size(rng);
                    if (s > max_length)
                        max_length = std::ranges::size(rng);
                    if (s != last_length)
                        all_same_length = false;
                    last_length = s;
                }

                write_type_descriptor(desc, max_length); // size refers to size per-sample

                if (all_same_length)
                {
                    for (auto & rng : values)
                        write_range_impl(rng, desc);
                }
                else
                {
                    //TODO pad individual range so they all have same length
                    error("implement me");
                }
            }
            else
            {
                write_type_descriptor1(desc); // size refers to size per-sample, which is one, because this isn't range
                write_range_impl(values, desc);
            }

        };

        if constexpr (detail::is_dynamic_vector_type<value_t>)
            std::visit(func, value);
        else
            func(value);
    }

    //!\brief Overload for GENOTYPES; genotypes_bcf_style.
    template <std::ranges::forward_range range_t>
        requires(detail::genotype_bcf_style_writer_concept<std::ranges::range_reference_t<range_t>>)
    void write_field(vtag_t<field::genotypes> /**/, range_t & range)
    {
        for (auto && genotype : range)
            write_genotypes_element(genotype);
    }

    //!\brief Overload for GENOTYPES; tuple of pairs.
    template <typename... elem_ts>
        requires(detail::genotype_bcf_style_writer_concept<elem_ts> &&...)
    void write_field(vtag_t<field::genotypes> /**/, std::tuple<elem_ts...> & tup) // TODO add const version
    {
        auto func = [&](auto & ... field) { (write_genotypes_element(field), ...); };
        std::apply(func, tup);
    }

    //TODO vcf-style
    //!\}

    //!\brief Write the record (supports const and non-const lvalue ref).
    void write_record_impl(auto & record)
    {
        using field_ids = typename std::remove_cvref_t<decltype(record)>::field_ids;

        if (!header_has_been_written)
        {
            if (header == nullptr)
            {
                bool set = false;

                if constexpr (field_ids::contains(field::_private))
                {
                    if (var_io::header const * ptr = get<field::_private>(record).header_ptr; ptr != nullptr)
                    {
                        set_header(*ptr);
                        set = true;
                    }
                }

                if (!set)
                {
                    throw std::runtime_error{
                      "You need to call set_header() on the writer/format before writing a "
                      "record."};
                }
            }

            std::string vcf_header = header->to_plaintext();

            // BINARY HEADER
            it->write_range("BCF");
            it->write_as_binary(uint8_t{2});
            it->write_as_binary(uint8_t{2});
            uint32_t l_text = static_cast<uint32_t>(vcf_header.size());
            it->write_as_binary(l_text);

            it->write_range(vcf_header);
            it = '\0';
            header_has_been_written = true;

            /* compute the smallest int type for the dictionaries */
            header_dictionary_desc = smallest_int_desc(header->max_idx());
        }

        static_assert(field_ids::contains(field::chrom), "The record must contain the CHROM field.");
        static_assert(field_ids::contains(field::pos), "The record must contain the POS field.");
        static_assert(field_ids::contains(field::ref), "The record must contain the REF field.");


        /* IMPLEMENTATION NOTE:
         * A problem when writing BCF is that the first fields of the record hold the size of the record
         * (or parts thereof). These are however not known to us, yet, because we haven't transformed the input
         * data to the respective binary representation, yet.
         * We therefor use an intermediate buffer stream implemented as an ostringstream and write placeholder
         * values initially. Later, after the full record is written to the buffer stream, we come back and
         * replace the placeholders with the actual values. Since we are using a stringstream, we know that
         * the full record is inside the buffer (in contrast to using the regular stream which might have been
         * flushed to disk already).
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
        set_core_chrom(get<field::chrom>(record));
        set_core_pos(get<field::pos>(record));

        set_core_rlen(get<field::ref>(record));

        if constexpr (field_ids::contains(field::qual))
            set_core_qual(get<field::qual>(record));

        if constexpr (field_ids::contains(field::info))
            set_core_n_info(get<field::info>(record));

        if constexpr (field_ids::contains(field::alt))
            set_core_n_allele(get<field::alt>(record));

        record_core.n_sample = header->column_labels.size() > 9 ? header->column_labels.size() - 9 : 0;

        if constexpr (field_ids::contains(field::genotypes))
            set_core_n_fmt(get<field::genotypes>(record));

        // write record core
        it->write_as_binary(record_core);

        /* After this point, the order of writers is important! */

        if constexpr (field_ids::contains(field::id))
            write_field(vtag<field::id>, get<field::id>(record));
        else
            write_field(vtag<field::id>, var_io::missing_value<std::string_view>);

        write_field(vtag<field::ref>, get<field::ref>(record));

        if constexpr (field_ids::contains(field::alt))
            write_field(vtag<field::alt>, get<field::alt>(record));
        else
            write_field(vtag<field::alt>, var_io::missing_value<std::string_view>);

        if constexpr (field_ids::contains(field::filter))
            write_field(vtag<field::filter>, get<field::filter>(record));
        else
            write_field(vtag<field::filter>, var_io::missing_value<std::string_view>);

        if constexpr (field_ids::contains(field::info))
            write_field(vtag<field::info>, get<field::info>(record));
        else
            write_field(vtag<field::info>, var_io::missing_value<std::string_view>);

        assert(streambuf_exposer->pptr() - streambuf_exposer->pbase() - this_record_offset > 0);
        l_shared_tmp = streambuf_exposer->pptr() - streambuf_exposer->pbase() - this_record_offset;

        if (header->column_labels.size() > 8)
        {
            if constexpr (field_ids::contains(field::genotypes))
                write_field(vtag<field::genotypes>, get<field::genotypes>(record));
            else
                write_field(vtag<field::genotypes>, std::tuple<>{});
        }

        l_indiv_tmp = streambuf_exposer->pptr() - streambuf_exposer->pbase() - this_record_offset - l_shared_tmp;

        // now we get a pointer to where l_shared and l_indiv are stored in the buffer and write them
        uint32_t * iptr = reinterpret_cast<uint32_t*>(streambuf_exposer->pbase() + this_record_offset);
        *iptr = l_shared_tmp;
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
     * \brief These are all private to prevent wrong instantiation.
     * \{
     */
    format_output_handler()                              = default;            //!< Defaulted.
    format_output_handler(format_output_handler const &) = delete;             //!< Deleted.
    format_output_handler(format_output_handler &&)      = default;            //!< Defaulted.
    ~format_output_handler()                             = default;            //!< Defaulted.
    format_output_handler & operator=(format_output_handler const &) = delete; //!< Deleted.
    format_output_handler & operator=(format_output_handler &&) = default;     //!< Defaulted.

    /*!\brief Construct with an options object.
     * \param[in,out] str The output stream.
     * \param[in] options An object with options for the output handler.
     * \details
     *
     * The options argument is typically bio::var_io::writer_options, but any object with a subset of similarly named
     * members is also accepted. See bio::format_output_handler<bcf> for the supported options and defaults.
     */
    format_output_handler(std::ostream & str, auto const & options) : base_t{str}
    {
        // extract options
        if constexpr (requires { (bool)options.compress_integers; })
            compress_integers = options.compress_integers;
    }

    //!\brief Construct with only an output stream.
    format_output_handler(std::ostream & str) : format_output_handler(str, 1) {}
    //!\}

    //!\brief Get the header.
    var_io::header const & get_header() const
    {
        if (header == nullptr)
            throw std::runtime_error{"Attempting to read header, but no header was set."};

        return *header;
    }

    //!\brief Set the header.
    void set_header(var_io::header const & hdr)
    {
        header = {&hdr, [](var_io::header const *) {}};
    }
    //!\overload
    void set_header(var_io::header const && hdr)
    {
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
    template <typename field_types, typename field_ids>
    void write_record(record<field_types, field_ids> const & record)
    {
        write_record_impl(record);
    }

    //!\overload
    template <typename field_types, typename field_ids>
    void write_record(record<field_types, field_ids> & record)
    {
        write_record_impl(record);
    }
};

} // namespace bio
