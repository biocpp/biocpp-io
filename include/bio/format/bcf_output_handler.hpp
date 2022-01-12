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

#include <bio/detail/magic_get.hpp>
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

    stream_buffer_exposer<char> streambuf_exposer{buffer_stream.rdbuf()};
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
                return bcf_type_descriptor::int32:
            case 16:
                return bcf_type_descriptor::int16;
            default:
                return bcf_type_descriptor::int8;
        }
    }

    detail::bcf_type_descriptor smallest_int_desc(std::signed_integral auto const num)
    {
        return smallest_int_desc(static_cast<uint64_t>(std::abs(num)));
    }

    detail::bcf_type_descriptor get_type_descriptor(auto & data, bool shrink_int = compress_integers)
    {
        using data_t = std::remove_cvref_t<decltype(data)>;

        detail::bcf_type_descriptor desc;

        if constexpr (std::arithmetic<data_t>)
        {
            desc = detail::type_2_bcf_type_descriptor<data_t>;

            if constexpr (std::integral<data_t> && sizeof(data_t) > 1)
                if (shrink_int)
                    desc = smallest_int_desc(num);
        }
        else if constexpr (detail::deliberate_alphabet<data_t>)
        {
            desc = detail::bcf_type_descriptor::char8;
        }
        else if constexpr (std::ranges::range<data_t> && detail::char_range<std::ranges::range_reference_t<data_t>)
        {
            desc = detail::bcf_type_descriptor::char8;
        }
        else if constexpr (std::ranges::range<data_t>)
        {
            using val_t = std::ranges::range_value_t<data_t>;

            if constexpr (std::arithmetic<val_t>)
            {
                desc = detail::type_2_bcf_type_descriptor<val_t>;

                if constexpr (std::integral<val_t> && sizeof(val_t) > 1)
                {
                    if (shrink_int)
                    {
                        val_t max = std::numeric_limits<val_t>::lowest();
                        for (val_t elem : data)
                        {
                            if (elem > max)
                            {
                                max = elem;
                                desc = smallest_int_desc(elem);
                                if (desc == detail::bcf_type_descriptor::int32) // this is max(type_descriptor)
                                    break;
                            }
                        }
                    }
                }
            }
            else if constexpr (detail::deliberate_alphabet<val_t>)
            {
                desc = detail::bcf_type_descriptor::char8;
            }
            else
            {
                static_assert(std::arithmetic<val_t>, "THIS SHOULDN'T HAPPEN");
            }
        }
        else
        {
            static_assert(std::arithmetic<val_t>, "THIS SHOULDN'T HAPPEN 2");
        }

        return desc;
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

    void write_single_impl(seqan3::arithmetic auto num, detail::bcf_type_descriptor const desc)
    {
        write_type_descriptor(desc, 1);
        switch (desc)
        {
            case bcf_type_descriptor::char8:
                it->write_as_binary(static_cast<char>(num));
                break;
            case bcf_type_descriptor::int8:
                it->write_as_binary(static_cast<int8_t>(num));
                break;
            case bcf_type_descriptor::int16:
                it->write_as_binary(static_cast<int16_t>(num));
                break;
            case bcf_type_descriptor::int32:
                it->write_as_binary(static_cast<int32_t>(num));
                break;
            case bcf_type_descriptor::float32:
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

        auto dump_or_convert = [&] <typename target_t> ()
        {
            if constexpr (std::same_as<target_t, elem_t> &&
                          std::ranges::contiguous_range<data_t> &&
                          std::ranges::sized_range<data_t>)
            {
                it->write_as_binary(data);
            }
            else
            {
                //TODO add something for deliberate_alphabet here
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
            case bcf_type_descriptor::char8:
                dump_or_convert<char8_t>();
                break;
            case bcf_type_descriptor::int8:
                dump_or_convert<int8_t>();
                break;
            case bcf_type_descriptor::int16:
                dump_or_convert<int16_t>();
                break;
            case bcf_type_descriptor::int32:
                dump_or_convert<int32_t>();
                break;
            case bcf_type_descriptor::float32:
                dump_or_convert<float>();
                break;
            default:
                error("Trying to write an unknown type");
        }
    }

    void write_typed_data(seqan3::arithmetic auto & num, bool shrink_int = compress_integers)
    {
        write_single_impl(num, get_type_descriptor(num, shrink_int));
    }

    void write_typed_data(detail::char_range auto & data, bool = compress_integers)
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


    template <std::ranges::forward_range data_t>
        requires seqan3::arithmetic<std::ranges::value_t<data_t>>
    void write_typed_data(data_t & data, bool shrink_int = compress_integers)
    {
        using data_t = decltype(data);
        using elem_t = std::ranges::range_value_t<data_t>;

        detail::bcf_type_descriptor desc = get_type_descriptor(data, shrink_int);

        write_type_descriptor(desc, std::ranges::distance(data));

        write_range_impl(data, desc);
    }

    template <std::ranges::forward_range data_t>
        requires detail::char_range<std::ranges::value_t<data_t>>
    void write_typed_data(data_t & data, bool = compress_integers)
    {
        //TODO
    }

    //!\brief Write chracter ranges.
    void write_field_aux(detail::char_range & range) { write_typed_range(range); }

    //!\brief Write alphabet ranges.
    template <std::ranges::forward_range rng_t>
        requires(detail::deliberate_alphabet<std::ranges::range_reference_t<rng_t>>)
    void write_field_aux(rng_t & range) { write_typed_range(range | seqan3::views::to_char); }

    //!\brief Write CStrings.
    void write_field_aux(char const * const cstr) { write_typed_range(std::string_view{cstr}); }

    //!\brief Write numbers.
    void write_field_aux(seqan3::arithmetic auto number) { write_types_single_data(number); }

#if 0
    //!\brief A range adaptor that gets the first element in a decomposable type.
    static constexpr auto views_get_first =
      std::views::transform([](auto & pair) -> decltype(auto) { return detail::get_first(pair); });

    //!\brief A range adaptor that gets the second element in a decomposable type.
    static constexpr auto views_get_second =
      std::views::transform([](auto & pair) -> decltype(auto) { return detail::get_second(pair); });

    //!\brief Get the size of a range or dynamic_vector_type.
    static size_t dyn_vec_size(auto & in)
    {
        if constexpr (detail::is_dynamic_vector_type<std::remove_cvref_t<decltype(in)>>)
            return std::visit([](auto & val) { return std::ranges::size(val); }, in);
        else
            return std::ranges::size(in);
    }
    //!\brief Write the elements of the rangeor tuple, char-delimited.
    void write_delimited(std::ranges::input_range auto && range, char const delim, auto && func)
    {
        if (std::ranges::empty(range))
            it = '.';
        else
        {
            auto b = std::ranges::begin(range);
            auto e = std::ranges::end(range);
            func(*b);
            ++b;
            for (; b != e; ++b)
            {
                it = delim;
                func(*b);
            }
        }
    }

    //!\overload
    void write_delimited(std::ranges::input_range auto && range, char const delim)
    {
        if (std::ranges::empty(range))
            it = '.';
        else
        {
            auto b = std::ranges::begin(range);
            auto e = std::ranges::end(range);
            write_field_aux(*b);
            ++b;
            for (; b != e; ++b)
            {
                it = delim;
                write_field_aux(*b);
            }
        }
    }

    //!\overload
    void write_delimited(auto && tup, char const delim, auto && func)
    {
        if constexpr (std::tuple_size_v<std::remove_cvref_t<decltype(tup)>> == 0)
            it = '.';
        else
        {
            auto pack_for_each = [&](auto &&... args)
            {
                bool first_elem = true;
                (((first_elem ? (first_elem = false, it) : it = delim), func(std::forward<decltype(args)>(args))), ...);
            };
            std::apply(pack_for_each, std::forward<decltype(tup)>(tup));
        }
    }

    //!\brief Write variant or a type that is given inplace of a variant.
    void write_variant(auto const & var)
    {
        auto visitor = [&](auto const & param)
        {
            using param_t = std::remove_cvref_t<decltype(param)>;
            if constexpr (!std::same_as<param_t, bool>) // flags don't have any values
            {
                if constexpr (std::ranges::input_range<param_t> && !detail::char_range<param_t>)
                {
                    write_delimited(param, ',');
                }
                else
                {
                    write_field_aux(param);
                }
            }
        };

        if constexpr (detail::is_dynamic_type<std::remove_cvref_t<decltype(var)>>)
            std::visit(visitor, var);
        else
            visitor(var);
    }

    //!\brief Write variant or a type that is given inplace of a variant; possibly verify.
    void write_variant(auto const & var, var_io::dynamic_type_id const type_id)
    {
        if constexpr (detail::is_dynamic_type<std::remove_cvref_t<decltype(var)>>)
        {
            if (static_cast<size_t>(type_id) != var.index())
                throw format_error{"The variant was not in the proper state."}; // TODO improve text
        }
        else
        {
            // TODO verify that type corresponds to the type_id given?
        }

        write_variant(var);
    }

    //!\brief Write a vector variant or a type that is given in place of one.
    void write_vector_variant(auto const &                var,
                              std::vector<size_t> const & lengths,
                              size_t const                i_sample,
                              size_t const                i_field)
    {
        // var is always a vector and sometimes vector-of-vector
        auto visitor = [&](auto & param)
        {
            using param_t     = std::remove_cvref_t<decltype(param)>;
            using param_ref_t = std::remove_cvref_t<std::ranges::range_reference_t<param_t>>;

            if (i_sample < param.size()) // param.size() is equal to lengths[i_field]
            {
                if (i_field > 0)
                    it = ':';

                if constexpr (std::ranges::input_range<param_ref_t> && !detail::char_range<param_ref_t>)
                    write_delimited(param[i_sample], ',');
                else
                    write_field_aux(param[i_sample]);
            }
            else
            {
                // when this field and all following field for this sample are empty, omit all of them
                bool is_trailing = true;
                for (size_t k = i_field; k < std::ranges::size(lengths); ++k)
                {
                    if (i_sample < lengths[k])
                    {
                        is_trailing = false;
                        break;
                    }
                }

                if (!is_trailing)
                {
                    if (i_field > 0)
                        it = ':';
                    it = '.';
                }
            }
        };

        if constexpr (detail::is_dynamic_vector_type<std::remove_cvref_t<decltype(var)>>)
            std::visit(visitor, var);
        else
            visitor(var);
    }

    //!\brief Implementation for writing the ID field.
    void write_id(auto const &                         header_container,
                  var_io::header::idx_to_pos_t const & idx_to_pos_map,
                  auto const &                         field)
    {
        using field_t = std::remove_cvref_t<decltype(field)>;
        if constexpr (std::integral<field_t>) // field is index
        {
            size_t pos = idx_to_pos_map.at(field);
            // TODO I don't think we need the following anymore
            if (pos >= header_container.size())
            {
                throw format_error{"The given numeric ID has no corresponding entry in the header."};
            }

            write_field_aux(header_container[pos].id);
        }
        else // probably string or string_view; write as-is
        {
            write_field_aux(field);
        }
    }

    //!\brief Implementation for writing the info field.
    void write_info_pair(auto && pair) // constraints checked in parent
    {
        using pair_t = decltype(pair);
        using key_t  = detail::first_elem_t<pair_t>;

        auto & [key, val] = pair;

        write_id(header->infos, header->idx_to_info_pos(), key);

        size_t pos = -1;

        if constexpr (std::integral<key_t>)
            pos = header->idx_to_info_pos().at(key);
        else
            pos = header->string_to_info_pos().at(key);

        var_io::dynamic_type_id type_id = header->infos[pos].type;

        if (type_id != var_io::dynamic_type_id::flag) // all fields that aren't flags have second part
        {
            it = '=';
            write_variant(val, type_id);
        }
    }

    //!\}

    /*!\name Writing individual fields - defaults (step 3)
     * \{
     */
    //!\brief This overrides default behaviour.
    template <std::ranges::input_range rng_t>
        requires(std::convertible_to<std::ranges::range_reference_t<rng_t>, char>)
    void write_field_aux(rng_t & range)
    {
        if (std::ranges::empty(range))
            it = '.';
        else
            it->write_range(range);
    }

    //!\brief This overrides default behaviour.
    void write_field_aux(seqan3::arithmetic auto const number)
    {
        using field_t = std::remove_cvref_t<decltype(number)>;
        if (number == var_io::missing_value<field_t>)
            it = '.';
        else
            it->write_number(number);
    }

    //!\}
#endif
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
    void set_core_chrom(detail::char_range auto & field)
    {
        std::string_view const buf = field; // TODO make this more generic
        if (auto it = header->string_to_contig_pos()[buf]; it == header->string_to_contig_pos().end())
            error("The contig '", field, "' is not present in the header.");
        else
            record_core.chrom = header->contigs[*it].idx;
    }

    //!\brief Overload for POS.
    void set_core_pos(std::integral auto & field)
    {
        record_core.pos = static_cast<int32_t>(field);
    }

    //!\brief Overload for rlen.
    void set_core_rlen(auto & field)
    {
        record_core.rlen = static_cast<int32_t>(std::ranges::distance(field));
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
        else if constexpr (std::ranges::forward_range<field_t>)
            record_core.n_info = std::ranges::distance(field);
        else
            static_assert(std::ranges::forward_range<field_t>, "Wrong type for INFO field."); //TODO improve diagnostic
    }

    //!\brief Overload for n_allele.
    void set_core_n_allele(auto & field)
    {
        record_core.n_allele = 1; // for REF

        using field_t = decltype(field);
        if constexpr (std::ranges::forward_range<field_t> &&
                      std::ranges::forward_range<std::ranges::range_reference_t<field_t>>)
        {
            record_core.n_allele += std::ranges::distance(field);
        }
        else // assuming single string
        {
            record_core.n_allele++;
        }
    }

    //!\brief Overload for n_fmt.
    void set_core_n_fmt(auto & field)
    {
        if constexpr (std::ranges::forward_range<decltype(field)>) // genotype_bcf_style
            record_core.n_fmt = std::ranges::distance(field);
        else // genotypes_vcf_style
            record_core.n_fmt = std::ranges::distance(detail::get_first(field));
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
                auto it = header->string_to_filter_pos()[text_id];

                if (it == header->string_to_filter_pos().end())
                    error("The filter '", text_id, "' is not present in the header.");

                return header->filters[*it].idx;
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
        using id_t = detail::first_elem_t<decltype(info_element)>;
        using value_t = detail::second_elem_t<decltype(info_element)>;

        /* ID */
        int32_t idx = 0;
        if constexpr (std::integral<id_t>)
        {
            assert((int64_t)detail::get_first(info_element) <= (int64_t)header->max_idx());
            idx = detail::get_first(info_element);
        }
        else
        {
            auto it = header->string_to_info_pos()[detail::get_first(info_element)];
            if (it == header->string_to_info_pos().end())
                error("The info '", field, "' is not present in the header.");

            idx = header->infos[*it].idx;
        }
        write_single_impl(idx, header_dictionary_desc);

        /* VALUE */
        if constexpr (detail::is_dynamic_type<value_t>)
        {
            auto func = [&] (auto & param) { write_typed_data(param); };
            std::visit(func, detail::get_second(info_element));
        }
        else
        {
            write_typed_data(detail::get_second(info_element));
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
        auto func = [&](auto const & field) { write_info_element(field); };
        std::apply(func, tup);
    }

    void write_genotypes_element(auto & genotype)
    {
        auto & [ id, value ] = genotype;

        /* ID */
        int32_t idx = 0;
        if constexpr (std::integral<id_t>)
        {
            assert((int64_t)id <= (int64_t)header->max_idx());
            idx = id;
        }
        else
        {
            auto it = header->string_to_genotype_pos()[id];
            if (it == header->string_to_genotype_pos().end())
                error("The genotype '", id, "' is not present in the header.");

            idx = header->infos[*it].idx;
        }
        write_single_impl(idx, header_dictionary_desc);

        /* value */
        auto func = [](auto & value)
        {
            if (std::ranges::size(value) != 0 && std::ranges::size(value) != record_core.n_sample)
            {
                error("There are ", std::ranges::size(value), " values in the genotype vector, "
                      "but there must be exactly one per sample or none at all.");
            }

            size_t max_length = 1;
            if constexpr (std::ranges::range<decltype(value)>)
            {
                auto to_size = [] (auto & r) { return std::ranges::size(r); };
                max_length = std::ranges::max(value, {}, to_size);

                //TODO what if value is empty?
                detail::bcf_type_descriptor desc = get_type_descriptor(*std::ranges::begin(value),
                                                                       false /* no shrinking here */);

                //TODO if shrink, manually look for max over concatenated input




            }


        };




    }

    //!\brief Overload for GENOTYPES; genotypes_bcf_style.
    template <std::ranges::forward_range range_t>
        requires(detail::genotype_bcf_style_writer_concept<std::ranges::range_reference_t<range_t>>)
    void write_field(vtag_t<field::genotypes> /**/, range_t & range)
    {
        for (auto && genotype : range)
            write_genotype_bcf_style_impl(genotype);
    }


        if (header->column_labels.size() <= 8)
            return;

        /* format field */
        auto func = [this](auto const & field) { write_id(header->formats, header->idx_to_format_pos(), field); };
        write_delimited(range | views_get_first, ':', func);

        if (header->column_labels.size() <= 9)
            return;

        it = '\t';

        /* sample fields */
        size_t n_samples = header->column_labels.size() - 9;

        std::vector<size_t> lengths;
        std::ranges::copy(range | views_get_second |
                            std::views::transform([](auto & var) { return dyn_vec_size(var); }),
                          std::back_insert_iterator{lengths});

        for (size_t i_sample = 0; i_sample < n_samples; ++i_sample) // for every sample
        {
            for (size_t i_field = 0; i_field < std::ranges::size(range); ++i_field) // for every field
                write_vector_variant(detail::get_second(range[i_field]), lengths, i_sample, i_field);

            if (i_sample < n_samples - 1)
                it = '\t';
        }
    }

    //!\brief Overload for GENOTYPES; genotypes_bcf_style; nonvariant
    template <typename... elem_ts>
        requires(detail::genotype_bcf_style_writer_concept<std::remove_cvref_t<elem_ts>> &&...)
    void write_field(vtag_t<field::genotypes> /**/, std::tuple<elem_ts...> & tup)
    {
        if (header->column_labels.size() <= 8)
            return;

        /* format field */
        auto print_format = [&](auto & pair)
        { write_id(header->formats, header->idx_to_format_pos(), detail::get_first(pair)); };
        write_delimited(tup, ':', print_format);

        if (header->column_labels.size() <= 9)
            return;

        it = '\t';

        /* sample fields */
        size_t n_samples = header->column_labels.size() - 9;

        std::vector<size_t> lengths;
        auto                add_sizes = [&](auto const &... pairs)
        { (lengths.push_back(dyn_vec_size(detail::get_second(pairs))), ...); };
        std::apply(add_sizes, tup);

        for (size_t i_sample = 0; i_sample < n_samples; ++i_sample) // for every sample
        {
            auto write_fields = [&]<typename tup_t, size_t... i_field>(tup_t const & tup,
                                                                       std::index_sequence<i_field...>)
            {
                (write_vector_variant(detail::get_second(std::get<i_field>(tup)), lengths, i_sample, i_field), ...);
            };

            write_fields(tup, std::make_index_sequence<sizeof...(elem_ts)>{});

            if (i_sample < n_samples - 1)
                it = '\t';
        }
    }

    //!\brief Overload for GENOTYPES; genotypes_vcf_style
    void write_field(vtag_t<field::genotypes> /**/, detail::genotypes_vcf_style_writer_concept auto & field)
    {
        if (header->column_labels.size() <= 8)
            return;

        auto & [format, samples] = field;

        /* format field */
        write_delimited(format, ':');
        it = '\t';

        if (header->column_labels.size() <= 9)
            return;

        /* samples */
        // functional programming for the win! [this works for tuples and ranges!]
        auto write_var    = [&](auto const & var) { write_variant(var); };
        auto write_sample = [&](auto const & sample) { write_delimited(sample, ':', write_var); };
        write_delimited(samples, '\t', write_sample);
    }
    //!\}
#endif
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
        this_record_offset = stream_buffer_exposer.pptr() - stream_buffer_exposer.pbase();

        it->write_number_as_binary(uint32_t{}); // l_shared
        uint32_t l_shared_tmp = 0;
        it->write_number_as_binary(uint32_t{}); // l_indiv
        uint32_t l_indiv_tmp = 0;

        /* this prepares the record_core and doesn't write anything til the end */
        record_core = {}; // reset this
        set_core_chrom(get<field::chrom>(record));
        set_core_pos(get<field::pos>(record));

        set_core_rlen(get<field::ref>(record));

        if constexpr (field_ids::contains(field::qual))
            set_core_qual(vtag<field::qual>, get<field::qual>(record));

        if constexpr (field_ids::contains(field::info))
            set_core_n_info(get<field::info>(record));

        if constexpr (field_ids::contains(field::alt))
            set_core_n_allele(get<field::alt>(record));

        record_core.n_sample = header->column_labels.size() > 9 ? header->column_labels.size() - 9 : 0;

        if constexpr (field_ids::contains(field::genotypes))
            set_core_n_fmt(vtag<field::genotypes>, get<field::genotypes>(record));

        // write record core
        it->write_as_binary(record_core);

        /* After this point, the order of writers is important! */

        if constexpr (field_ids::contains(field::id))
            write_field(vtag<field::id>, get<field::id>(record));
        else
            write_field(vtag<field::id>, detail::missing_value<std::string_view>);

        write_field(vtag<field::ref>, get<field::ref>(record));

        if constexpr (field_ids::contains(field::alt))
            write_field(vtag<field::alt>, get<field::alt>(record));
        else
            write_field(vtag<field::alt>, detail::missing_value<std::string_view>);

        if constexpr (field_ids::contains(field::filter))
            write_field(vtag<field::filter>, get<field::filter>(record));
        else
            write_field(vtag<field::filter>, detail::missing_value<std::string_view>);

        if constexpr (field_ids::contains(field::info))
            write_field(vtag<field::info>, get<field::info>(record));
        else
            write_field(vtag<field::info>, detail::missing_value<std::string_view>);

        l_shared_tmp = stream_buffer_exposer.pptr() - stream_buffer_exposer.ebase() - this_record_offset;

        if (header->column_labels.size() > 8)
        {
            if constexpr (field_ids::contains(field::genotypes))
                write_field(vtag<field::genotypes>, get<field::genotypes>(record));
            else
                write_field(vtag<field::genotypes>, std::tuple<>{});
        }

        l_indiv_tmp = stream_buffer_exposer.pptr() - stream_buffer_exposer.pbase() - this_record_offset - l_shared_tmp;

        // now we get a pointer to where l_shared and l_indiv are stored in the buffer and write them
        uint32_t * iptr = reinterpret_cast<uint32_t*>(stream_buffer_exposer.pbase() + this_record_offset);
        *iptr = l_shared_tmp;
        ++iptr;
        *iptr = l_indiv_tmp;

        // possibly turn this into an option or take from the transparent_istream options
        size_t min_write_size = 10 * 1024 * 1024; // 10MB
        if (stream_buffer_exposer.pptr() - stream_buffer_exposer.pbase() > min_write_size)
        {
            // write data from the buffer_stream into the actual stream
            // this *should* be unbuffered for large min_write_sizes
            stream.rdbuf().sputn(stream_buffer_exposer.pbase(), stream_buffer_exposer.pptr() - stream_buffer_exposer.pbase());

            // reset the buffer_stream(first argument sets pbase AND pptr):
            stream_buffer_exposer.setp(stream_buffer_exposer.pbase(), stream_buffer_exposer.epptr());
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
        if constexpr (requires { (bool)options.write_IDX; })
            write_IDX = options.write_IDX;

        if constexpr (requires { (bool)options.windows_eol; })
            windows_eol = options.windows_eol;
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
