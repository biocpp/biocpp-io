// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::format_output_handler<vcf>.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bio/detail/magic_get.hpp>
#include <bio/format/format_output_handler.hpp>
#include <bio/format/vcf.hpp>
#include <bio/record.hpp>
#include <bio/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/var_io/header.hpp>
#include <bio/var_io/misc.hpp>
#include <bio/var_io/writer_options.hpp>

namespace bio
{

/*!\brief Format output handler for the VCF format (bio::vcf).
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
 * | Member          | Type    | Default | Description                                                       |
 * |-----------------|---------|---------|-------------------------------------------------------------------|
 * |`write_IDX`      |`bool`   | `false` | Whether IDX values (used internally, see BCF spec) are printed.   |
 * |`windows_eol`    |`bool`   | `false` | Whether old-Windows style carriage return characters are printed. |
 *
 * ### Performance
 *
 * TODO after genotype redesign
 */
template <>
class format_output_handler<vcf> : public format_output_handler_base<format_output_handler<vcf>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The base class.
    using base_t = format_output_handler_base<format_output_handler<vcf>>;
    //!\brief Befriend the base class so we can instantiate.
    friend base_t;

    using base_t::it;
    using base_t::stream;
    using base_t::write_field;
    using base_t::write_field_aux;
    //!\}

    /*!\name State
     * \{
     */
    //!\brief Can be used to infer if the object is in moved-from state.
    detail::move_tracker move_tracker;
    //!\brief Whether the header has been written or not.
    bool                 header_has_been_written = false;

    //!\brief Pointer to header that can be owning or non-owning.
    std::unique_ptr<var_io::header const, void (*)(var_io::header const *)> header = {nullptr,
                                                                                      [](var_io::header const *) {}};
    //!\}

    /*!\name Options
     * \{
     */
    //!\brief Write legacy Windows line-endings including carriage return.
    bool windows_eol = false;
    //!\brief Write BCF's IDX fields from into the header->
    bool write_IDX   = false;
    //!\}

    /*!\name Arbitrary helpers
     * \{
     */
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
            if (!detail::type_id_is_compatible(type_id, var_io::dynamic_type_id{static_cast<int32_t>(var.index())}))
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

    /*!\name Field writers
     * \{
     */
    //!\brief Overload for CHROM and numeric IDs (text IDs are handled by defaults).
    void write_field(vtag_t<field::chrom> /**/, auto & field)
    {
        write_id(header->contigs, header->idx_to_contig_pos(), field);
    }

    // POS, ID, REF all handled by defaults

    //!\brief Overload for ALT that is range-of-range (single-range ALT is handled by defaults).
    template <std::ranges::input_range rng_t>
        requires std::ranges::input_range<std::ranges::range_reference_t<rng_t>> // TOOD and requires
                                                                                 // write_field_aux(value)
    void write_field(vtag_t<field::alt> /**/, rng_t & range) { write_delimited(range, ','); }

    // QUAL is handled by defaults

    //!\brief Overload for FILTER; single string is handled by default; single-numeric by this overload.
    void write_field(vtag_t<field::filter> /**/, auto & field)
    {
        write_id(header->filters, header->idx_to_filter_pos(), field);
    }

    //!\brief Overload for FILTER; handles vector of numeric IDs and vector
    template <std::ranges::input_range rng_t>
        requires(!std::same_as<std::ranges::range_value_t<rng_t>, char>)
    void write_field(vtag_t<field::filter> /**/, rng_t & range)
    {
        auto func = [this](auto const & val) { write_field(vtag<field::filter>, val); };
        write_delimited(range, ';', func);
    }

    //!\brief Overload for INFO; range of pairs.
    template <std::ranges::input_range rng_t>
        requires(detail::info_element_writer_concept<std::ranges::range_reference_t<rng_t>>)
    void write_field(vtag_t<field::info> /**/, rng_t & range)
    {
        auto func = [&](auto const & field) { write_info_pair(field); };
        write_delimited(range, ';', func);
    }

    //!\brief Overload for INFO; range of pairs.
    template <typename... elem_ts>
        requires(detail::info_element_writer_concept<elem_ts> &&...)
    void write_field(vtag_t<field::info> /**/, std::tuple<elem_ts...> & tup) // TODO add const version
    {
        auto func = [&](auto const & field) { write_info_pair(field); };
        write_delimited(tup, ';', func);
    }

    //!\brief Overload for GENOTYPES; genotypes_bcf_style.
    template <std::ranges::forward_range range_t>
        requires(detail::genotype_bcf_style_writer_concept<std::ranges::range_reference_t<range_t>>)
    void write_field(vtag_t<field::genotypes> /**/, range_t & range)
    {
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
            auto write_fields =
              [&]<typename tup_t, size_t... i_field>(tup_t const & tup, std::index_sequence<i_field...>)
            { (write_vector_variant(detail::get_second(std::get<i_field>(tup)), lengths, i_sample, i_field), ...); };

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

    //!\brief Write the header.
    void write_header() noexcept(false)
    {
        if (header == nullptr)
            throw missing_header_error{"The VCF output handler needs a header but no header was set."};

        if (write_IDX)
            it->write_range(header->to_plaintext());
        else
            it->write_range(header->to_plaintext_without_idx());

        header_has_been_written = true;
    }

    //!\brief Write the record (supports const and non-const lvalue ref).
    void write_record_impl(auto & record)
    {
        using field_ids = typename std::remove_cvref_t<decltype(record)>::field_ids;

        if (!header_has_been_written)
        {
            if (header == nullptr)
            {
                if constexpr (field_ids::contains(field::_private))
                {
                    if (var_io::header const * ptr = get<field::_private>(record).header_ptr; ptr != nullptr)
                        set_header(*ptr);
                }
            }

            write_header();
        }

        static_assert(field_ids::contains(field::chrom), "The record must contain the CHROM field.");
        write_field(vtag<field::chrom>, get<field::chrom>(record));
        it = '\t';

        static_assert(field_ids::contains(field::pos), "The record must contain the POS field.");
        write_field(vtag<field::pos>, get<field::pos>(record));
        it = '\t';

        if constexpr (field_ids::contains(field::id))
            write_field(vtag<field::id>, get<field::id>(record));
        else
            it = '.';
        it = '\t';

        static_assert(field_ids::contains(field::ref), "The record must contain the REF field.");
        write_field(vtag<field::ref>, get<field::ref>(record));
        it = '\t';

        if constexpr (field_ids::contains(field::alt))
            write_field(vtag<field::alt>, get<field::alt>(record));
        else
            it = '.';
        it = '\t';

        if constexpr (field_ids::contains(field::qual))
            write_field(vtag<field::qual>, get<field::qual>(record));
        else
            it = '.';
        it = '\t';

        if constexpr (field_ids::contains(field::filter))
            write_field(vtag<field::filter>, get<field::filter>(record));
        else
            it = '.';
        it = '\t';

        if constexpr (field_ids::contains(field::info))
            write_field(vtag<field::info>, get<field::info>(record));
        else
            it = '.';

        if (header->column_labels.size() > 8)
        {
            if constexpr (field_ids::contains(field::genotypes))
            {
                it = '\t';
                write_field(vtag<field::genotypes>, get<field::genotypes>(record));
            }
            else
            {
                for (size_t i = 8; i < header->column_labels.size(); ++i)
                {
                    it = '\t';
                    it = '.';
                }
            }
        }

        it->write_end_of_line(windows_eol);
    }

public:
    /*!\name Constructors, destructor and assignment.
     * \brief These are all private to prevent wrong instantiation.
     * \{
     */
    format_output_handler()                                          = delete;  //!< Defaulted.
    format_output_handler(format_output_handler const &)             = delete;  //!< Deleted.
    format_output_handler(format_output_handler &&)                  = default; //!< Defaulted.
    format_output_handler & operator=(format_output_handler const &) = delete;  //!< Deleted.
    format_output_handler & operator=(format_output_handler &&)      = default; //!< Defaulted.

    /*!\brief Construct with an options object.
     * \param[in,out] str The output stream.
     * \param[in] options An object with options for the output handler.
     * \details
     *
     * The options argument is typically bio::var_io::writer_options, but any object with a subset of similarly named
     * members is also accepted. See bio::format_output_handler<vcf> for the supported options and defaults.
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

    //!\brief The destructor writes the header if necessary and cleans up.
    ~format_output_handler() noexcept(false)
    {
        // never throw if the stack is unwinding
        if (std::uncaught_exceptions() > 0)
            return;

        // no cleanup is needed if we are in moved-from state
        if (move_tracker.moved_from)
            return;

        // if no records were written, the header also wasn't written, but needs to be:
        if (!header_has_been_written)
            write_header();
    }
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
