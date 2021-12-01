// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the bio::format_output_handler_base.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/views/elements.hpp>

#include <bio/format/format_output_handler.hpp>
#include <bio/format/vcf.hpp>
#include <bio/record.hpp>
#include <bio/stream/iterator.hpp>
#include <bio/utility.hpp>
#include <bio/var_io/header.hpp>
#include <bio/var_io/misc.hpp>
#include <bio/var_io/writer_options.hpp>

namespace bio
{

template <>
class output_format_handler<vcf> : public output_format_handler_base<output_format_handler<vcf>>
{
private:
    /*!\name CRTP stuff
     * \{
     */
    using base_t = output_format_handler_base<output_format_handler<vcf>>;
    friend base_t;

    using base_t::it;
    using base_t::stream;
    using base_t::write_field;
    using base_t::write_field_aux;
    //!\}

    /*!\name State
     * \{
     */
    //!\brief Whether the header has been written or not.
    bool                                                                    header_has_been_written = false;
    //!\brief Pointer to header that can be owning or non-owning.
    std::unique_ptr<var_io::header const, void (*)(var_io::header const *)> header                  = {nullptr,
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

    //!\brief Write the elements of the range, comma-delimited. This function is not selected automatically.
    template <std::ranges::input_range rng_t, typename func_t>
    void write_delimited(rng_t && range, char const delim, func_t && func)
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
    template <std::ranges::input_range rng_t>
    void write_delimited(rng_t && range, char const delim)
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

    void write_variant(detail::is_dynamic_type auto const & var)
    {
        auto visitor = [&](auto const & param)
        {
            using param_t = std::remove_cvref_t<decltype(param)>;
            if constexpr (!std::same_as<param_t, bool>) // flags don't have any values
            {
                if constexpr (std::ranges::input_range<param_t> && !std::same_as<param_t, std::string> &&
                              !std::same_as<param_t, std::string_view>)
                {
                    write_delimited(param, ',');
                }
                else
                {
                    write_field_aux(param);
                }
            }
        };
        std::visit(visitor, var);
    }

    void write_variant(detail::is_dynamic_type auto const & var, var_io::dynamic_type_id const type_id)
    {
        if (static_cast<size_t>(type_id) != var.index())
            throw format_error{"The variant was not in the proper state."}; // TODO improve text

        write_variant(var);
    }

    template <typename header_container_t, typename field_t>
    void write_id(header_container_t const &           header_container,
                  var_io::header::idx_to_pos_t const & idx_to_pos_map,
                  field_t const &                      field)
    {
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

    template <typename key_t, typename val_t>
    void write_info_pair(std::pair<key_t, val_t> const & pair)
    {
        write_id(header->infos, header->idx_to_info_pos(), pair.first);

        if constexpr (detail::is_dynamic_type<val_t>) // all fields that aren't flags have second part
        {
            size_t pos = -1;

            if constexpr (std::integral<key_t>)
                pos = header->idx_to_info_pos().at(pair.first);
            else
                pos = header->string_to_info_pos().at(pair.first);

            var_io::dynamic_type_id type_id = header->infos[pos].type;

            if (type_id != var_io::dynamic_type_id::flag)
            {
                it = '=';
                write_variant(pair.second, type_id);
            }
        }
        else
        {
            if (!std::ranges::empty(pair.second))
            {
                it = '=';
                write_field_aux(pair.second);
            }
        }
    }

    //!\}

    /*!\name Writing individual fields - defaults (step 3)
     * \{
     */
    //!\brief This overrides default behaviour.
    template <std::ranges::input_range rng_t>
        requires std::same_as < std::ranges::range_reference_t<rng_t>,
    char > void write_field_aux(rng_t & range)
    {
        if (std::ranges::empty(range))
            it = '.';
        else
            it->write_range(range);
    }

    //!\brief This overrides default behaviour.
    template <seqan3::arithmetic num_t>
    void write_field_aux(num_t const number)
    {
        if (number == var_io::missing_value<num_t>)
            it = '.';
        else
            it->write_number(number);
    }

    //!\}

    /*!\name Field writers
     * \{
     */
    //!\brief Overload for CHROM and numeric IDs (text IDs are handled by defaults).
    void write_field(tag_t<field::chrom> /**/, auto & field)
    {
        write_id(header->contigs, header->idx_to_contig_pos(), field);
    }

    // POS, ID, REF all handled by defaults

    //!\brief Overload for ALT that is range-of-range (single-range ALT is handled by defaults).
    template <std::ranges::input_range rng_t>
        requires std::ranges::input_range<std::ranges::range_reference_t<rng_t>> // TOOD and requires
                                                                                 // write_field_aux(value)
    void write_field(tag_t<field::alt> /**/, rng_t & range) { write_delimited(range, ','); }

    // QUAL is handled by defaults

    //!\brief Overload for FILTER; single string is handled by default; single-numeric by this overload.
    void write_field(tag_t<field::filter> /**/, auto & field)
    {
        write_id(header->filters, header->idx_to_filter_pos(), field);
    }

    //!\brief Overload for FILTER; handles vector of numeric IDs and vector
    template <std::ranges::input_range rng_t>
        requires(!std::same_as<std::ranges::range_value_t<rng_t>, char>)
    void write_field(tag_t<field::filter> /**/, rng_t & range)
    {
        auto func = [this](auto const & val) { write_field(tag<field::filter>, val); };
        write_delimited(range, ';', func);
    }

    //!\brief Overload for INFO; vector of pairs.
    template <std::ranges::input_range rng_t>
        requires(tuple_like<std::ranges::range_reference_t<rng_t>>) // TODO we actually require std::pair
    void write_field(tag_t<field::info> /**/, rng_t & range)
    {
        auto func = [this](auto const & field) { write_info_pair(field); };
        write_delimited(range, ';', func);
    }

    //!\brief Overload for GENOTYPES; genotypes_as_strings
    template <std::ranges::input_range range_t>
        requires(std::ranges::input_range<std::ranges::range_reference_t<range_t>> &&
                   std::same_as<char, std::remove_cvref_t<range_innermost_value_t<range_t>>>)
    void write_field(tag_t<field::genotypes> /**/, range_t & range) { write_delimited(range, '\t'); }

    //!\brief Overload for GENOTYPES; genotypes_bcf_style
    template <std::ranges::input_range range_t>
        requires(std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<range_t>>,
                              var_io::genotype_bcf_style<ownership::shallow>> ||
                 std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<range_t>>,
                              var_io::genotype_bcf_style<ownership::deep>>)
    void write_field(tag_t<field::genotypes> /**/, range_t & range)
    {
        if (header->column_labels.size() <= 8)
            return;

        /* format field */
        auto func = [this](auto const & field) { write_id(header->formats, header->idx_to_format_pos(), field); };
        write_delimited(views::elements<0>(range), ':', func);

        if (header->column_labels.size() <= 9)
            return;

        it = '\t';

        /* sample fields */
        size_t n_samples = header->column_labels.size() - 9;

        std::vector<size_t> lengths;
        std::ranges::copy(
          views::elements<1>(range) |
            std::views::transform([](auto const & var)
                                  { return std::visit([](auto const & vec) { return vec.size(); }, var); }),
          std::back_insert_iterator{lengths});

        for (size_t i = 0; i < n_samples; ++i) // for every sample
        {
            for (size_t j = 0; j < std::ranges::size(range); ++j) // for every field
            {
                auto visitor = [&](auto const & param) // this is always a vector and sometimes vector-of-vector
                {
                    using param_t     = std::remove_cvref_t<decltype(param)>;
                    using param_ref_t = std::remove_cvref_t<std::ranges::range_reference_t<param_t>>;

                    if (i < param.size()) // param.size() is equal to lengths[j]
                    {
                        if (j > 0)
                            it = ':';

                        if constexpr (std::ranges::input_range<param_ref_t> &&
                                      !std::same_as<param_ref_t, std::string> &&
                                      !std::same_as<param_ref_t, std::string_view>)
                        {
                            write_delimited(param[i], ',');
                        }
                        else
                        {
                            write_field_aux(param[i]);
                        }
                    }
                    else
                    {
                        // when this field and all following field for this sample are empty, omit all of them
                        bool is_trailing = true;
                        for (size_t k = j; k < std::ranges::size(range); ++k)
                        {
                            if (i < lengths[k])
                            {
                                is_trailing = false;
                                break;
                            }
                        }

                        if (!is_trailing)
                        {
                            if (j > 0)
                                it = ':';
                            it = '.';
                        }
                    }
                };

                std::visit(visitor, range[j].second);
            }

            if (i < n_samples - 1)
                it = '\t';
        }
    }

    //!\brief Overload for GENOTYPES; genotypes_vcf_style
    template <typename field_t>
        requires(std::same_as<std::remove_cvref_t<field_t>, var_io::genotypes_vcf_style<ownership::shallow>> ||
                 std::same_as<std::remove_cvref_t<field_t>, var_io::genotypes_vcf_style<ownership::deep>>)
    void write_field(tag_t<field::genotypes> /**/, field_t & field)
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
        auto write_var    = [&](auto const & var) { write_variant(var); };
        auto write_sample = [&](auto const & sample) { write_delimited(sample, ':', write_var); };
        write_delimited(samples, '\t', write_sample);
    }

    //!\brief Overload for GENOTYPES; delete the generic ones.
    //     template <typename arg_t>
    //     void write_field(tag_t<field::genotypes> /**/, std::any) = delete; // prevent writing of string or number
    //!\}

    //!\brief Write the record (supports const and non-const lvalue ref).
    void write_record_impl(auto & record)
    {
        using field_ids = typename std::remove_cvref_t<decltype(record)>::field_ids;

        if (!header_has_been_written)
        {
            if (header == nullptr)
            {
                bool failed = false;

                if constexpr (field_ids::contains(field::_private))
                {
                    if (var_io::header const * ptr = get<field::_private>(record).header_ptr; ptr != nullptr)
                        set_header(*ptr);
                    else
                        failed = true;
                }
                else
                {
                    failed = true;
                }

                if (failed)
                {
                    throw std::runtime_error{
                      "You need to call set_header() on the writer/format before writing a "
                      "record."};
                }
            }

            if (write_IDX)
                it->write_range(header->to_plaintext());
            else
                it->write_range(header->to_plaintext_without_idx());

            header_has_been_written = true;
        }

        static_assert(field_ids::contains(field::chrom), "The record must contain the CHROM field.");
        write_field(tag<field::chrom>, get<field::chrom>(record));
        it = '\t';

        static_assert(field_ids::contains(field::pos), "The record must contain the POS field.");
        write_field(tag<field::pos>, get<field::pos>(record));
        it = '\t';

        if constexpr (field_ids::contains(field::id))
            write_field(tag<field::id>, get<field::id>(record));
        else
            it = '.';
        it = '\t';

        static_assert(field_ids::contains(field::ref), "The record must contain the REF field.");
        write_field(tag<field::ref>, get<field::ref>(record));
        it = '\t';

        if constexpr (field_ids::contains(field::alt))
            write_field(tag<field::alt>, get<field::alt>(record));
        else
            it = '.';
        it = '\t';

        if constexpr (field_ids::contains(field::qual))
            write_field(tag<field::qual>, get<field::qual>(record));
        else
            it = '.';
        it = '\t';

        if constexpr (field_ids::contains(field::filter))
            write_field(tag<field::filter>, get<field::filter>(record));
        else
            it = '.';
        it = '\t';

        if constexpr (field_ids::contains(field::info))
            write_field(tag<field::info>, get<field::info>(record));
        else
            it = '.';

        if (header->column_labels.size() > 8)
        {
            if constexpr (field_ids::contains(field::genotypes))
            {
                it = '\t';
                write_field(tag<field::genotypes>, get<field::genotypes>(record));
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
    output_format_handler()                              = default;
    output_format_handler(output_format_handler const &) = delete;
    output_format_handler(output_format_handler &&)      = default;
    output_format_handler & operator=(output_format_handler const &) = delete;
    output_format_handler & operator=(output_format_handler &&) = default;

    template <typename other_options_t>
    output_format_handler(std::ostream & str, other_options_t const &) : base_t{str}
    {}

    template <typename... types>
    output_format_handler(std::ostream & str, var_io::writer_options<types...> const & options) : base_t{str}
    {
        /* options */
        write_IDX   = options.write_IDX;
        windows_eol = options.windows_eol;
    }

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
