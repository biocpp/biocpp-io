// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::format_output_handler<bed>.
 * \author Joshua Kim <kim_j AT molgen.mpg.de>
 */

#pragma once

#include <bio/detail/magic_get.hpp>
#include <bio/format/format_output_handler.hpp>
#include <bio/format/bed.hpp>
#include <bio/record.hpp>
#include <bio/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/ann_io/header.hpp>
#include <bio/ann_io/misc.hpp>
#include <bio/ann_io/writer_options.hpp>

namespace bio
{

/*!\brief Format output handler for the BED format (bio::bed).
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
 * |`windows_eol`    |`bool`   | `false` | Whether old-Windows style carriage return characters are printed. |
 *
 * ### Performance
 *
 * TODO after genotype redesign
 */
template <>
class format_output_handler<bed> : public format_output_handler_base<format_output_handler<bed>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The base class.
    using base_t = format_output_handler_base<format_output_handler<bed>>;
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
    std::unique_ptr<ann_io::header const, void (*)(ann_io::header const *)> header = {nullptr,
                                                                                      [](ann_io::header const *) {}};
    //!\}

    /*!\name Options
     * \{
     */
    //!\brief Write legacy Windows line-endings including carriage return.
    bool windows_eol = false;
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

    //!\brief Write the elements of the range or tuple, char-delimited.
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
    //!\}

    //!\brief Write the header.
    void write_header()
    {
        if (header != nullptr)
            it->write_range(header->to_plaintext());

        header_has_been_written = true;
    }

    //!\brief Write the record (supports const and non-const lvalue ref).
    void write_record_impl(auto & record)
    {
        using field_ids = typename std::remove_cvref_t<decltype(record)>::field_ids;

        if (!header_has_been_written)
            write_header();

        static_assert(field_ids::contains(field::chrom), "The record must contain the chrom field.");
        write_field(vtag<field::chrom>, get<field::chrom>(record));
        it = '\t';

        static_assert(field_ids::contains(field::chromStart), "The record must contain the chromStart field.");
        write_field(vtag<field::chromStart>, get<field::chromStart>(record));
        it = '\t';

        static_assert(field_ids::contains(field::chromEnd), "The record must contain the chromEnd field.");
        write_field(vtag<field::chromEnd>, get<field::chromEnd>(record));

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
     * The options argument is typically bio::ann_io::writer_options, but any object with a subset of similarly named
     * members is also accepted. See bio::format_output_handler<vcf> for the supported options and defaults.
     */
    format_output_handler(std::ostream & str, auto const & options) : base_t{str}
    {
        // extract options
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
    ann_io::header const & get_header() const
    {
        if (header == nullptr)
            throw missing_header_error{"Attempting to read header, but no header was set."};

        return *header;
    }

    //!\brief Set the header.
    void set_header(ann_io::header const & hdr)
    {
        header = {&hdr, [](ann_io::header const *) {}};
    }
    //!\overload
    void set_header(ann_io::header const && hdr)
    {
        header = {new ann_io::header(std::move(hdr)), [](ann_io::header const * ptr) { delete ptr; }};
    }
    //!\overload
    void set_header(ann_io::header & hdr)
    {
        set_header(std::as_const(hdr));
    }
    //!\overload
    void set_header(ann_io::header && hdr)
    {
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
