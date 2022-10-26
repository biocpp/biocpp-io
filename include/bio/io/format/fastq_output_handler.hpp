// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::format_output_handler<fastq>.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bio/io/format/fastq.hpp>
#include <bio/io/format/format_output_handler.hpp>
#include <bio/io/seq/record.hpp>
#include <bio/io/stream/detail/fast_streambuf_iterator.hpp>

namespace bio::io
{

/*!\brief Format output handler for the FastQ format (bio::io::fastq).
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
 * |`double_id`      |`bool`   | `false` | Whether the ID is also written to the "+"-line.                   |
 * |`windows_eol`    |`bool`   | `false` | Whether old-Windows style carriage return characters are printed. |
 *
 */
template <>
class format_output_handler<fastq> : public format_output_handler_base<format_output_handler<fastq>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The base class.
    using base_t = format_output_handler_base<format_output_handler<fastq>>;
    //!\brief Befriend the base class so we can instantiate.
    friend base_t;

    using base_t::it;
    using base_t::stream;
    using base_t::write_field_aux;
    //!\}

    /*!\name Options
     * \{
     */
    //!\brief Write ID also to "+" line.
    bool double_id = false;

    //!\brief Write legacy Windows line-endings including carriage return.
    bool windows_eol = false;
    //!\}

    //!\brief Write the record (supports const and non-const lvalue ref).
    void write_record_impl(auto & record)
    {
        using record_t = std::remove_cvref_t<decltype(record)>;

        /* ID */
        static_assert(meta::different_from<typename record_t::id_t, meta::ignore_t>,
                      "The record must contain the ID field.");
        it = '@';
        write_field_aux(record.id);
        it->write_end_of_line(windows_eol);

        /* SEQ */
        static_assert(meta::different_from<typename record_t::seq_t, meta::ignore_t>,
                      "The record must contain the SEQ field.");
        write_field_aux(record.seq);
        it->write_end_of_line(windows_eol);

        /* + */
        it = '+';
        if (double_id)
            write_field_aux(record.id);
        it->write_end_of_line(windows_eol);

        /* QUAL */
        static_assert(meta::different_from<typename record_t::qual_t, meta::ignore_t>,
                      "The record must contain the QUAL field.");

        if constexpr (std::ranges::sized_range<typename record_t::seq_t> &&
                      std::ranges::sized_range<typename record_t::qual_t>)
        {
            if (std::ranges::size(record.seq) != std::ranges::size(record.qual))
                throw format_error{"The SEQ and QUAL fields must have the same length."};
        }

        write_field_aux(record.qual);
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
    ~format_output_handler()                                         = default; //!< Defaulted.

    /*!\brief Construct with an options object.
     * \param[in,out] str The output stream.
     * \param[in] options An object with options for the output handler.
     * \details
     *
     * The options argument is typically bio::io::seq::writer_options, but any object with a subset of similarly
     * named members is also accepted. See bio::io::format_output_handler<fastq> for the supported options and defaults.
     */
    format_output_handler(std::ostream & str, auto const & options) : base_t{str}
    {
        // extract options
        if constexpr (requires { (bool)options.double_id; })
            double_id = options.double_id;
        if constexpr (requires { (bool)options.windows_eol; })
            windows_eol = options.windows_eol;
    }

    //!\brief Construct with only an output stream.
    format_output_handler(std::ostream & str) : format_output_handler(str, 1) {}
    //!\}

    //!\brief Write the record.
    template <typename... field_types>
    void write_record(seq::record<field_types...> const & record)
    {
        write_record_impl(record);
    }

    //!\overload
    template <typename... field_types>
    void write_record(seq::record<field_types...> & record)
    {
        write_record_impl(record);
    }
};

} // namespace bio::io
