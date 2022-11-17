// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::format_output_handler<fasta>.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bio/io/format/fasta.hpp>
#include <bio/io/format/format_output_handler.hpp>
#include <bio/io/seq/record.hpp>
#include <bio/io/stream/detail/fast_streambuf_iterator.hpp>

namespace bio::io
{

/*!\brief Format output handler for the FastA format (bio::io::fasta).
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
 * | Member              | Type    | Default | Description                                                       |
 * |---------------------|---------|---------|-------------------------------------------------------------------|
 * |`max_seq_line_length`|`size_t` | 0       | Whether to split sequence lines after N characters.               |
 * |`windows_eol`        |`bool`   | `false` | Whether old-Windows style carriage return characters are printed. |
 *
 */
template <>
class format_output_handler<fasta> : public format_output_handler_base<format_output_handler<fasta>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The base class.
    using base_t = format_output_handler_base<format_output_handler<fasta>>;
    //!\brief Befriend the base class so we can instantiate.
    friend base_t;

    using base_t::it;
    using base_t::write_field_aux;
    //!\}

    //!\brief Insert empty line between records. [this is not an option]
    bool insert_newline = false;

    /*!\name Options
     * \{
     */
    //!\brief Break seq-lines after N characters.
    size_t max_seq_line_length = 0;

    //!\brief Write legacy Windows line-endings including carriage return.
    bool windows_eol = false;
    //!\}

    //!\brief Write the record (supports const and non-const lvalue ref).
    void write_record_impl(auto & record)
    {
        using record_t = std::remove_cvref_t<decltype(record)>;

        if (insert_newline) // TODO make an option out of this?
            it->write_end_of_line(windows_eol);

        /* ID */
        static_assert(meta::different_from<typename record_t::id_t, meta::ignore_t>,
                      "The record must contain the ID field.");
        it = '>';
        write_field_aux(record.id);
        it->write_end_of_line(windows_eol);

        /* SEQ */
        static_assert(meta::different_from<typename record_t::seq_t, meta::ignore_t>,
                      "The record must contain the SEQ field.");

        if (max_seq_line_length > 0)
        {
#if 0 // in SeqAn3, the following was slow; benchmark again!
            write_field_aux(record.id | views::to_char | views::interleave(max_seq_line_length, windows_eol ? std::string_view{"\r\n"} : std::string_view{"\n"}));
#else
            auto char_sequence = record.seq | views::to_char;
            auto cit           = std::ranges::begin(char_sequence);
            auto end           = std::ranges::end(char_sequence);

            while (cit != end)
            {
                /* Note: This solution is slightly suboptimal for sized but non-random-access ranges.*/
                auto   current_end = cit;
                size_t steps       = std::ranges::advance(current_end, max_seq_line_length, end);
                using subrange_t =
                  std::ranges::subrange<decltype(cit), decltype(cit), std::ranges::subrange_kind::sized>;
                cit = it.write_range(subrange_t{cit, current_end, (max_seq_line_length - steps)});

                it->write_end_of_line(windows_eol);
            }
#endif
        }
        else
        {
            write_field_aux(record.seq);
            it->write_end_of_line(windows_eol);
        }

        insert_newline = !std::ranges::empty(record.seq);
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
     * named members is also accepted. See bio::io::format_output_handler<fasta> for the supported options and defaults.
     */
    format_output_handler(std::ostream & str, auto const & options) : base_t{str}
    {
        // extract options
        if constexpr (requires { (size_t) options.max_seq_line_length; })
            max_seq_line_length = options.max_seq_line_length;
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
