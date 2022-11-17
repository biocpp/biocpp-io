// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::seq::writer.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <filesystem>
#include <iosfwd>

#include <bio/io/detail/writer_base.hpp>
#include <bio/io/format/fasta_output_handler.hpp>
#include <bio/io/format/fastq_output_handler.hpp>
#include <bio/io/seq/record.hpp>
#include <bio/io/seq/writer_options.hpp>

namespace bio::io::seq
{

// ----------------------------------------------------------------------------
// writer
// ----------------------------------------------------------------------------

/*!\brief A class for writing sequence files, e.g. FastA and FastQ.
 * \tparam option_args_t Arguments that are forwarded to bio::io::seq::writer_options.
 * \ingroup seq
 *
 * \details
 *
 * ### Introduction
 *
 * Sequence files are the most commonly used files in Bioinformatics. They provide sequences together
 * with a (unique) identifier and (sometimes) a quality string.
 *
 * This writer supports the following formats:
 *
 *   1. FastA (see also bio::io::fasta)
 *   2. FastQ (see also bio::io::fastq)
 *
 *
 * ### Creating a writer
 *
 * Creating a writer from a filename and with default options:
 *
 * \snippet test/snippet/seq/seq_writer.cpp creation
 *
 * ### Writing a record
 *
 * Create a record using the bio::io::seq::record_dna and setting the members:
 *
 * \snippet test/snippet/seq/seq_writer.cpp simple_usage_file
 *
 * ### Writing a record without having a record
 *
 * You can use #emplace_back() to write the fields directly without creating a record first:
 *
 * \snippet test/snippet/seq/seq_writer.cpp emplace_back
 *
 * This is especially helpful if your fields exist in other separate data structures already.
 * Be aware that it is easier to mess up the order of the arguments this way.
 * The order/composition is exactly as in bio::io::seq::record!
 *
 * ### Specifying options
 *
 * Create a writer that does not split the sequence fields over multiple lines:
 *
 * \snippet test/snippet/seq/seq_writer.cpp options
 *
 * For more advanced options, see bio::io::seq::writer_options.
 *
 * ### Combining reading and writing
 *
 * This simple snippet demonstrates how to pipe from a reader into a writer (transparently converting):
 *
 * \snippet test/snippet/seq/seq_writer.cpp inout
 *
 * Views can be used to modify or filter the records before being written:
 *
 * \snippet test/snippet/seq/seq_writer.cpp inout2
 *
 * A more "traditional" programming style with the same semantic would look like this:
 *
 * \snippet test/snippet/seq/seq_writer.cpp inout3
 *
 */
template <typename... option_args_t>
class writer : public writer_base<writer<option_args_t...>, writer_options<option_args_t...>>
{
private:
    //!\brief The base class.
    using base_t      = writer_base<writer<option_args_t...>, writer_options<option_args_t...>>;
    //!\brief Inherit the format_type definition.
    using format_type = typename base_t::format_type;
    /* Implementation note
     * format_type is "inherited" as private here to avoid appearing twice in the documentation.
     * Its actual visibility is public because it is public in the base class.
     */

    //!\brief Make the format handler visible.
    using base_t::format_handler;
    //!\brief Make the init_state handler visible.
    using base_t::init_state;

    using base_t::write_record;

public:
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    writer()                           = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    writer(writer const &)             = delete;
    //!\brief Move construction is defaulted.
    writer(writer &&)                  = default;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    writer & operator=(writer const &) = delete;
    //!\brief Move assignment is defaulted.
    writer & operator=(writer &&)      = default;
    //!\brief Defaulted.
    ~writer()                          = default;

    using base_t::operator=;

    // clang-format off
    //!\copydoc bio::io::writer_base::writer_base(std::filesystem::path const & filename, format_type const & fmt, options_t const & opt = options_t{})
    // clang-format on
    writer(std::filesystem::path const &            filename,
           format_type const &                      fmt,
           writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{filename, fmt, opt}
    {}

    //!\overload
    explicit writer(std::filesystem::path const &            filename,
                    writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{filename, opt}
    {}

    // clang-format off
    //!\copydoc bio::io::writer_base::writer_base(std::ostream & str, format_type const & fmt, options_t const & opt = options_t{})
    // clang-format on
    writer(std::ostream &                           str,
           format_type const &                      fmt,
           writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{str, fmt, opt}
    {}

    //!\overload
    template <movable_ostream temporary_stream_t>
        //!\cond REQ
        requires(!std::is_lvalue_reference_v<temporary_stream_t>)
    //!\endcond
    writer(temporary_stream_t &&                    str,
           format_type const &                      fmt,
           writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{std::move(str), fmt, opt}
    {}

    /*!\brief Write a record to the file by passing individual fields.
     * \param[in] id   The ID parameter.
     * \param[in] seq  The SEQ parameter.
     * \param[in] qual The QUAL parameter.
     *
     * \details
     *
     * This function accepts exactly 3 arguments and assumes that these
     * correspond exactly to the members of the bio::io::seq::record class.
     *
     * ### Example
     *
     * \snippet test/snippet/seq/seq_writer.cpp emplace_back
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     */
    void emplace_back(auto && id, auto && seq, auto && qual)
    {
        push_back(tie_record(id, seq, qual));
    }

    /*!\brief Write a record to the file.
     * \tparam field_types Types of the fields in the record.
     * \tparam field_ids   IDs of the fields in the record.
     * \param[in] r        The record to write.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     */
    template <typename... member_types>
    void push_back(record<member_types...> const & r)
    {
        static_assert(detail::record_write_concept_checker(std::type_identity<record<member_types...>>{}));
        write_record(r);
    }

    //!\overload
    template <typename... member_types>
    void push_back(record<member_types...> & r)
    {
        static_assert(detail::record_write_concept_checker(std::type_identity<record<member_types...>>{}));
        write_record(r); // pass as non-const to allow parsing views that are not const-iterable
    }

    //!\overload
    template <typename... member_types>
    void push_back(record<member_types...> && r)
    {
        static_assert(detail::record_write_concept_checker(std::type_identity<record<member_types...>>{}));
        write_record(r); // pass as non-const to allow parsing views that are not const-iterable
    }
};

} // namespace bio::io::seq
