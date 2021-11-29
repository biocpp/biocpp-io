// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::seq_io::reader and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/misc.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/utility/views/to.hpp>

#include <bio/detail/reader_base.hpp>
#include <bio/format/fasta_input_handler.hpp>
#include <bio/seq_io/reader_options.hpp>

namespace bio::seq_io
{

// ----------------------------------------------------------------------------
// reader
// ----------------------------------------------------------------------------

/*!\brief A class for reading sequence files, e.g. FASTA, FASTQ.
 * \ingroup seq_io
 *
 * \details
 *
 * ### Introduction
 *
 * Sequence files are the most generic and common biological files. Well-known formats include
 * FastA and FastQ, but sometimes you may also be interested in treating SAM or BAM files as sequence
 * files, discarding the alignment.
 *
 * The Sequence I/O reader supports reading the following fields:
 *
 *   1. bio::field::seq
 *   2. bio::field::id
 *   3. bio::field::qual
 *
 * And it supports the following formats:
 *
 *   1. FASTA (see also bio::fasta)
 *
 * ### Simple usage
 *
 * Iterate over a sequence file via the reader and print the record's contents:
 *
 * \snippet test/snippet/seq_io/snippet_reader.cpp simple_usage_file
 *
 * Read from standard input instead of a file:
 *
 * \snippet test/snippet/seq_io/snippet_reader.cpp simple_usage_stream
 *
 * ### Decomposed iteration
 *
 * The records can be decomposed on-the-fly using
 * [structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding):
 *
 * \snippet test/snippet/seq_io/snippet_reader.cpp decomposed
 *
 * Note that the order of the fields is defined by bio::seq_io::default_field_ids and independent
 * of the names you give to the binding.
 *
 * ### Views on files
 *
 * This iterates over the first five records where the sequence length is at least 10:
 *
 * \snippet test/snippet/seq_io/snippet_reader.cpp views
 *
 * ### Specifying options
 *
 * This snippet demonstrates how to read sequence data as protein data and have the IDs truncated
 * at the first whitespace:
 * \snippet test/snippet/seq_io/snippet_reader.cpp options
 *
 * For more advanced options, see bio::seq_io::reader_options.
 */
template <typename... option_args_t>
class reader : public reader_base<reader_options<option_args_t...>>
{
private:
    //!\brief The base class.
    using base_t      = reader_base<reader_options<option_args_t...>>;
    //!\brief Inherit the format_type definition.
    using format_type = typename base_t::format_type;
    /* Implementation note
     * format_type is "inherited" as private here to avoid appearing twice in the documentation.
     * Its actual visibility is public because it is public in the base class.
     */
public:
    // clang-format off
    //!\copydoc bio::reader_base::reader_base(std::filesystem::path const & filename, format_type const & fmt, options_t const & opt = options_t{})
    // clang-format on
    reader(std::filesystem::path const &            filename,
           format_type const &                      fmt,
           reader_options<option_args_t...> const & opt = reader_options<option_args_t...>{}) :
      base_t{filename, fmt, opt}
    {}

    //!\overload
    explicit reader(std::filesystem::path const &            filename,
                    reader_options<option_args_t...> const & opt = reader_options<option_args_t...>{}) :
      base_t{filename, opt}
    {}

    // clang-format off
    //!\copydoc bio::reader_base::reader_base(std::istream & str, format_type const & fmt, options_t const & opt = options_t{})
    // clang-format on
    reader(std::istream &                           str,
           format_type const &                      fmt,
           reader_options<option_args_t...> const & opt = reader_options<option_args_t...>{}) :
      base_t{str, fmt, opt}
    {}

    //!\overload
    template <movable_istream temporary_stream_t>
        //!\cond REQ
        requires(!std::is_lvalue_reference_v<temporary_stream_t>)
    //!\endcond
    reader(temporary_stream_t &&                    str,
           format_type const &                      fmt,
           reader_options<option_args_t...> const & opt = reader_options<option_args_t...>{}) :
      base_t{std::move(str), fmt, opt}
    {}
};

} // namespace bio::seq_io
