// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::fasta.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/io.hpp>

namespace bio::io
{

/*!\brief       The FastA format.
 * \ingroup     format
 *
 * \details
 *
 * This is the FastA format tag. If you want to read FastA files, use bio::io::seq::reader, and if you want
 * to write FastA files, use bio::io::seq::writer.
 *
 * ### Introduction
 *
 * FastA is the de-facto-standard for sequence storage in bionformatics. See the
 * [article on wikipedia](https://en.wikipedia.org/wiki/FASTA_format) for a an in-depth description of the format.
 *
 * ### Fields
 *
 * The FastA format provides the fields bio::io::detail::field::seq and bio::io::detail::field::id. Both fields are
 * required when writing.
 *
 * ### Implementation notes
 *
 * When reading the ID-line, the identifier (either `;` or `>`) is stripped.
 *
 * This implementation supports the following less known and optional features of the format:
 *
 *   * ID lines beginning with `;` instead of `>`.
 *   * Line breaks and other whitespace characters in any part of the sequence.
 *   * Character counts within the sequence (they are simply ignored/stripped).
 *
 * The following optional features are currently **not supported:**
 *
 *   * Multiple comment lines (starting with either `;` or `>`); only one ID line before the sequence line is accepted.
 *
 */
struct fasta
{
    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions{
      {"fasta"},
      {"fa"},
      {"fna"},
      {"ffn"},
      {"faa"},
      {"frn"},
      {"fas"},
    };
};

} // namespace bio::io
