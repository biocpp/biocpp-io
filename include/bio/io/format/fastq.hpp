// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::fastq.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/io/platform.hpp>

namespace bio
{

/*!\brief       The FastQ format.
 * \ingroup     format
 *
 * \details
 *
 * This is the FastQ format tag. If you want to read FastQ files, use bio::seq_io::reader, and if you want
 * to write FastQ files, use bio::seq_io::writer.
 *
 * ### Introduction
 *
 * FastQ is the de-facto-standard for read data in bionformatics. See the
 * [article on wikipedia](https://en.wikipedia.org/wiki/FASTQ_format) for a an in-depth description of the format.
 *
 * ### Fields
 *
 * The FastQ format provides the fields bio::field::seq, bio::field::id and bio::field::qual.
 * All fields are required when writing.
 *
 * ### Implementation notes
 *
 * The following optional features are supported by the implementation:
 *
 *   * Second ID line as third line (after `+`).
 *
 * The following optional features are *not* supported by the implementation:
 *
 *   * Linebreaks within any field.
 */
struct fastq
{
    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions{"fastq", "fq"};
};

} // namespace bio
