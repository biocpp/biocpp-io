// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::seq::writer_options.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bio/meta/tag/ttag.hpp>

#include <bio/io/format/fasta.hpp>
#include <bio/io/format/fastq.hpp>
#include <bio/io/stream/transparent_ostream.hpp>

namespace bio::io::seq
{

/*!\brief Options that can be used to configure the behaviour of bio::io::seq::writer.
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup seq
 *
 */
template <typename formats_t = meta::type_list<fasta, fastq>>
struct writer_options
{
    /*!\brief Split the sequence over multiple lines (FastA-only).
     *
     * \details
     *
     * **FastA-only**
     *
     * When set to a non-zero value, this option results in the SEQ field
     * being split over multiple lines (each of the specified length).
     */
    size_t max_seq_line_length = 70;

    /*!\brief Write two ID-lines (FastQ-only).
     *
     * \details
     *
     * **FastQ-only**
     *
     * When set, this option results in the ID being written to, both, the
     * first line in the record (starting with "@"), and the third line
     * (starting with "+").
     */
    bool double_id = false;

    /*!\brief The formats that output files can take; a bio::meta::ttag over the types.
     *
     * \details
     *
     * See bio::io::seq::writer for an overview of the the supported formats.
     */
    formats_t formats = meta::ttag<fasta, fastq>;

    //!\brief Options that are passed on to the internal stream oject.
    transparent_ostream_options stream_options{};

    /*!\brief Write legacy Windows line-endings including carriage return.
     *
     * \details
     *
     * This option results in old Windows-style line-endings ("\r\n"). Since Windows supports the typical UNIX
     * line-endigns ("\n") nowadays, this option is is highly discouraged.
     */
    bool windows_eol = false;

private:
    static_assert(meta::detail::is_type_list<formats_t>, "formats must be a bio::meta::ttag / bio::meta::type_list.");
};

} // namespace bio::io::seq
