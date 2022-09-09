// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::seq_io::reader_options and various field_types definitions.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/alphabet/adaptation/char.hpp>
#include <bio/alphabet/aminoacid/aa27.hpp>
#include <bio/alphabet/concept.hpp>
#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/quality/phred63.hpp>
#include <bio/meta/tag/ttag.hpp>
#include <bio/meta/type_list/traits.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/detail/concept.hpp>
#include <bio/io/detail/misc.hpp>
#include <bio/io/detail/range.hpp>
#include <bio/io/format/fasta.hpp>
#include <bio/io/format/fastq.hpp>
#include <bio/io/misc.hpp>
#include <bio/io/seq_io/record.hpp>
#include <bio/io/stream/transparent_istream.hpp>

namespace bio::io::seq_io
{

/*!\brief Options that can be used to configure the behaviour of bio::io::seq_io::reader.
 * \ingroup seq_io
 * \tparam record_t      Type of the record member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup seq_io
 * \details
 *
 * The options can be used to configure the behaviour of the reader. An important option is
 * the record type. By default, the reader assumes that sequence data is DNA, but this can be changed
 * to any other bio::alphabet::alphabet or the plain `char` type which accepts any values.
 *
 * Furthermore, different container and view types can be chosen. See bio::io::seq_io::record
 * for more details.
 *
 * ### Example
 *
 * Options can be easily set via [designated
 * initialisers](https://en.cppreference.com/w/cpp/language/aggregate_initialization).
 *
 * To switch from DNA reading (the default) to protein reading, activate truncating of IDs, and set the maximum
 * number of decoder threads (affect only bgzipped files), do the following:
 *
 * \snippet test/snippet/seq_io/seq_io_reader_options.cpp example_simple
 *
 * It is not required to specify all options; default values are documented.
 *
 * Please be aware that those options that you modify need to be set in the correct order -- **which is
 * alphabetical order** for all option classes in this library.
 *
 * Typically, the options are set as part of the bio::io::seq_io::reader construction. See the respective
 * documentation page for more information.
 *
 * ### Example (advanced)
 *
 * This code makes FASTA the only legal format and creates records with only the ID field as
 * a std::string:
 *
 * \snippet test/snippet/seq_io/seq_io_reader_options.cpp example_advanced
 */
template <typename record_t = record_dna_shallow, typename formats_t = meta::type_list<fasta, fastq>>
struct reader_options
{
    /*!\brief The formats that input files can take; a bio::meta::ttag over the types.
     *
     * \details
     *
     * See bio::io::seq_io::reader for an overview of the the supported formats.
     */
    formats_t formats = meta::ttag<fasta, fastq>;

    /*!\brief A record that can store the fields read from disk.
     * \details
     *
     * See bio::io::seq_io::record for details.
     */
    record_t record{};

    //!\brief Options that are passed on to the internal stream oject.
    transparent_istream_options stream_options{};

    //!\brief Truncate IDs at first whitespace.
    bool truncate_ids = false;

private:
    static_assert(meta::detail::is_type_list<formats_t>, "formats must be a bio::meta::ttag / bio::meta::type_list.");

    static_assert(detail::record_read_concept_checker(std::type_identity<record_t>{}));
};

} // namespace bio::io::seq_io
