// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::var_io::reader_options.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include <bio/alphabet/adaptation/char.hpp> // make sure that some concept checks don't fail when using strings
#include <bio/meta/tag/ttag.hpp>
#include <bio/meta/type_list/traits.hpp>

#include <bio/io/detail/misc.hpp>
#include <bio/io/format/bcf.hpp>
#include <bio/io/format/vcf.hpp>
#include <bio/io/genomic_region.hpp>
#include <bio/io/stream/transparent_istream.hpp>
#include <bio/io/stream/transparent_ostream.hpp>
#include <bio/io/var_io/header.hpp>
#include <bio/io/var_io/misc.hpp>
#include <bio/io/var_io/record.hpp>

namespace bio::io::var_io
{

//-----------------------------------------------------------------------------
// reader_options
//-----------------------------------------------------------------------------

/*!\brief Options that can be used to configure the behaviour of bio::io::var_io::reader.
 * \tparam record_t  Type of the record member (usually deduced).
 * \tparam formats_t Type of the formats member (usually deduced).
 * \ingroup var_io
 *
 * \details
 *
 * This object can be configured in a similar way as bio::io::seq_io::reader_options.
 * If you are new to the way options are set in this library, have a look at
 * bio::io::seq_io::reader first.
 *
 * ### Example
 *
 * We use value-based options and designated intialisers to simplify configuration.
 * This example shows how to read only a subset of the available fields and manually specify their type:
 *
 * \snippet test/snippet/var_io/var_io_reader_options.cpp field_types_expert
 *
 * Reading fewer fields than available may provide a noticeable speed-up since only the
 * requested fields are actually parsed.
 */
template <typename formats_t = meta::type_list<vcf, bcf>, typename record_t = record_default_shallow>
struct reader_options
{
    /*!\brief The formats that input files can take; a bio::meta::ttag over the types.
     *
     * \details
     *
     * See bio::io::var_io::reader for an overview of the the supported formats.
     */
    formats_t formats = meta::ttag<vcf, bcf>;

    //!\brief Whether to print non-critical file format warnings.
    bool print_warnings = true;

    //!\brief The record data structure; equals bio::io::var_io::record_default_shallow by default.
    record_t record{};

    /*!\name Region filtering
     * \brief These options allow filtering the file for a sub-region.
     * \{
     */
    /*!\brief Only display records that overlap with the given region (ignored if default-initialised).
     * \details
     *
     * The region must be specified as a **0-based, half-open interval** (even though VCF is 1-based).
     *
     * This option works with and without an index, but reading without an index is disabled by default.
     * Indexes are not usable when the format is uncompressed VCF or if the input is read from STDIN.
     *
     * ### Region filter with an index
     *
     * An index file allows skipping regions of the file on-disk. This is usually a lot faster than filtering
     * without an index. The filename is detected automatically but can also be specified manually (see
     * #region_index_file).
     *
     * ### Region filter without an index
     *
     * The file will be scanned linearly for records that overlap the region.
     * Only minimal parts of the record are parsed to compare the positions, so this is faster than
     * using a `std::views::filter` on the file.
     * This option is disabled by default, so an exception is thrown when a region is given but no
     * index file is found. To allow region filtering without an index, set #region_index_optional to
     * `true`.
     */
    genomic_region<ownership::deep> region{};

    /*!\brief Path to the index file [optional, auto-detected if not specified].
     * \details
     *
     * This option is ignored if no #region is specified, if the file format is uncompressed VCF or the data
     * are read from standard input.
     *
     * If no path is given, the library will look for `FILE.tbi` where `FILE` is the filename of the input
     * file. CSI-indexes are not (yet) supported.
     */
    std::filesystem::path region_index_file{};

    //!\brief Allow linear-time region filtering when no index file is available.
    bool region_index_optional = false;
    //!\}

    //!\brief Options that are passed on to the internal stream oject.
    transparent_istream_options stream_options{};

private:
    static_assert(meta::detail::is_type_list<formats_t>, "formats must be a bio::meta::ttag / bio::meta::type_list.");

    static_assert(detail::record_read_concept_checker(std::type_identity<record_t>{}));
};

} // namespace bio::io::var_io
