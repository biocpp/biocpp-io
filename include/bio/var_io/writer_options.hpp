// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::var_io::writer_options.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bio/format/vcf.hpp>
#include <bio/stream/transparent_ostream.hpp>
#include <bio/var_io/misc.hpp>

namespace bio::var_io
{

/*!\brief Options that can be used to configure the behaviour of bio::var_io::writer.
 * \tparam field_ids_t   Type of the field_ids member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup var_io
 *
 * \details
 *
 * TODO describe how to easily initialise this
 */
template <typename field_ids_t = decltype(default_field_ids), typename formats_t = type_list<vcf>>
struct writer_options
{
    //!\brief The fields that shall be contained in each record; a seqan3::tag over seqan3::field.
    field_ids_t field_ids{};

    /*!\brief The formats that output files can take; a seqan3::type_tag over the types.
     *
     * \details
     *
     * See bio::var_io::writer for an overview of the the supported formats.
     */
    formats_t formats{};

    //!\brief Whether to print non-critical file format warnings.
    bool print_warnings = true;

    //!\brief Options that are passed on to the internal stream oject.
    transparent_ostream_options stream_options{};

    /*!\brief Write legacy Windows line-endings including carriage return.
     *
     * \details
     *
     * This option results in old Windows-style line-endings ("\r\n"). Since Windows supports the typical UNIX
     * line-endigns ("\n") nowadays, this option is is highly discouraged.
     *
     * Binary formats always ignore this option.
     */
    bool windows_eol = false;

    /*!\brief Write IDX fields in the header.
     *
     * \details
     *
     * TODO
     */
    bool write_IDX = false;

    // TODO static_assert
};

} // namespace bio::var_io
