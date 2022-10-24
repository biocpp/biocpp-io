// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::var::writer_options.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bio/meta/tag/ttag.hpp>

#include <bio/io/format/bcf.hpp>
#include <bio/io/format/vcf.hpp>
#include <bio/io/stream/transparent_ostream.hpp>
#include <bio/io/var/misc.hpp>
#include <bio/io/var/record.hpp>

namespace bio::io::var
{

/*!\brief Options that can be used to configure the behaviour of bio::io::var::writer.
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup var
 *
 * \details
 *
 * TODO describe how to easily initialise this
 */
template <typename formats_t = meta::type_list<bcf, vcf>>
struct writer_options
{
    /*!\brief Try to use types smaller than 32bit to represent integers.
     *
     * \details
     *
     * **BCF-only**
     *
     * TODO
     *
     */
    bool compress_integers = true;

    /*!\brief The formats that output files can take; a bio::meta::ttag over the types.
     *
     * \details
     *
     * See bio::io::var::writer for an overview of the the supported formats.
     */
    formats_t formats = meta::ttag<bcf, vcf>;

    //!\brief Options that are passed on to the internal stream oject.
    transparent_ostream_options stream_options{};

    /*!\brief Write legacy Windows line-endings including carriage return.
     *
     * \details
     *
     * **VCF-only**
     *
     * This option results in old Windows-style line-endings ("\r\n"). Since Windows supports the typical UNIX
     * line-endigns ("\n") nowadays, this option is is highly discouraged.
     */
    bool windows_eol = false;

    /*!\brief Write IDX fields in the header.
     *
     * \details
     *
     * According to the specification for VCF/BCF, entries in the header may be given an IDX value
     * which uniquely identifies that entry¹.
     * Similar to bcftools, the BIO implementation for BCF always writes IDX values.
     * For VCF the default is to not write them (although they are always used internally).
     * This switch turns on writing for VCF, too.
     *
     * ¹ There are two sets of IDX values: one for contigs and one for INFO, FILTER and FORMAT entries (combined).
     *
     * This option is always assumed to be true for bio::io::bcf.
     */
    bool write_IDX = false;

    /*!\brief Verify types when writing.
     *
     * \details
     *
     * **BCF-only**
     *
     * By default this implementation takes your data and transforms into the respective BCF encoding, e.g.
     * a vector of floats (or doubles) is always written as a vector of floats. This is independent of what
     * the header says the field should be.
     *
     * This option activates a check that verifies compatibility of the type information in the header
     * and the data you provide.
     */
    bool verify_header_types = false;

private:
    static_assert(meta::detail::is_type_list<formats_t>, "formats must be a bio::meta::ttag / bio::meta::type_list.");
};

} // namespace bio::io::var
