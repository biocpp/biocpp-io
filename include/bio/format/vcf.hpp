// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::vcf.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/platform.hpp>

namespace bio
{

/*!\brief The Variant Call Format (VCF) format.
 * \ingroup format
 *
 * \details
 *
 * This is the VCF format tag. If you want to read VCF files, use bio::var_io::reader, and if you want
 * to write VCF files, use bio::var_io::writer.
 *
 * ### Fields
 *
 * The format consists of the following fields:
 *
 *   1. bio::field::chrom
 *   2. bio::field::pos
 *   3. bio::field::id
 *   4. bio::field::ref
 *   5. bio::field::alt
 *   6. bio::field::qual
 *   7. bio::field::filter
 *   8. bio::field::info
 *   9. bio::field::genotypes
 *
 * See bio::var_io::reader and bio::var_io::writer for more details.
 *
 * ### Implementation
 *
 * The implementation targets [version 4.3 of the VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf).
 * However, reading version 4.2 should be possible, too.
 * HTSLib (BCFTools) deviates from the specification in several places. Because it is the de-facto standard, this
 * implementation emulates HTSLib behaviour whenever necessary (even if that means violating the spec).
 *
 * Little testing has been done on handling structural variants and breakend strings, but in theory the values
 * should be parsed correctly (as strings). No testing has been done on gVCF files, but in theory all values should
 * be parsed correctly.
 *
 *
 * Please report any issues you find.
 */
struct vcf
{
    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions{{"vcf"}};
};

} // namespace bio
