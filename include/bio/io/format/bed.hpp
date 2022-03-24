// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::bed.
 * \author Joshua Kim <kim_j AT molgen.mpg.de>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/platform.hpp>

namespace bio
{

/*!\brief The Browser Extensible Data (BED) format.
 * \ingroup format
 *
 * \details
 *
 * This is the BED format tag. If you want to read BED files, use bio::ann_io::reader, and if you want
 * to write BED files, use bio::ann_io::writer.
 *
 * ### Fields
 *
 * The format consists of the following fields:
 *
 *   1. bio::field::chrom
 *   2. bio::field::chromStart
 *   3. bio::field::chromEnd
 *
 * See bio::ann_io::reader and bio::ann_io::writer for more details.
 *
 * ### Implementation
 *
 * The implementation target [version 4.3 of the VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf).
 * However, reading version 4.2 should be possible, too.
 * Little testing has been done on handling structural variants and breakend strings, but in theory the values
 * should be parsed correctly (as strings).
 *
 * No testing has been done on gVCF files, but in theory all values should be parsed correctly.
 *
 * Please report any issues you find.
 */
struct bed
{
    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions{{"bed"}};
};

} // namespace bio
