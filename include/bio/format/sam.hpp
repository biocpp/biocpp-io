// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the bio::format_sam.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <string>
#include <vector>

#include <bio/platform.hpp>

namespace bio
{

/*!\brief       The sam format.
 * \ingroup     format
 *
 * \details
 *
 * This is the sam format tag. If you want to read sam files, use bio::map_io::reader, and if you want
 * to write sam files, use bio::map_io::writer.
 *
 * ### Introduction
 *
 * sam is the de-facto-standard for mapping storage in bionformatics. See the
 * [article on wikipedia](todo) for a an in-depth description of the format.
 *
 * ### Fields
 *
 * todo
 *
 * ### Implementation notes
 *
 * todo
 *
 */
struct sam
{
    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions{{"sam"}};
};

} // namespace bio
