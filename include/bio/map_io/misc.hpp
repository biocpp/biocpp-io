// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides miscellaneous content for sequence IO.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <bio/misc.hpp>

#include <bio/record.hpp>

namespace bio::map_io
{

//!\brief Default fields for bio::map_io::reader_options.
//!\ingroup map_io
inline constexpr auto default_field_ids = vtag<field::qname,
                                               field::flag,
                                               field::rname,
                                               field::pos,
                                               field::mapq,
                                               field::cigar,
                                               field::rnext,
                                               field::pnext,
                                               field::tlen,
                                               field::seq,
                                               field::qual,
                                               field::tags,
                                               field::_private>;

} // namespace bio::map_io
