// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides miscellaneous content for sequence IO.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bio/misc.hpp>

#include <bio/record.hpp>

namespace bio::seq_io
{

//!\brief Default fields for seqan3::seq_io::reader_options.
//!\ingroup seq_io
inline constinit auto default_field_ids = vtag<field::id, field::seq, field::qual>;

} // namespace bio::seq_io
