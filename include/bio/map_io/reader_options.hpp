// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::map_io::reader_options and various pre-defined field_types.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/utility/type_list/traits.hpp>

#include <bio/format/sam.hpp>
#include <bio/map_io/header.hpp>
#include <bio/map_io/misc.hpp>
#include <bio/map_io/sam_flag.hpp>
#include <bio/map_io/sam_tag_dictionary.hpp>
#include <bio/stream/transparent_istream.hpp>
#include <bio/stream/transparent_ostream.hpp>

namespace bio::map_io
{

/*!\name Pre-defined field types
 * \brief These can be used to configure the behaviour of the bio::map_io::reader va bio::map_io::reader_options.
 * \{
 */

/*!\brief The default field types for variant io.
 *!\ingroup map_io
 *
 * \details
 *
 * These traits define a record type with minimal memory allocations for all input formats.
 * It is the recommended record type when iterating ("streaming") over files that ca be any variant IO format.
 *
 * The "style" of the record resembles the BCF specification, i.e. contigs, FILTERs and INFO identifiers are
 * represented as numbers (not strings); and the genotypes are encoded by-genotype (not by-sample).
 * See bio::map_io::genotypes_bcf_style for more information on the latter.
 *
 * \warning Shallow types
 *
 * These records are not self-contained, i.e. they depend on caches and will become invalid when the reader moves to
 * the next record.
 * Since some elements in the record are views, it may not be possible and/or safe to change all values.
 */
template <ownership own = ownership::shallow>
inline constexpr auto field_types =
  ttag<std::string_view,                                                       // field::qname,
       sam_flag,                                                               // field::flag,
       std::string_view,                                                       // field::rname,
       int32_t,                                                                // field::pos,
       uint8_t,                                                                // field::mapq,
       std::vector<seqan3::cigar>,                                             // field::cigar,
       std::string_view,                                                       // field::rnext,
       int32_t,                                                                // field::pnext,
       int32_t,                                                                // field::tlen,
       decltype(std::string_view{} | seqan3::views::char_to<seqan3::dna5>),    // field::seq
       decltype(std::string_view{} | seqan3::views::char_to<seqan3::phred42>), // field::qual
       sam_tag_dictionary,                                                     // field::tags
       record_private_data>;                                                   // field::_private
//!\}

/*!\brief Options that can be used to configure the behaviour of bio::map_io::reader.
 * \tparam field_ids_t   Type of the field_ids member (usually deduced).
 * \tparam field_types_t Type of the field_types member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup map_io
 *
 * \details
 *
 * TODO describe how to easily initialise this
 */
template <typename field_ids_t   = std::remove_cvref_t<decltype(default_field_ids)>,
          typename field_types_t = std::remove_cvref_t<decltype(field_types<ownership::shallow>)>,
          typename formats_t     = seqan3::type_list<sam>>
struct reader_options
{
    //!\brief The fields that shall be contained in each record; a seqan3::tag over seqan3::field.
    field_ids_t field_ids{};

    /*!\brief The types corresponding to each field; a seqan3::type_tag over the types.
     *
     * \details
     *
     * See bio::map_io::reader for an overview of the supported field/type combinations.
     */
    field_types_t field_types{};

    /*!\brief The formats that input files can take; a seqan3::type_tag over the types.
     *
     * \details
     *
     * See bio::map_io::reader for an overview of the the supported formats.
     */
    formats_t formats{};

    //!\brief Whether to print non-critical file format warnings.
    bool print_warnings = true;

    //!\brief Options that are passed on to the internal stream oject.
    transparent_istream_options stream_options{};

    // TODO static_assert
};

} // namespace bio::map_io
