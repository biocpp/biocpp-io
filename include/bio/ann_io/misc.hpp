// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::var_io::tag_dictionary class and auxiliaries.
 * \author Joshua Kim <kim_j AT molgen.mpg.de>
 */

 #include <bio/misc.hpp>
 #include <bio/record.hpp>

namespace bio::ann_io
{
//-----------------------------------------------------------------------------
// default_field_ids
//-----------------------------------------------------------------------------

//!\brief Default fields for bio::var_io::reader_options.
//!\ingroup var_io
inline constinit auto default_field_ids = vtag<field::chrom,
                                               field::chromStart,
                                               field::chromEnd>;

//-----------------------------------------------------------------------------
// Pre-defined field types (reader)
//-----------------------------------------------------------------------------

/*!\name Pre-defined field types
* \brief These can be used to configure the behaviour of the bio::var_io::reader via bio::var_io::reader_options.
* \{
*/
/*!\brief The default field types for variant io.
*!\ingroup var_io
*
* \details
*
* These traits define a record type with minimal memory allocations for all input formats.
* It is the recommended record type when iterating ("streaming") over files that ca be any variant IO format.
*
* The "style" of the record resembles the VCF specification, i.e. contigs, FILTERs and INFO identifiers are
* represented as string/string_views. **However,** the genotypes are encoded by-genotype (BCF-style) and not by-sample
* (VCF-style) for performance reasons.
*
* See bio::var_io::genotypes_bcf_style for more information on the latter.
*/
template <ownership own = ownership::shallow>
inline constinit auto field_types =
ttag<std::string_view,                          // field::chrom,
     uint32_t,                                  // field::chromStart,
     uint32_t>;                                 // field::chromEnd

//!\brief Deep version of bio::var_io::field_types.
//!\ingroup var_io
template <>
inline constinit auto field_types<ownership::deep> =
ttag<std::string,                                     // field::chrom,
     uint32_t,                                        // field::chromStart,
     uint32_t>;                                       // field::chromEnd

}
