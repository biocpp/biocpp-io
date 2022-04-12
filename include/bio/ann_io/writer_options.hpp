// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::ann_io::writer_options.
 * \author Joshua Kim <kim_j AT molgen.mpg.de>
 */

#pragma once

#include <bio/format/bed.hpp>
#include <bio/stream/transparent_ostream.hpp>
#include <bio/ann_io/misc.hpp>

// namespace bio::detail
// {
//
// template <typename t>
// concept ann_io_legal_type_aux =
//   std::same_as<t, char> || std::signed_integral<t> || std::floating_point<t> || std::same_as < std::decay_t<t>,
// char const * > ;
//
// /*!\interface bio::detail::ann_io_legal_type <>
//  * \tparam t The type to check.
//  * \brief A type that is similar to one of the alternatives of bio::ann_io::info_element_value_type
//  */
// //!\cond CONCEPT_DEF
// template <typename t>
// concept ann_io_legal_type = ann_io_legal_type_aux<std::remove_cvref_t<t>> || std::same_as<t const &, bool const &> ||
//   (std::ranges::forward_range<t> && (ann_io_legal_type_aux<std::remove_cvref_t<std::ranges::range_reference_t<t>>> ||
//                                      (std::ranges::forward_range<std::ranges::range_reference_t<t>> &&
//                                       std::same_as<char const &, std::ranges::range_reference_t<t> const &>)));
// //!\endcond
//
// /*!\interface bio::detail::ann_io_legal_vector_type <>
//  * \tparam t The type to check.
//  * \brief A type that is similar to one of the alternatives of bio::ann_io::info_element_value_type
//  */
// //!\cond CONCEPT_DEF
// template <typename t>
// concept ann_io_legal_vector_type =
//   std::ranges::forward_range<t> && ann_io_legal_type<std::ranges::range_reference_t<t>> &&
//   !std::same_as<bool const &, std::ranges::range_reference_t<t>>;
// //!\endcond
//
// /*!\interface bio::detail::ann_io_legal_or_dynamic <>
//  * \tparam t The type to check.
//  * \brief A type that is similar to one of the alternatives of bio::ann_io::info_element_value_type
//  */
// //!\cond CONCEPT_DEF
// template <typename t>
// concept ann_io_legal_or_dynamic = ann_io_legal_type<t> || is_info_element_value_type<t>;
// //!\endcond
//
// /*!\interface bio::detail::ann_io_vector_legal_or_dynamic <>
//  * \tparam t The type to check.
//  * \brief A type that is similar to one of the alternatives of bio::ann_io::info_element_value_type
//  */
// //!\cond CONCEPT_DEF
// template <typename t>
// concept ann_io_vector_legal_or_dynamic = ann_io_legal_vector_type<t> || is_genotype_element_value_type<t>;
// //!\endcond
//
// /*!\interface bio::detail::info_element_writer_concept <>
//  * \tparam t The type to check.
//  * \brief Types "similar" to bio::ann_io::info_element / bio::ann_io::info_element_bcf.
//  */
// //!\cond CONCEPT_DEF
// template <typename t>
// concept info_element_writer_concept = detail::decomposable_into_two<t> &&
//   (detail::char_range_or_cstring<detail::first_elem_t<t>> ||
//    std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::ann_io_legal_or_dynamic<detail::second_elem_t<t>>;
// //!\endcond
//
// /*!\interface bio::detail::genotype_writer_concept <>
//  * \tparam t The type to check.
//  * \brief Types "similar" to bio::ann_io::genotype_element / bio::ann_io::genotype_element_bcf.
//  */
// //!\cond CONCEPT_DEF
// template <typename t>
// concept genotype_writer_concept = detail::decomposable_into_two<t> &&
//   (detail::char_range_or_cstring<detail::first_elem_t<t>> ||
//    std::same_as<int32_t, detail::first_elem_t<t>>)&&detail::ann_io_vector_legal_or_dynamic<detail::second_elem_t<t>>;
// //!\endcond
//
// } // namespace bio::detail

namespace bio::ann_io
{

/*!\brief Options that can be used to configure the behaviour of bio::ann_io::writer.
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup ann_io
 *
 * \details
 *
 * TODO describe how to easily initialise this
 */
template <typename formats_t = seqan3::type_list<bed>>
struct writer_options
{
    /*!\brief The formats that output files can take; a bio::ttag over the types.
     *
     * \details
     *
     * See bio::ann_io::writer for an overview of the the supported formats.
     */
    formats_t formats = ttag<bed>;

    //!\brief Options that are passed on to the internal stream oject.
    transparent_ostream_options stream_options{};

    /*!\brief Write legacy Windows line-endings including carriage return.
     *
     * \details
     *
     * This option results in old Windows-style line-endings ("\r\n"). Since Windows supports the typical UNIX
     * line-endigns ("\n") nowadays, this option is is highly discouraged.
     */
    bool windows_eol = false;

private:
    static_assert(detail::is_type_list<formats_t>, "formats must be a bio::ttag / seqan3::type_list.");
};

} // namespace bio::ann_io
