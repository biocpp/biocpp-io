// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides miscellaneous utilities.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <algorithm>
#include <concepts>
#include <filesystem>
#include <ranges>
#include <string>

#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/type_list/detail/type_list_algorithm.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

#include <bio/exception.hpp>

namespace bio::detail
{

/*!\addtogroup bio
 * \{
 */

/*!\brief Sets the file format according to the file name extension.
 * \param[out] format    The format to set.
 * \param[in]  file_name The file name to extract the extension from.
 *
 * \throws seqan3::unhandled_extension_error If the extension in file_name does
 *         not occur in any valid extensions of the formats specified in the
 *         \p format_variant_type template argument list.
 */
void set_format(auto & format, std::filesystem::path const & file_name)
{
    using format_variant_type = std::remove_cvref_t<decltype(format)>;
    using valid_formats       = seqan3::detail::transfer_template_args_onto_t<format_variant_type, seqan3::type_list>;

    bool        format_found = false;
    std::string extension    = file_name.extension().string();
    if (extension.size() > 1)
    {
        extension = extension.substr(1); // drop leading "."
        seqan3::detail::for_each<valid_formats>(
          [&](auto fmt)
          {
              using fm_type = typename decltype(fmt)::type; // remove type_identity wrapper

              for (auto const & ext : fm_type::file_extensions)
              {
                  if (std::ranges::equal(ext, extension))
                  {
                      format.template emplace<fm_type>();
                      format_found = true;
                      return;
                  }
              }
          });
    }

    if (!format_found)
        throw unhandled_extension_error("No valid format found for this extension.");
}

//!\brief Wrapper to create an overload set of multiple functors.
template <typename... functors>
struct overloaded : functors...
{
    using functors::operator()...;
};

//!\brief Deduction guide for bio::detail::overloaded.
template <typename... functors>
overloaded(functors...) -> overloaded<functors...>;

/*!\brief Pass this function a constrained functor that accepts one argument and returns std::true_type.
 * \details
 *
 * See e.g. bio::seq_io::reader_options to see how this is used.
 */
constexpr bool lazy_concept_checker(auto fun)
{
    auto fallback = []<typename T = int>(auto) { return std::false_type{}; };
    using ret_t   = decltype(detail::overloaded{fallback, fun}(1));
    return ret_t::value;
}

//!\brief Returns the size of the argument, either a range or a tuple.
constexpr size_t range_or_tuple_size(std::ranges::forward_range auto && r)
{
    return std::ranges::distance(r);
}

//!\overload
template <typename tuple_t>
    requires requires { typename std::tuple_size<std::remove_cvref_t<tuple_t>>::type; }
constexpr size_t range_or_tuple_size(tuple_t)
{
    return std::tuple_size_v<tuple_t>;
}
//!\}

} // namespace bio::detail
