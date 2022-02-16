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

//!\brief Wrapper to create an overload set of multiple functors.
template <typename... functors>
struct overloaded : functors...
{
    using functors::operator()...;
};

//!\brief Deduction guide for bio::detail::overloaded.
template <typename... functors>
overloaded(functors...) -> overloaded<functors...>;
//!\}

//!\brief Can be included as a member to infer whether parent is in moved-from state.
//!\ingroup bio
struct move_tracker
{
    //!\brief Defaulted.
    move_tracker() = default;

    //!\brief Sets moved-from state.
    move_tracker(move_tracker && rhs) { rhs.moved_from = true; }
    //!\brief Sets moved-from state.
    move_tracker & operator=(move_tracker && rhs)
    {
        rhs.moved_from = true;
        return *this;
    }

    move_tracker(move_tracker const & rhs)             = default; //!< Defaulted.
    move_tracker & operator=(move_tracker const & rhs) = default; //!< Defaulted.

    //!\brief The state.
    bool moved_from = false;
};

} // namespace bio::detail
