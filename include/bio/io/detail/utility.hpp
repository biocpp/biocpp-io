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

#include <bio/meta/tuple.hpp>
#include <bio/meta/type_list/function.hpp>
#include <bio/meta/type_list/type_list.hpp>
#include <bio/meta/type_traits/template_inspection.hpp>

#include <bio/io/exception.hpp>

namespace bio::io::detail
{

/*!\addtogroup io
 * \{
 */

//=============================================================================
// move_tracker
//=============================================================================

//!\brief Can be included as a member to infer whether parent is in moved-from state.
//!\ingroup io
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

//=============================================================================
// copy_blocker
//=============================================================================

//!\brief Can be included as a member to prevent the parent from being copyable/movable.
//!\ingroup io
template <bool shallow>
struct copy_blocker
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    copy_blocker()                            = default; //!<Defaulted.
    copy_blocker(copy_blocker &&)             = default; //!<Defaulted.
    copy_blocker & operator=(copy_blocker &&) = default; //!<Defaulted.

    copy_blocker(copy_blocker const &) requires(!shallow) = default; //!<Conditionally defaulted.
    copy_blocker(copy_blocker const &) requires(shallow)  = delete;  //!<USE DEEP RECORDS IF YOU WANT TO COPY!
    copy_blocker & operator=(copy_blocker const &) requires(!shallow) = default; //!<Conditionally defaulted.
    copy_blocker & operator=(copy_blocker const &) requires(shallow) = delete; //!<USE DEEP RECORDS IF YOU WANT TO COPY!
    //!\}

    //!\brief Comparison operator that always compares as equal.
    friend auto operator<=>(copy_blocker const & lhs, copy_blocker const & rhs) = default;
};

//=============================================================================
// clear()
//=============================================================================

//!\brief Helper for clearing objects that provide such functionality.
//!\ingroup io
void clear(auto && arg) requires(requires { arg.clear(); })
{
    arg.clear();
}

//!\overload
void clear(auto && arg)
{
    arg = {};
}

//=============================================================================
// is_shallow_v
//=============================================================================

//!\brief A type that is a an lvalue reference or a copyable view.
template <typename t>
inline constexpr bool is_shallow_v = std::is_lvalue_reference_v<t> ||
                                     (std::ranges::view<t> && std::copy_constructible<t>);

//!\brief A type that is a an lvalue reference or a copyable view. [specialisation for std::tuple]
template <typename... args_t>
inline constexpr bool is_shallow_v<std::tuple<args_t...>> = (is_shallow_v<args_t> || ...);

//!\brief A type that is a an lvalue reference or a copyable view. [specialisation for std::pair]
template <typename... args_t>
inline constexpr bool is_shallow_v<std::pair<args_t...>> = (is_shallow_v<args_t> || ...);

//!\brief A type that is a an lvalue reference or a copyable view. [specialisation for bio::meta::tuple]
template <typename... args_t>
inline constexpr bool is_shallow_v<bio::meta::tuple<args_t...>> = (is_shallow_v<args_t> || ...);

//!\brief A type that is a an lvalue reference or a copyable view. [specialisation for containers]
template <std::ranges::range container_t>
    requires(!std::ranges::view<container_t>)
inline constexpr bool is_shallow_v<container_t> = is_shallow_v<std::ranges::range_value_t<container_t>>;

} // namespace bio::io::detail
