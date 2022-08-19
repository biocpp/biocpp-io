// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <concepts>
#include <tuple>

#include <bio/meta/type_list/type_list.hpp>

#include <bio/io/platform.hpp>

/*!\file
 * \brief Provides various minor utilities.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

namespace bio::io
{

//-----------------------------------------------------------------------------
// ownership
//-----------------------------------------------------------------------------

/*!\brief An enum used as an argument for templates that switch between owning and non-owning behaviour.
 * \ingroup io
 * \details
 *
 * Typically used to configure a class template to have members that are vectors/strings VS members that are views.
 * The "shallow" version of such a class is typically cheap to copy (no dynamic memory) while the "deep" version
 * is expensive to copy (holds dynamic memory).
 *
 * See \ref shallow_vs_deep on what this means in practice.
 */
enum class ownership
{
    shallow, //!< Cheap to copy.
    deep     //!< Expensive to copy.
};

} // namespace bio::io
