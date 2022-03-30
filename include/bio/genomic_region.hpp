// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2022, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <cassert>
#include <limits>

#include <bio/misc.hpp>

/*!\file
 * \brief Provides bio::genomic_region.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

namespace bio
{

/*!\brief Represents an interval on a chromosome or contig.
 * \tparam own Ownership of the #chrom member. See \ref shallow_vs_deep for more details.
 * \ingroup bio
 * \details
 *
 * This structure can be used to describe regions on chromosomes, contigs or other sequences.
 *
 * All member functions assume half-open intervals, and all member functions assume that
 * the invariant `end >= beg` holds. The results are unspecified if either is not true.
 *
 * ### Interval notations
 *
 * 0-based half-open (`[beg, end)`) intervals are the default in this library, but in some
 * places 1-based, closed intervals (`[beg, end]`) may be represented in this data structure. They will be explicitly
 * noted as such.
 */
template <ownership own = ownership::deep>
struct genomic_region
{
    //!\brief Type of the #chrom member.
    using string_t = std::conditional_t<own == ownership::deep, std::string, std::string_view>;

    //!\brief The chromosome or contig identifier .
    string_t chrom;
    //!\brief Beginning of the interval.
    int64_t  beg = 0;
    //!\brief End of the interval; must be >= #beg.
    int64_t  end = std::numeric_limits<int64_t>::max();

    //!\brief Defaulted comparison operators.
    friend bool operator==(genomic_region const &, genomic_region const &)  = default;
    //!\brief Defaulted comparison operators.
    friend auto operator<=>(genomic_region const &, genomic_region const &) = default;

    /*!\brief Checks whether this region lies before, over or beyond the given point.
     * \param chrom The contig/chromosome string for the query position.
     * \param pos   The query position.
     * \returns std::weak_ordering denoting the relative position of the point.
     * \details
     *
     * Returns:
     *
     *   * std::weak_ordering::less if this region ends before the point.
     *   * std::weak_ordering::equivalent if the point is inside the region.
     *   * std::weak_ordering::greater if this region begins beyond the point.
     *
     * The strings are compared lexicographically. The interval is assumed to be half-open, i.e. \p pos == #end will
     * result in std::weak_ordering::greater.
     */
    std::weak_ordering relative_to(std::string_view const chrom, int64_t const pos) const
    {
        assert(beg <= end);
        return std::tuple{std::string_view{this->chrom}, (pos < end)} <=> std::tuple{chrom, (pos >= beg)};
    }

    /*!\brief Checks whether this region lies before or beyond the given region or whether they overlap.
     * \param rhs The region to compare.
     * \returns std::weak_ordering denoting the relative position of the this region compared to rhs.
     * \details
     *
     * Returns:
     *
     *   * std::weak_ordering::less if this region ends before \p rhs begins.
     *   * std::weak_ordering::equivalent if they overlap.
     *   * std::weak_ordering::greater if this region begins before after \p rhs ends.
     *
     * The strings are compared lexicographically. The interval is assumed to be half-open.
     */
    template <ownership own_ = ownership::deep>
    std::weak_ordering relative_to(genomic_region<own_> const & rhs) const
    {
        assert(beg <= end);
        assert(rhs.beg <= rhs.end);
        return std::tuple{std::string_view{chrom}, (end > rhs.beg)} <=>
               std::tuple{std::string_view{rhs.chrom}, (beg < rhs.end)};
    }

    /*!\brief Computes the distance/overlap of two regions.
     * \param rhs The region to compare with.
     * \returns See below.
     * \details
     *
     * Returns:
     *
     *  * A negative value (the overlap) if the regions overlap.
     *  * A positive value (the distance) if the regions do not overlap but are on the same chromosome.
     *  * std::numeric_limits<int64_t>::max() if the regions have different #chrom values.
     *
     * The strings are compared lexicographically for identity. The interval is assumed to be half-open.
     */
    template <ownership own_ = ownership::deep>
    int64_t distance(genomic_region<own_> const & rhs) const
    {
        assert(beg <= end);
        assert(rhs.beg <= rhs.end);
        return chrom == rhs.chrom ? std::max(beg, rhs.beg) - std::min(end, rhs.end)
                                  : std::numeric_limits<int64_t>::max();
    }
};

} // namespace bio
