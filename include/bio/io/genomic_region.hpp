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
#include <string>
#include <tuple>

#include <bio/io.hpp>

/*!\file
 * \brief Provides bio::io::genomic_region.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

namespace bio::io
{

/*!\brief Represents an interval on a chromosome or contig.
 * \tparam own Ownership of the #chrom member. See \ref shallow_vs_deep for more details.
 * \ingroup io
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
struct genomic_region
{
    //!\brief The chromosome or contig identifier.
    std::string chrom;
    //!\brief Beginning of the interval.
    int64_t     beg = 0;
    //!\brief End of the interval; must be >= #beg.
    int64_t     end = std::numeric_limits<int64_t>::max();

    //!\brief Defaulted comparison operators.
    friend bool operator==(genomic_region const &, genomic_region const &)  = default;
    //!\brief Defaulted comparison operators.
    friend auto operator<=>(genomic_region const &, genomic_region const &) = default;

    /*!\brief Checks whether the given region lies before, on or beyond the given point.
     * \param lchrom The contig/chromosome string for the query region.
     * \param lbeg   The query region begin position.
     * \param lend   The query region end position.
     * \param rchrom The contig/chromosome string for the point.
     * \param rpos   The position of the point.
     *
     * \returns std::weak_ordering denoting the relative position of the point.
     * \details
     *
     * Returns:
     *
     *   * std::weak_ordering::less if this region ends before the point.
     *   * std::weak_ordering::equivalent if the point is inside the region.
     *   * std::weak_ordering::greater if this region begins beyond the point.
     *
     * The strings are compared lexicographically. The interval is assumed to be half-open, i.e. \p rpos == \p lend will
     * result in std::weak_ordering::greater.
     */
    static constexpr std::weak_ordering relative_to(std::string_view const lchrom,
                                                    int64_t const          lbeg,
                                                    int64_t const          lend,
                                                    std::string_view const rchrom,
                                                    int64_t const          rpos)
    {
        assert(lbeg <= lend);
        return std::tuple{lchrom, (rpos < lend)} <=> std::tuple{rchrom, (rpos >= lbeg)};
    }

    //!\overload
    static inline std::weak_ordering relative_to(genomic_region const & lhs,
                                                 std::string_view const rchrom,
                                                 int64_t const          rpos)
    {
        return relative_to(lhs.chrom, lhs.beg, lhs.end, rchrom, rpos);
    }

    /*!\brief Checks whether the first region lies before or beyond the second region or whether they overlap.
     * \param lchrom The contig/chromosome string for the left region.
     * \param lbeg   The left region's begin position.
     * \param lend   The left region's end position.
     * \param rchrom The contig/chromosome string for the right region.
     * \param rbeg   The right region's begin position.
     * \param rend   The right region's end position.
     * \returns std::weak_ordering denoting the relative position of the left region compared to right region.
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
    static constexpr std::weak_ordering relative_to(std::string_view const lchrom,
                                                    int64_t const          lbeg,
                                                    int64_t const          lend,
                                                    std::string_view const rchrom,
                                                    int64_t const          rbeg,
                                                    int64_t const          rend)
    {
        assert(lbeg <= lend);
        assert(rbeg <= rend);
        return std::tuple{lchrom, (lend > rbeg)} <=> std::tuple{rchrom, (lbeg < rend)};
    }

    //!\overload
    static inline std::weak_ordering relative_to(genomic_region const & lhs, genomic_region const & rhs)
    {
        return relative_to(lhs.chrom, lhs.beg, lhs.end, rhs.chrom, rhs.beg, rhs.end);
    }

    /*!\brief Computes the distance/overlap of two regions.
     * \param lchrom The contig/chromosome string for the left region.
     * \param lbeg   The left region's begin position.
     * \param lend   The left region's end position.
     * \param rchrom The contig/chromosome string for the right region.
     * \param rbeg   The right region's begin position.
     * \param rend   The right region's end position.
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
    static constexpr int64_t distance(std::string_view const lchrom,
                                      int64_t const          lbeg,
                                      int64_t const          lend,
                                      std::string_view const rchrom,
                                      int64_t const          rbeg,
                                      int64_t const          rend)
    {
        assert(lbeg <= lend);
        assert(rbeg <= rend);
        return lchrom == rchrom ? std::max(lbeg, rbeg) - std::min(lend, rend) : std::numeric_limits<int64_t>::max();
    }

    //!\overload
    static inline int64_t distance(genomic_region const & lhs, genomic_region const & rhs)
    {
        return distance(lhs.chrom, lhs.beg, lhs.end, rhs.chrom, rhs.beg, rhs.end);
    }
};

} // namespace bio::io
