// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2022, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the utilities for Tabix indexes.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <bio/detail/views_eager_split.hpp>
#include <bio/genomic_region.hpp>
#include <bio/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/stream/transparent_istream.hpp>
#include <bio/stream/transparent_ostream.hpp>

namespace bio::detail
{

//!\brief Converts a BGZip "virtual offset" into the on-disk offset to the beginning of the GZ-block and the distance
//! within that block after decompression.
constexpr std::pair<uint64_t, uint16_t> decode_bgz_virtual_offset(uint64_t const in)
{
    uint64_t const comp_offset   = in >> 16;
    uint16_t const uncomp_offset = static_cast<uint16_t>(in);

    return {comp_offset, uncomp_offset};
}

//!\brief Tabix index support
struct tabix_index
{
    //!\brief Magic bytes of the TABIX format (after decompression).
    static constexpr char magic_bytes[4] = {'T', 'B', 'I', 1};

    //!\brief The "header" or core data members.
    struct core_t
    {
        char    magic[4] = {}; //!< Magic string.
        int32_t n_ref    = 0;  //!< Number of reference sequences/indexes.
        int32_t format   = 0;  //!< Format.
        int32_t col_seq  = 0;  //!< Column for the sequence name.
        int32_t col_beg  = 0;  //!< Column for the start of a region.
        int32_t col_end  = 0;  //!< Column for the end of a region.
        int32_t meta     = 0;  //!< Leading character for comment lines.
        int32_t skip     = 0;  //!< # lines to skip at the beginning.
        int32_t l_nm     = 0;  //!< Length of concatenated sequence names.
    };

    //!\brief Chunk (n per bin).
    struct chunk_t
    {
        uint64_t cnk_beg = 0; //!< Begin position.
        uint64_t cnk_end = 0; //!< End position.
    };

    //!\brief Bin (n per index).
    struct bin_t
    {
        uint32_t             bin = 0; //!< Bin identifier.
        std::vector<chunk_t> chunks;  //!< Vector of chunks.
    };

    //!\brief Index (one per reference).
    struct index_t
    {
        std::vector<bin_t>                          bins;    //!< Bins.
        std::unordered_map<uint32_t, bin_t const *> bin_map; //!< Map of Bin-ID to Bin.
        std::vector<uint64_t>                       offsets; //!< Offsets for linear scan.
    };

    core_t                                       core;      //!< Core Fields.
    std::vector<std::string>                     names;     //!< Reference sequence names.
    std::unordered_map<std::string_view, size_t> names_map; //!< Map of reference name to index number.
    std::vector<index_t>                         indexes;   //!< Indexes (one per reference).

    std::optional<uint64_t> n_no_coor{}; //!< Number of unmapped reads (optional).

    //!\brief Read an index from disk.
    void read(std::filesystem::path const & path);
    //!\brief Write an index to disk.
    void write(std::filesystem::path const & path);

    //!\brief List of theoretically overlapping bin-numbers (independent of this data structure).
    static inline void reg2bins(uint32_t beg, uint32_t end, std::vector<uint32_t> & bin_numbers)
    {
        if (beg >= end)
            return;
        if (end >= 1u << 29)
            end = 1u << 29;
        --end;

        uint32_t k = 0;
        bin_numbers.push_back(k);

        for (k = 1 + (beg >> 26); k <= 1 + (end >> 26); ++k)
            bin_numbers.push_back(k);
        for (k = 9 + (beg >> 23); k <= 9 + (end >> 23); ++k)
            bin_numbers.push_back(k);
        for (k = 73 + (beg >> 20); k <= 73 + (end >> 20); ++k)
            bin_numbers.push_back(k);
        for (k = 585 + (beg >> 17); k <= 585 + (end >> 17); ++k)
            bin_numbers.push_back(k);
        for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k)
            bin_numbers.push_back(k);
    }

    //!\brief Create the list of chunks that potentially overlap the desired region.
    template <ownership own>
    std::vector<std::pair<uint64_t, uint64_t>> reg2chunks(genomic_region<own> const & reg);
};

inline void tabix_index::read(std::filesystem::path const & path)
{
    transparent_istream                    istream{path};
    detail::fast_istreambuf_iterator<char> it{istream};
    unexpected_end_of_input                in_end{"Unexpected end of input while trying to read Tabix index."};

    /* read core */
    it.read_as_binary(core);

    if (!std::ranges::equal(core.magic, magic_bytes))
        throw format_error{"This is not a tabix index."};

    // TODO sanity check the core values

    /* read names */
    std::string names_buffer;
    names_buffer.resize(core.l_nm - 1); // don't reading trailing \0
    it.read_n_chars_into(core.l_nm - 1, names_buffer.data());
    for (std::string_view name : names_buffer | detail::eager_split('\0'))
        names.push_back(static_cast<std::string>(name));
    ++it; // skip trailing \0

    size_t c = 0;
    for (std::string_view const name : names)
        names_map[name] = c++;

    /* read indexes*/
    indexes.resize(core.n_ref);
    for (int32_t i = 0; i < core.n_ref; ++i)
    {
        int32_t n_bin = 0;
        it.read_as_binary(n_bin);

        index_t & index = indexes[i];
        index.bins.resize(n_bin);

        for (int32_t j = 0; j < n_bin; ++j)
        {
            bin_t & bin = index.bins[j];
            it.read_as_binary(bin.bin);

            int32_t n_chunk = 0;
            it.read_as_binary(n_chunk);
            bin.chunks.resize(n_chunk);
            it.read_n_chars_into(sizeof(chunk_t) * n_chunk, reinterpret_cast<char *>(bin.chunks.data()));
        }

        for (bin_t const & bin : index.bins)
            index.bin_map[bin.bin] = &bin;

        int32_t n_intv = 0;
        it.read_as_binary(n_intv);
        index.offsets.resize(n_intv);
        it.read_n_chars_into(sizeof(uint64_t) * n_intv, reinterpret_cast<char *>(index.offsets.data()));
    }

    if (it != std::default_sentinel)
    {
        n_no_coor = uint64_t{};
        it.read_as_binary(*n_no_coor);
    }
}

inline void tabix_index::write(std::filesystem::path const & path)
{
    transparent_ostream                    ostream{path, {.compression = compression_format::bgzf}};
    detail::fast_ostreambuf_iterator<char> it{ostream};
    unexpected_end_of_input                in_end{"Unexpected end of input while trying to read Tabix index."};

    /* write core */
    it.write_as_binary(core);

    /* write names */
    std::string names_buffer;
    for (std::string_view name : names)
    {
        it.write_range(name);
        it = '\0';
    }

    /* write indexes*/
    for (index_t const & i : indexes)
    {
        int32_t n_bin = i.bins.size();
        it.write_as_binary(n_bin);

        for (bin_t const & b : i.bins)
        {
            it.write_as_binary(b.bin);

            int32_t n_chunk = b.chunks.size();
            it.write_as_binary(n_chunk);

            std::string_view chunk_data{reinterpret_cast<char const *>(b.chunks.data()),
                                        sizeof(chunk_t) * b.chunks.size()};
            it.write_range(chunk_data);
        }

        int32_t n_intv = i.offsets.size();
        it.write_as_binary(n_intv);

        std::string_view offsets_data{reinterpret_cast<char const *>(i.offsets.data()),
                                      sizeof(uint64_t) * i.offsets.size()};
        it.write_range(offsets_data);
    }

    if (n_no_coor.has_value())
        it.write_as_binary(*n_no_coor);
}

template <ownership own>
inline std::vector<std::pair<uint64_t, uint64_t>> tabix_index::reg2chunks(genomic_region<own> const & reg)
{
    std::vector<std::pair<uint64_t, uint64_t>> ret;

    size_t n_index = 0;
    if (auto it = names_map.find(reg.chrom); it != names_map.end())
        n_index = it->second;
    else
        throw bio_error{"FOO"};

    index_t const & index = indexes[n_index];

    /* linear index evaluation */
    uint64_t virtual_offset_lower_bound = 0;
    size_t   linear_interval_i          = reg.beg >> 14; // 16kb
    if (linear_interval_i >= index.offsets.size())
        return {};
    else
        virtual_offset_lower_bound = index.offsets[linear_interval_i];

    /* binning index evaluation */
    std::vector<uint32_t> bin_numbers;
    reg2bins(reg.beg, reg.end, bin_numbers);

    for (uint32_t const bin_number : bin_numbers)
    {
        if (auto it = index.bin_map.find(bin_number); it != index.bin_map.end())
        {
            bin_t const & b = *it->second;

            for (chunk_t const & c : b.chunks)
                if (c.cnk_end > virtual_offset_lower_bound)
                    ret.emplace_back(c.cnk_beg, c.cnk_end);
        }
    }

    return ret;
}

} // namespace bio::detail
