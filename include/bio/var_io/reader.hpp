// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::var_io::reader.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <filesystem>

#include <bio/detail/reader_base.hpp>
#include <bio/format/bcf_input_handler.hpp>
#include <bio/format/vcf_input_handler.hpp>
#include <bio/var_io/header.hpp>
#include <bio/var_io/reader_options.hpp>

namespace bio::var_io
{

// ----------------------------------------------------------------------------
// reader
// ----------------------------------------------------------------------------

/*!\brief A class for reading variant files, e.g. VCF, BCF, GVCF.
 * \tparam options_t A specialisation of bio::var_io::reader_options.
 * \ingroup var_io
 *
 * \details
 *
 * ### Introduction
 *
 * Variant files are files that contain sequence variation information. Well-known formats include
 * VCF and BCF.
 *
 * The Variant I/O reader supports reading the following fields:
 *
 *   1. bio::field::chrom
 *   2. bio::field::pos
 *   3. bio::field::id
 *   4. bio::field::ref
 *   5. bio::field::alt
 *   6. bio::field::qual
 *   7. bio::field::filter
 *   8. bio::field::info
 *   9. bio::field::genotypes
 *
 * These fields correspond to the order and names defined in the VCF specification. The types and values that
 * are returned by default also correspond to VCF specification (i.e. 1-based positions, string as strings and not
 * as numbers) **with one exception:** the genotypes are not grouped by sample (as in the VCF format) but by
 * genotype field (as in the BCF format).
 * This results in  a notably better performance when reading BCF files. See below for information on how to change
 * this.
 *
 * This reader supports the following formats:
 *
 *   1. VCF (see also bio::vcf)
 *   2. BCF (see also bio::bcf)
 *
 * If you only need to read VCF and not BCF and you do not want to parse the fields into high-level data structures
 * (and simply use them as strings), you can use bio::plain_io::reader instead of this reader.
 *
 * ### Simple usage
 *
 * Iterate over a variant file via the reader and print "CHROM:POS:REF:ALT" for each record:
 *
 * \snippet test/snippet/var_io/var_io_reader.cpp simple_usage_file
 *
 * Read from standard input instead of a file:
 *
 * \snippet test/snippet/var_io/var_io_reader.cpp simple_usage_stream
 *
 * ### Accessing more complex fields
 *
 * TODO
 *
 * ### Views on readers
 *
 * Print information for the first five records where quality is better than 23:
 *
 * \snippet test/snippet/var_io/var_io_reader.cpp views
 *
 * ### Specifying options
 *
 * TODO
 *
 * For more advanced options, see bio::var_io::reader_options.
 */
template <typename... option_args_t>
class reader : public reader_base<reader_options<option_args_t...>>
{
private:
    //!\brief The base class.
    using base_t      = reader_base<reader_options<option_args_t...>>;
    //!\brief Inherit the format_type definition.
    using format_type = typename base_t::format_type;
    /* Implementation note
     * format_type is "inherited" as private here to avoid appearing twice in the documentation.
     * Its actual visibility is public because it is public in the base class.
     */
    //!\brief Make the format handler visible.
    using base_t::format_handler;

    //!\brief A pointer to the header inside the format.
    bio::var_io::header const * header_ptr = nullptr;

public:
    // clang-format off
    //!\copydoc bio::reader_base::reader_base(std::filesystem::path const & filename, format_type const & fmt, options_t const & opt = options_t{})
    // clang-format on
    reader(std::filesystem::path const &            filename,
           format_type const &                      fmt,
           reader_options<option_args_t...> const & opt = reader_options<option_args_t...>{}) :
      base_t{filename, fmt, opt}
    {}

    //!\overload
    explicit reader(std::filesystem::path const &            filename,
                    reader_options<option_args_t...> const & opt = reader_options<option_args_t...>{}) :
      base_t{filename, opt}
    {}

    // clang-format off
    //!\copydoc bio::reader_base::reader_base(std::istream & str, format_type const & fmt, options_t const & opt = options_t{})
    // clang-format on
    reader(std::istream &                           str,
           format_type const &                      fmt,
           reader_options<option_args_t...> const & opt = reader_options<option_args_t...>{}) :
      base_t{str, fmt, opt}
    {}

    //!\overload
    template <movable_istream temporary_stream_t>
        //!\cond REQ
        requires(!std::is_lvalue_reference_v<temporary_stream_t>)
    //!\endcond
    reader(temporary_stream_t &&                    str,
           format_type const &                      fmt,
           reader_options<option_args_t...> const & opt = reader_options<option_args_t...>{}) :
      base_t{std::move(str), fmt, opt}
    {}

    //!\brief Access the header.
    bio::var_io::header const & header()
    {
        if (header_ptr == nullptr)
        {
            // ensure that the format_handler is created
            this->begin();

            header_ptr = std::visit([](auto const & handler) { return &handler.get_header(); }, format_handler);
        }

        return *header_ptr;
    }
};

} // namespace bio::var_io
