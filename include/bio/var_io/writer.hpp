// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::var_io::writer.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <iosfwd>

#include <bio/detail/writer_base.hpp>
#include <bio/format/vcf_output_handler.hpp>
#include <bio/var_io/header.hpp>
#include <bio/var_io/writer_options.hpp>
#include <seqan3/std/filesystem>

namespace bio::var_io
{

// ----------------------------------------------------------------------------
// writer
// ----------------------------------------------------------------------------

/*!\brief A class for writing variant files, e.g. VCF, BCF, GVCF.
 * \tparam options_t A specialisation of bio::var_io::writer_options.
 * \ingroup var_io
 *
 * \details
 *
 * TODO
 */
template <typename... option_args_t>
class writer : public writer_base<writer_options<option_args_t...>>
{
private:
    using base_t = writer_base<writer_options<option_args_t...>>;
    using typename base_t::format_type;

    using base_t::format_handler;
    using base_t::init_state;

public:
    // need these for CTAD
    explicit writer(std::filesystem::path const &            filename,
                    writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{filename, opt}
    {}

    writer(std::filesystem::path const &            filename,
           format_type const &                      fmt,
           writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{filename, fmt, opt}
    {}

    writer(std::ostream &                           stream,
           format_type const &                      fmt,
           writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{stream, fmt, opt}
    {}

    template <typename temporary_stream_t>
        requires(movable_ostream<temporary_stream_t> && !std::is_lvalue_reference_v<temporary_stream_t>)
    writer(temporary_stream_t &&                    stream,
           format_type const &                      fmt,
           writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{std::move(stream), fmt, opt}
    {}

    //!\brief Get the header used by the format.
    bio::var_io::header const & header()
    {
        return std::visit([](auto const & handler) { return handler.get_header(); }, format_handler);
    }

    //!\brief Set the header to the given value.
    template <typename header_t>
        requires std::same_as<bio::var_io::header, std::remove_cvref_t<header_t>>
    void set_header(header_t && hdr)
    {
        if (!init_state)
            throw std::logic_error{"You cannot change the header after I/O has happened."};

        std::visit([&hdr](auto & handler) { handler.set_header(std::forward<header_t>(hdr)); }, format_handler);
    }
};

} // namespace bio::var_io
