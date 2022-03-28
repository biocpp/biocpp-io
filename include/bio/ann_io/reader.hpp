// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::ann_io::reader.
 * \author Joshua Kim <kim_j AT molgen.mpg.de>
 */

 #pragma once

 #include <filesystem>

 #include <bio/detail/reader_base.hpp>
 #include <bio/ann_io/reader_options.hpp>
 #include <bio/ann_io/header.hpp>
 #include <bio/format/bed_input_handler.hpp>

namespace bio::ann_io
{
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
    bio::ann_io::header const * header_ptr = nullptr;
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
    bio::ann_io::header const & header()
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
}
