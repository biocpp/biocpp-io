// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::detail::compression_stream trait.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <filesystem>
#include <string_view>
#include <vector>

#include <bio/io/stream/compression.hpp>
#ifdef BIO_HAS_BZIP2
#    include <bio/io/stream/detail/bz2_istream.hpp>
#    include <bio/io/stream/detail/bz2_ostream.hpp>
#endif
#ifdef BIO_HAS_ZLIB
#    include <bio/io/stream/detail/bgzf_istream.hpp>
#    include <bio/io/stream/detail/bgzf_ostream.hpp>
#    include <bio/io/stream/detail/bgzf_stream_util.hpp>
#    include <bio/io/stream/detail/gz_istream.hpp>
#    include <bio/io/stream/detail/gz_ostream.hpp>
#endif

namespace bio::io::detail
{

/*!\brief Stream aliases for the respective compression formats.
 * \tparam format The bio::io::compression_format whose traits are provided.
 * \ingroup io
 *
 * \details
 *
 * Note that this is **not** a customisation point. You may not specialise it.
 */
template <compression_format format>
struct compression_stream
{
    //!\brief The default input stream corresponding to the format.
    using istream = void;

    //!\brief The default output stream corresponding to the format.
    using ostream = void;
};

/*!\brief Traits for the bio::io::compression_format::bgzf.
 * \see bio::io::compression_traits
 */
template <>
struct compression_stream<compression_format::bgzf> : compression_stream<compression_format::none>
{
#ifdef BIO_HAS_ZLIB
    //!\copydoc bio::io::compression_traits<compression_format::none>::basic_istream
    using istream = contrib::bgzf_istream;

    //!\copydoc bio::io::compression_traits<compression_format::none>::basic_ostream
    using ostream = contrib::bgzf_ostream;
#endif
};

/*!\brief Traits for the bio::io::compression_format::gz.
 * \see bio::io::compression_traits
 */
template <>
struct compression_stream<compression_format::gz> : compression_stream<compression_format::none>
{
#ifdef BIO_HAS_ZLIB
    //!\copydoc bio::io::compression_traits<compression_format::none>::basic_istream
    using istream = contrib::gz_istream;

    //!\copydoc bio::io::compression_traits<compression_format::none>::basic_ostream
    using ostream = contrib::gz_ostream;
#endif
};

/*!\brief Traits for the bio::io::compression_format::bz2.
 * \see bio::io::compression_traits
 */
template <>
struct compression_stream<compression_format::bz2> : compression_stream<compression_format::none>
{
#ifdef BIO_HAS_BZIP2
    //!\copydoc bio::io::compression_traits<compression_format::none>::basic_istream
    using istream = contrib::bz2_istream;

    //!\copydoc bio::io::compression_traits<compression_format::none>::basic_ostream
    using ostream = contrib::bz2_ostream;
#endif
};

/*!\brief Traits for the bio::io::compression_format::zstd.
 * \see bio::io::compression_traits
 */
template <>
struct compression_stream<compression_format::zstd> : compression_stream<compression_format::none>
{
    // NOT YET SUPPORTED BY SEQAN
};

//-------------------------------------------------------------------------------
// make_istream
//-------------------------------------------------------------------------------

//!\brief Creates an input stream by forwarding the given arguments.
template <compression_format format, typename... args_t>
std::istream * make_istream(args_t &&... args)
{
    if constexpr (compression_traits<format>::available)
    {
        return new typename compression_stream<format>::istream{std::forward<args_t>(args)...};
    }
    else
    {
        throw file_open_error{"The file is ",
                              compression_traits<format>::as_string,
                              "-compressed, but SeqAn3 is built without support for this format."};
        return nullptr;
    }
}

//-------------------------------------------------------------------------------
// make_ostream
//-------------------------------------------------------------------------------

//!\brief Creates an output stream by forwarding the given arguments.
template <compression_format format, typename... args_t>
std::ostream * make_ostream(args_t &&... args)
{
    if constexpr (compression_traits<format>::available)
    {
        return new typename compression_stream<format>::ostream{std::forward<args_t>(args)...};
    }
    else
    {
        throw file_open_error{compression_traits<format>::as_string,
                              "-compression was selected, but SeqAn3 is built without support for this format."};
        return nullptr;
    }
}

} // namespace bio::io::detail
