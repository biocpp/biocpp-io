// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::transparent_istream.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <span>
#include <thread>

#include <bio/io/exception.hpp>
#include <bio/io/stream/compression.hpp>
#include <bio/io/stream/concept.hpp>
#include <bio/io/stream/detail/make_stream.hpp>

namespace bio::io
{

//!\brief Options that can be provided to bio::io::transparent_istream.
//!\ingroup stream
struct transparent_istream_options
{
    //!\brief Size of the buffer used when opening a file from a filename.
    size_t buffer1_size = 1024 * 1024;
    //!\brief Size of the buffer used for the decompression stream.
    size_t buffer2_size = 1024 * 1024 * 4;

    /*!\brief Which decompressor to use.
     *
     * \details
     *
     * In almost all cases, you will want to have this auto-detected. You can explicitly set this to
     * bio::io::compression_format::gz when opening BGZF-compressed files to use the regular, single-threaded
     * GZ-decompressor, but it is simpler to just set threads to 1.
     */
    compression_format compression = compression_format::detect;

    /*!\brief Maximum number of threads to use for decompression.
     *
     * \details
     *
     * This value is currently only relevant for BGZF compressed streams/files. Note that these threads refer to the
     * total number of used threads, i.e. a value of 4 means that three extra threads are spawned. A value of 1 will
     * result in the regular GZ decompressor being used.
     *
     * The default value for this is 4 or "available CPUs" if that is less than 4. The reason is that for the
     * default compression levels of the GZip blocks, performance degrades with more than 4 threads. If you use files
     * compressed with stronger settings, more threads might be useful.
     *
     * **4 threads can provide up to 3.5x speed-up.**
     */
    size_t threads = std::max<size_t>(1, std::min<size_t>(4, std::thread::hardware_concurrency()));
};

/*!\brief A std::istream that automatically detects compressed streams and transparently decompresses them.
 * \ingroup stream
 * \details
 *
 * This is a c++ iostream compatible type that transparently decompresses compressed streams, i.e. you can open
 * uncompressed or compressed files with this stream and will always get uncompressed data. The format of compression
 * is detected automatically by default.
 *
 * See bio::io::compression_format for a list of currently supported formats.
 *
 * A filename may be provided, in which case this stream behaves like a file stream. Or an existing input stream can be
 * given which is then wrapped by this stream.
 *
 * ### Example
 *
 * \snippet test/snippet/stream/transparent_istream.cpp construction
 *
 * Do not use extra threads for decompression:
 *
 * \snippet test/snippet/stream/transparent_istream.cpp decompression
 */
class transparent_istream : public std::basic_istream<char>
{
private:
    /*TODO
     * evaluate whether to use custom sized buffer on both streams
     * evaluate whether to use custom sized buffer strings passed in (what happens to stuff in old buffer?)
     * tie in threads
     */

    //!\brief The options.
    transparent_istream_options options_;
    //!\brief The stream buffer.
    std::vector<char>           stream1_buffer;
    //!\brief The stream buffer.
    std::vector<char>           stream2_buffer;
    //!\brief Filename (if stream was opened from path).
    std::filesystem::path       filename_;
    //!\brief Filename after possible compression extensions have been removed.
    std::filesystem::path       truncated_filename_;
    //!\brief The decompressor selected.
    compression_format          selected_compression{};

    //!\brief The type of the internal stream pointers. Allows dynamically setting ownership management.
    using stream_ptr_t = std::unique_ptr<std::basic_istream<char>, std::function<void(std::basic_istream<char> *)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_istream<char> *) {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_istream<char> * ptr) { delete ptr; }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};

    //!\brief This function reads the magic bytes from the stream and adds a decompression layer if necessary.
    void set_secondary_stream()
    {
        assert(primary_stream->good());

        /* detect compression format */
        std::string magic_header = detail::read_magic_header(*primary_stream);

        if (options_.compression == compression_format::detect)
        {
            selected_compression = detail::detect_format_from_magic_header(magic_header);
        }
        else
        {
            selected_compression = options_.compression;
            if (!detail::header_matches_dyn(selected_compression, magic_header))
                throw file_open_error{"The file has a different compression format than the one selected."};
        }

        // Thread handling
        if (selected_compression == compression_format::bgzf && options_.threads == 1)
            selected_compression = compression_format::gz;

        std::span<std::string> file_extensions{};
        std::istream *         sec = nullptr;
        switch (selected_compression)
        {
            case compression_format::bgzf:
                // "- 1" because bgzf spawns **additional** threads, but user sets total
                sec             = detail::make_istream<compression_format::bgzf>(*primary_stream, options_.threads - 1);
                file_extensions = compression_traits<compression_format::bgzf>::file_extensions;
                break;
            case compression_format::gz:
                sec             = detail::make_istream<compression_format::gz>(*primary_stream);
                file_extensions = compression_traits<compression_format::gz>::file_extensions;
                break;
            case compression_format::bz2:
                sec             = detail::make_istream<compression_format::bz2>(*primary_stream);
                file_extensions = compression_traits<compression_format::bz2>::file_extensions;
                break;
            case compression_format::zstd:
                sec             = detail::make_istream<compression_format::zstd>(*primary_stream);
                file_extensions = compression_traits<compression_format::zstd>::file_extensions;
                break;
            default:
                break;
        }

        if (sec == nullptr) // secondary stream will be primary stream
        {
            secondary_stream = stream_ptr_t{&*primary_stream, stream_deleter_noop};
        }
        else
        {
            secondary_stream = stream_ptr_t{sec, stream_deleter_default};

            // truncate the filename in truncated_filename_ to show that decompression has taken place
            if (filename_.has_extension())
            {
                std::string extension = filename_.extension().string().substr(1);
                if (std::ranges::find(file_extensions, extension) != std::ranges::end(file_extensions))
                    truncated_filename_.replace_extension();
            }
        }
    }

    //!\brief Initialise state of object.
    void init()
    {
        truncated_filename_ = filename_;

        // possibly add intermediate compression stream
        set_secondary_stream();

        // make this behave like the secondary stream
        this->rdbuf(secondary_stream->rdbuf());
    }

    //!\brief Reset secondary stream after seek on primary.
    void post_seek(compression_format const old_compression)
    {
        if (primary_stream->fail())
            throw bio_error{"Seek failed on input stream."};

        secondary_stream.reset();

        options_.compression = compression_format::detect;
        init();

        if (old_compression != selected_compression)
        {
            throw bio_error{
              "Cannot restart decompression after seek on compressed file.\nThe file is not compressed "
              "blockwise and/or the seek position is not the beginning of a block.\nOld: ",
              (size_t)old_compression,
              " new: ",
              (size_t)selected_compression};
        }
    }

public:
    /*!\name Associated types
     * \{
     */
    using char_type   = char;                   //!< The character type.
    using traits_type = std::char_traits<char>; //!< Character traits.
    using pos_type    = traits_type::pos_type;  //!< Position type (typically unsigned).
    using off_type    = traits_type::off_type;  //!< Offset type (typically signed).
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    transparent_istream() : std::basic_istream<char>{} {}                   //!< Call default constructor of base.
    transparent_istream(transparent_istream const &)             = delete;  //!< Deleted.
    transparent_istream & operator=(transparent_istream &&)      = default; //!< Defaulted.
    transparent_istream & operator=(transparent_istream const &) = delete;  //!< Defaulted.

    //!\brief Move construction swaps the members maually before setting the buffer s.t. it is not invalidated.
    transparent_istream(transparent_istream && rhs)
    {
        std::swap(options_, rhs.options_);
        std::swap(stream1_buffer, rhs.stream1_buffer);
        std::swap(stream2_buffer, rhs.stream2_buffer);
        std::swap(filename_, rhs.filename_);
        std::swap(truncated_filename_, rhs.truncated_filename_);
        std::swap(primary_stream, rhs.primary_stream);
        std::swap(secondary_stream, rhs.secondary_stream);

        this->set_rdbuf(secondary_stream->rdbuf());
    }

    /*!\brief Construct from a filename.
     * \param[in] filename  The filename to open.
     * \param[in] options   See bio::io::transparent_istream_options.
     */
    explicit transparent_istream(std::filesystem::path       filename,
                                 transparent_istream_options options = transparent_istream_options{}) :
      options_{std::move(options)},
      filename_{std::move(filename)},
      primary_stream{new std::ifstream{}, stream_deleter_default}
    {
        stream1_buffer.resize(options.buffer1_size);

        primary_stream->rdbuf()->pubsetbuf(stream1_buffer.data(), stream1_buffer.size());
        static_cast<std::basic_ifstream<char> *>(primary_stream.get())
          ->open(filename_, std::ios_base::in | std::ios::binary);

        if (!primary_stream->good())
            throw file_open_error{"Could not open file ", filename_.string(), " for reading."};

        init();
    }

    /*!\brief Construct from a stream.
     * \param[in] stream    The stream to wrap.
     * \param[in] options   See bio::io::transparent_istream_options.
     */
    explicit transparent_istream(std::basic_istream<char> &  stream,
                                 transparent_istream_options options = transparent_istream_options{}) :
      options_{std::move(options)}, primary_stream{&stream, stream_deleter_noop}
    {
        init();
    }

    //!\overload
    template <meta::different_from<transparent_istream> temporary_stream_t>
        requires(movable_istream<temporary_stream_t> && !std::is_lvalue_reference_v<temporary_stream_t>)
    explicit transparent_istream(temporary_stream_t &&       stream,
                                 transparent_istream_options options = transparent_istream_options{}) :
      options_{std::move(options)}, primary_stream{new temporary_stream_t{std::move(stream)}, stream_deleter_default}
    {
        init();
    }
    //!\}

    //!\brief The filename this object was created from; empty if this object was not created from a file.
    std::filesystem::path const & filename() { return filename_; }

    /*!\brief The filename this object was created from without compression-specific suffix.
     *
     * \details
     *
     * If this object was created from e.g. "foo.fasta.gz", #filename() will return the full name, but this function
     * will return only "foo.fasta". This is useful for determining the format "inside" the compression.
     *
     * If this object was not created from a file, an empty path is returned.
     */
    std::filesystem::path const & truncated_filename() { return truncated_filename_; }

    /*!\brief Seek on the primary stream and reset secondary stream.
     * \details
     *
     * The function enables seeking to the beginning of another block in a file using block-compression.
     * This is used with indexed I/O on BGZF compressed files, but other forms of block-based
     * compression should also be supported (the decompression stream is re-created after the seek).
     */
    transparent_istream & seekg_primary(pos_type const pos)
    {
        compression_format old_compression = selected_compression;

        primary_stream->seekg(pos);

        post_seek(old_compression);

        return *this;
    }

    //!\overload
    transparent_istream & seekg_primary(off_type const off, std::ios_base::seekdir dir)
    {
        compression_format old_compression = selected_compression;

        primary_stream->seekg(off, dir);

        post_seek(old_compression);

        return *this;
    }
};

} // namespace bio::io
