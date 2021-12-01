// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::transparent_ostream.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <span>
#include <thread>

#include <bio/exception.hpp>
#include <bio/stream/compression.hpp>
#include <bio/stream/concept.hpp>
#include <bio/stream/detail/make_stream.hpp>

namespace bio
{

//!\brief Options that can be provided to bio::transparent_ostream.
//!\ingroup stream
struct transparent_ostream_options
{
    //!\brief Size of the buffer used when opening a file from a filename.
    size_t buffer1_size = 1024 * 1024;
    //!\brief Size of the buffer used for the compression stream.
    size_t buffer2_size = 1024 * 1024 * 4;

    /*!\brief Which compressor to use.
     *
     * \details
     *
     * For ostream opened from filenames, the default is to detect the desired compression from the file's extension.
     *
     * For transparent ostreams constructed from existing stream objects compression_format::detect has no
     * effect as there is nothing to detect (it instead defaults to compression_format::none).
     *
     * Set this manually to force a desired compression format (this is especially useful when working on streams).
     */
    compression_format compression = compression_format::detect;

    /*!\brief The compression level to use by the algorithm.
     *
     * \details
     *
     * The default value is -1 which maps to the default value of the respective algorithm (6 for GZ/BGZF and 9 for
     * BZip2). ZLIB macros and numeric values between -1 and 9 are supported.
     */
    int compression_level = -1;

    /*!\brief Maximum number of threads to use for compression.
     *
     * \details
     *
     * This value is currently only relevant for BGZF compressed streams/files. Note that these threads refer to the
     * total number of used threads, i.e. a value of 4 means that three extra threads are spawned.
     *
     * The default value for this is 8 or "available CPUs" if that is less than 8. The reason is that for the
     * default compression levels of the GZip blocks, there is only little speed-up after 8 threads.
     *
     * TODO double-check max-value
     *
     * ### Attention
     *
     * A value of 1 is currently not supported by the BGZF implementation!
     */
    size_t threads = std::max<size_t>(1, std::min<size_t>(8, std::thread::hardware_concurrency()));
};

/*!\brief A std::ostream that automatically detects compressed streams and transparently decompresses them.
 * \ingroup stream
 * \details
 *
 * This is a c++ iostream compatible type that transparently compresses streams, i.e. depending on run-time arguments
 * provided to this stream object, the data written to it will be compressed (or not).
 *
 * See bio::compression_format for a list of currently supported formats.
 *
 * A filename may be provided, in which case this stream behaves like a file stream and the format of compression is
 * detected from the filename (extension). Or an existing input stream can be given which is then wrapped by this
 * stream; in this case the type of compression needs to be specified via the options.
 *
 * ### Example
 *
 * ```cpp
 * std::string_view content = "FOOBAR";
 *
 * bio::transparent_ostream s1{"my_file.txt"};      // behaves like std::ofstream
 * s1 << content;
 *
 * bio::transparent_ostream s2{"my_file.txt.gz"};   // data is written GZ compressed
 * s2 << content;
 *
 * bio::transparent_ostream s3{std::cout};          // wrap standard input; no compression
 * s3 << content;
 *
 * bio::transparent_ostream s4{std::cout, { .compression = bio::compression_format::bgzf } }; // use BGZF compression
 * s4 << content;
 * ```
 *
 * Explicitly request BZGF compression and a total of two threads (one extra thread for compression):
 *
 * ```cpp
 * std::string_view content = "FOOBAR";
 *
 * bio::transparent_ostream s{std::cout,
 *                            { .compression = bio::compression_format::bgzf, .threads = 2 } };
 * s << content;
 * ```
 */
class transparent_ostream : public std::basic_ostream<char>
{
private:
    /*TODO
     * evaluate whether to use custom sized buffer on both streams
     * evaluate whether to use custom sized buffer strings passed in (what happens to stuff in old buffer?)
     */

    //!\brief The options.
    transparent_ostream_options options_;
    //!\brief The stream buffer.
    std::vector<char>           stream1_buffer;
    //!\brief The stream buffer.
    std::vector<char>           stream2_buffer;
    //!\brief Filename (if stream was opened from path).
    std::filesystem::path       filename_;
    //!\brief Filename after possible compression extensions have been removed.
    std::filesystem::path       truncated_filename_;

    //!\brief The type of the internal stream pointers. Allows dynamically setting ownership management.
    using stream_ptr_t = std::unique_ptr<std::basic_ostream<char>, std::function<void(std::basic_ostream<char> *)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_ostream<char> *) {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_ostream<char> * ptr) { delete ptr; }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};

    //!\brief This function reads the magic bytes from the stream and adds a decompression layer if necessary.
    void set_secondary_stream()
    {
        assert(primary_stream->good());

        /* detect compression format */
        if (options_.compression == compression_format::detect)
        {
            if (filename_.empty()) // constructed from ostream → there is nothing to detect
                options_.compression = compression_format::none;
            else
                options_.compression = detail::detect_format_from_filename(filename_);
        }

        // Thread handling
        if (options_.compression == compression_format::bgzf)
        {
            // TODO catch threads == 0 as error
            if (options_.threads == 1) // TODO this needs a real resolution
                throw file_open_error{"BGZF compression with only one thread is currently not supported."};
            else
                --options_.threads; // bgzf spawns **additional** threads, but user sets total
        }

        std::span<std::string> file_extensions{};
        std::ostream *         sec = nullptr;
        switch (options_.compression)
        {
            case compression_format::bgzf:
                sec             = detail::make_ostream<compression_format::bgzf>(*primary_stream,
                                                                     options_.threads,
                                                                     static_cast<size_t>(8ul),
                                                                     options_.compression_level);
                file_extensions = compression_traits<compression_format::bgzf>::file_extensions;
                break;
            case compression_format::gz:
                sec = detail::make_ostream<compression_format::gz>(*primary_stream, options_.compression_level);
                file_extensions = compression_traits<compression_format::gz>::file_extensions;
                break;
            case compression_format::bz2:
                sec = detail::make_ostream<compression_format::bz2>(*primary_stream, options_.compression_level);
                file_extensions = compression_traits<compression_format::bz2>::file_extensions;
                break;
            case compression_format::zstd:
                sec = detail::make_ostream<compression_format::zstd>(*primary_stream, options_.compression_level);
                file_extensions = compression_traits<compression_format::zstd>::file_extensions;
                break;
            default:
                break;
        }

        if (sec == nullptr)
            secondary_stream = stream_ptr_t{&*primary_stream, stream_deleter_noop};
        else
            secondary_stream = stream_ptr_t{sec, stream_deleter_default};

        // truncate the filename in truncated_filename_ to show that compression has taken place
        if (filename_.has_extension())
        {
            std::string extension = filename_.extension().string().substr(1);
            if (std::ranges::find(file_extensions, extension) != std::ranges::end(file_extensions))
                truncated_filename_.replace_extension();
        }
    }

    //!\brief Initialise state of object.
    void init()
    {
        truncated_filename_ = filename_;

        // possibly add intermediate compression stream
        set_secondary_stream();
        assert(secondary_stream != nullptr);

        // make this behave like the secondary stream
        this->rdbuf(secondary_stream->rdbuf());
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Manually defined default constructor that behaves as expected.
    transparent_ostream() : std::basic_ostream<char> {}
    {}                                                                     //!< Call default constructor of base.
    transparent_ostream(transparent_ostream const &) = delete;             //!< Defaulted.
    transparent_ostream & operator=(transparent_ostream const &) = delete; //!< Defaulted.
    // TODO double check that this works:
    transparent_ostream & operator=(transparent_ostream &&) = default; //!< Defaulted.

    //!\brief Manually defined move constructor that behaves as expected.
    transparent_ostream(transparent_ostream && rhs)
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
     * \param[in] options   See bio::transparent_ostream_options.
     *
     * \details
     *
     * The compression format is auto-detected from the filename by default. It can manually be selected via the
     * options.
     *
     * The stream is opened in binary mode and provided with a buffer the size of options.buffer1_size.
     */
    explicit transparent_ostream(std::filesystem::path       filename,
                                 transparent_ostream_options options = transparent_ostream_options{}) :
      options_{std::move(options)},
      filename_{std::move(filename)},
      primary_stream{new std::ofstream{}, stream_deleter_default}
    {
        stream1_buffer.resize(options.buffer1_size);

        primary_stream->rdbuf()->pubsetbuf(stream1_buffer.data(), stream1_buffer.size());
        static_cast<std::basic_ofstream<char> *>(primary_stream.get())
          ->open(filename_, std::ios_base::out | std::ios::binary);

        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename_.string() + " for reading."};

        init();
    }

    /*!\brief Construct from a stream.
     * \param[in] stream  The stream to wrap.
     * \param[in] options See bio::transparent_ostream_options.
     *
     * \details
     *
     * The compression format is "none" by default. It can manually be selected via the options.
     */
    explicit transparent_ostream(std::ostream &              stream,
                                 transparent_ostream_options options = transparent_ostream_options{}) :
      options_{std::move(options)}, primary_stream{&stream, stream_deleter_noop}
    {
        init();
    }

    //!\overload
    template <typename temporary_stream_t>
        requires(!std::same_as<temporary_stream_t, transparent_ostream> && movable_ostream<temporary_stream_t> &&
                 !std::is_lvalue_reference_v<temporary_stream_t>)
    explicit transparent_ostream(temporary_stream_t &&       stream,
                                 transparent_ostream_options options = transparent_ostream_options{}) :
      options_{std::move(options)}, primary_stream{new temporary_stream_t{std::move(stream)}, stream_deleter_default}
    {
        init();
    }

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
};

} // namespace bio
