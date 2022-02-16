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
 * \tparam option_args_t Arguments that are forwarded to bio::var_io::writer_options.
 * \ingroup var_io
 *
 * \details
 *
 * ### Introduction
 *
 * Variant files are files that contain sequence variation information. Well-known formats include
 * VCF and BCF.
 *
 * The Variant I/O writer supports writing the following fields:
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
 * These fields correspond to the order and names defined in the VCF specification. The value conventions
 * also correspond to the VCF specification (i.e. 1-based positions) although many fields can
 * also be provided in BCF like types (e.g. IDX values instead of strings).
 * See below for the list of supported types and the semantic implications.
 *
 * This writer supports the following formats:
 *
 *   1. VCF (see also bio::vcf)
 *   2. BCF (see also bio::bcf) [TODO NOT YET]
 *
 * If you only need to write VCF and not BCF and you have all your column data as strings,
 * you can use bio::plain_io::writer instead of this writer (it will be easier to use and faster).
 *
 * ### Creating a writer
 *
 * This creates a writer with a valid header:
 *
 * \snippet test/snippet/var_io/var_io_writer.cpp creation
 *
 * **If you copy'n'paste from this example, make sure that columns in the last line are tab-separated and not
 * space-separated!**
 *
 * This example is used as "prefix" for some of the following examples.
 *
 * ### Writing a record
 *
 * Create a record using the bio::var_io::default_record and setting the members:
 *
 * \snippet test/snippet/var_io/var_io_writer.cpp simple_usage_file
 *
 * ### Writing a record without having a record
 *
 * You can use #emplace_back() to write the fields directly without creating a record first:
 *
 * \snippet test/snippet/var_io/var_io_writer.cpp emplace_back
 *
 * This is especially helpful if your fields exist in other separate data structures already.
 * Be aware that it is easier to mess up the order of the arguments this way.
 * The order/composition can specified by bio::vtag as first argument (see next example).
 * If it is omitted, it is equal to bio::var_io::default_field_ids.
 *
 * The #emplace_back() function can be used to write fewer fields:
 *
 * \snippet test/snippet/var_io/var_io_writer.cpp emplace_back2
 *
 * These three fields are required; the missing fields are replaced with ".".
 *
 * ### Specifying options
 *
 * This snippet demonstrates how to create VCF files with windows line-endings:
 *
 * \snippet test/snippet/var_io/var_io_writer.cpp options
 *
 * For more advanced options, see bio::var_io::writer_options.
 *
 * ### Combining reading and writing
 *
 * This simple snippet demonstrates how to pipe from a reader into a writer (transparently converting):
 *
 * \snippet test/snippet/var_io/var_io_writer.cpp inout
 *
 * Views can be used to modify or filter the records before being written:
 *
 * \snippet test/snippet/var_io/var_io_writer.cpp inout2
 *
 * A more "traditional" programming style with the same semantic would look like this:
 *
 * \snippet test/snippet/var_io/var_io_writer.cpp inout3
 *
 * ### Field type requirements
 *
 * TODO
 *
 */
template <typename... option_args_t>
class writer : public writer_base<writer_options<option_args_t...>>
{
private:
    //!\brief The base class.
    using base_t      = writer_base<writer_options<option_args_t...>>;
    //!\brief Inherit the format_type definition.
    using format_type = typename base_t::format_type;
    /* Implementation note
     * format_type is "inherited" as private here to avoid appearing twice in the documentation.
     * Its actual visibility is public because it is public in the base class.
     */

    //!\brief Make the format handler visible.
    using base_t::format_handler;
    //!\brief Make the init_state handler visible.
    using base_t::init_state;

public:
    // TODO wrap this, so we don't return reference to base
    using base_t::operator=;

    // clang-format off
    //!\copydoc bio::writer_base::writer_base(std::filesystem::path const & filename, format_type const & fmt, options_t const & opt = options_t{})
    // clang-format on
    writer(std::filesystem::path const &            filename,
           format_type const &                      fmt,
           writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{filename, fmt, opt}
    {}

    //!\overload
    explicit writer(std::filesystem::path const &            filename,
                    writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{filename, opt}
    {}

    // clang-format off
    //!\copydoc bio::writer_base::writer_base(std::ostream & str, format_type const & fmt, options_t const & opt = options_t{})
    // clang-format on
    writer(std::ostream &                           str,
           format_type const &                      fmt,
           writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{str, fmt, opt}
    {}

    //!\overload
    template <movable_ostream temporary_stream_t>
        //!\cond REQ
        requires(!std::is_lvalue_reference_v<temporary_stream_t>)
    //!\endcond
    writer(temporary_stream_t &&                    str,
           format_type const &                      fmt,
           writer_options<option_args_t...> const & opt = writer_options<option_args_t...>{}) :
      base_t{std::move(str), fmt, opt}
    {}

    //!\brief Destructor which can potentially throw.
    ~writer() noexcept(false) = default;

    // prevent the overload below from removing the overload from base_t
    using base_t::emplace_back;

    /*!\brief Write a record to the file by passing individual fields.
     * \param[in] args   The fields to be written.
     *
     * \details
     *
     * This function is the same as bio::writer_base::emplace_back, except that the field_ids can be
     * omitted. If the number of arguments 10, bio::var_io::default_field_ids is chosen; if it is
     */
    void emplace_back(auto &&... args)
    {
        static_assert(sizeof...(args) == default_field_ids.size || sizeof...(args) == default_field_ids.size - 1,
                      "emplace_back() has to be called with 9 or 10 arguments (including private_data).");

        if constexpr (sizeof...(args) == default_field_ids.size - 1)
        {
            // TODO replace this with some metaprogramming?
            base_t::emplace_back(vtag<field::chrom,
                                      field::pos,
                                      field::id,
                                      field::ref,
                                      field::alt,
                                      field::qual,
                                      field::filter,
                                      field::info,
                                      field::genotypes>,
                                 args...);
        }
        else
        {
            base_t::emplace_back(default_field_ids, args...);
        }
    }

    //!\brief Get the header used by the format.
    bio::var_io::header const & header()
    {
        return std::visit(
          detail::overloaded{[](std::monostate) {}, [](auto const & handler) { return handler.get_header(); }},
          format_handler);
    }

    //!\brief Set the header to the given value.
    template <typename header_t>
        requires std::same_as<bio::var_io::header, std::remove_cvref_t<header_t>>
    void set_header(header_t && hdr)
    {
        if (!init_state)
            throw std::logic_error{"You cannot change the header after I/O has happened."};

        std::visit(detail::overloaded{[](std::monostate) {},
                                      [&hdr](auto & handler) { handler.set_header(std::forward<header_t>(hdr)); }},
                   format_handler);
    }
};

} // namespace bio::var_io
