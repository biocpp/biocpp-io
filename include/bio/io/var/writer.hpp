// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::var::writer.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <filesystem>
#include <iosfwd>

#include <bio/io/detail/writer_base.hpp>
#include <bio/io/format/bcf_output_handler.hpp>
#include <bio/io/format/vcf_output_handler.hpp>
#include <bio/io/var/header.hpp>
#include <bio/io/var/writer_options.hpp>

namespace bio::io::var
{

// ----------------------------------------------------------------------------
// writer
// ----------------------------------------------------------------------------

/*!\brief A class for writing variant files, e.g. VCF, BCF, GVCF.
 * \tparam option_args_t Arguments that are forwarded to bio::io::var::writer_options.
 * \ingroup var
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
 *   1. bio::io::detail::field::chrom
 *   2. bio::io::detail::field::pos
 *   3. bio::io::detail::field::id
 *   4. bio::io::detail::field::ref
 *   5. bio::io::detail::field::alt
 *   6. bio::io::detail::field::qual
 *   7. bio::io::detail::field::filter
 *   8. bio::io::detail::field::info
 *   9. bio::io::detail::field::genotypes
 *
 * These fields correspond to the order and names defined in the VCF specification. The value conventions
 * also correspond to the VCF specification (i.e. 1-based positions) although many fields can
 * also be provided in BCF like types (e.g. IDX values instead of strings).
 * See below for the list of supported types and the semantic implications.
 *
 * This writer supports the following formats:
 *
 *   1. VCF (see also bio::io::vcf)
 *   2. BCF (see also bio::io::bcf)
 *
 * If you only need to write VCF and not BCF and you have all your column data as strings,
 * you can use bio::io::txt::writer instead of this writer (it will be easier to use and faster).
 *
 * ### Creating a writer
 *
 * This creates a writer with a valid header:
 *
 * \snippet test/snippet/var/var_writer.cpp creation
 *
 * **If you copy'n'paste from this example, make sure that columns in the last line are tab-separated and not
 * space-separated!**
 *
 * This example is used as "prefix" for some of the following examples.
 *
 * ### Writing a record
 *
 * Create a record using the bio::io::var::default_record and setting the members:
 *
 * \snippet test/snippet/var/var_writer.cpp simple_usage_file
 *
 * ### Writing a record without having a record
 *
 * You can use #emplace_back() to write the fields directly without creating a record first:
 *
 * \snippet test/snippet/var/var_writer.cpp emplace_back
 *
 * This is especially helpful if your fields exist in other separate data structures already.
 * Be aware that it is easier to mess up the order of the arguments this way.
 * The order/composition is exactly as in bio::io::var::record!
 *
 * ### Specifying options
 *
 * This snippet demonstrates how to create VCF files with windows line-endings:
 *
 * \snippet test/snippet/var/var_writer.cpp options
 *
 * For more advanced options, see bio::io::var::writer_options.
 *
 * ### Combining reading and writing
 *
 * This simple snippet demonstrates how to pipe from a reader into a writer (transparently converting):
 *
 * \snippet test/snippet/var/var_writer.cpp inout
 *
 * Views can be used to modify or filter the records before being written:
 *
 * \snippet test/snippet/var/var_writer.cpp inout2
 *
 * A more "traditional" programming style with the same semantic would look like this:
 *
 * \snippet test/snippet/var/var_writer.cpp inout3
 *
 * ### Field type requirements
 *
 * TODO
 *
 */
template <typename... option_args_t>
class writer : public writer_base<writer<option_args_t...>, writer_options<option_args_t...>>
{
private:
    //!\brief The base class.
    using base_t      = writer_base<writer<option_args_t...>, writer_options<option_args_t...>>;
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

    using base_t::write_record;

public:
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    writer()                           = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    writer(writer const &)             = delete;
    //!\brief Move construction is defaulted.
    writer(writer &&)                  = default;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    writer & operator=(writer const &) = delete;
    //!\brief Move assignment is defaulted.
    writer & operator=(writer &&)      = default;
    //!\brief Destructor which can potentially throw.

    using base_t::operator=;

    // clang-format off
    //!\copydoc bio::io::writer_base::writer_base(std::filesystem::path const & filename, format_type const & fmt, options_t const & opt = options_t{})
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
    //!\copydoc bio::io::writer_base::writer_base(std::ostream & str, format_type const & fmt, options_t const & opt = options_t{})
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

    /*!\brief Write a record to the file by passing individual fields.
     * \param[in] chrom     The CHROM parameter.
     * \param[in] pos       The POS parameter.
     * \param[in] id        The ID parameter.
     * \param[in] ref       The REF parameter.
     * \param[in] alt       The ALT parameter.
     * \param[in] qual      The QUAL parameter.
     * \param[in] filter    The FILTER parameter.
     * \param[in] info      The INFO parameter.
     * \param[in] genotypes The GENOTYPES parameter.
     *
     * \details
     *
     * This function accepts exactly 9 arguments and assumes that these
     * correspond exactly to the members of the bio::io::var::record class.
     *
     * ### Example
     *
     * \snippet test/snippet/var/var_writer.cpp emplace_back
     * TODO add snippet!
     */
    void emplace_back(auto && chrom,
                      auto && pos,
                      auto && id,
                      auto && ref,
                      auto && alt,
                      auto && qual,
                      auto && filter,
                      auto && info,
                      auto && genotypes)
    {
        auto record = tie_record(chrom, pos, id, ref, alt, qual, filter, info, genotypes);
        // TODO check concept
        write_record(record);
    }

    /*!\brief Write a record to the file.
     * \tparam field_types Types of the fields in the record.
     * \tparam field_ids   IDs of the fields in the record.
     * \param[in] r        The record to write.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     */
    template <typename... member_types>
    void push_back(record<member_types...> const & r)
    {
        write_record(r);
        // TODO check concept
    }

    //!\overload
    template <typename... member_types>
    void push_back(record<member_types...> & r)
    {
        write_record(r); // pass as non-const to allow parsing views that are not const-iterable
        // TODO check concept
    }

    //!\overload
    template <typename... member_types>
    void push_back(record<member_types...> && r)
    {
        write_record(r); // pass as non-const to allow parsing views that are not const-iterable
        // TODO check concept
    }

    //!\brief Get the header used by the format.
    bio::io::var::header const & header()
    {
        return std::visit(meta::overloaded{[](std::monostate) {},
                                           [](auto const & handler) { return handler.get_header(); }},
                          format_handler);
    }

    //!\brief Set the header to the given value.
    template <typename header_t>
        requires std::same_as<bio::io::var::header, std::remove_cvref_t<header_t>>
    void set_header(header_t && hdr)
    {
        if (!init_state)
            throw bio_error{"You cannot change the header after I/O has happened."};

        std::visit(meta::overloaded{[](std::monostate) {},
                                    [&hdr](auto & handler) { handler.set_header(std::forward<header_t>(hdr)); }},
                   format_handler);
    }
};

} // namespace bio::io::var
