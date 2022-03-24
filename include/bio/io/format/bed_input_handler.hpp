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

#include <bio/plain_io/reader.hpp>

namespace bio
{
template <>
class format_input_handler<bed> : public format_input_handler_base<format_input_handler<bed>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The type of the CRTP base class.
    using base_t = format_input_handler_base<format_input_handler<bed>>;
    using base_t::parse_field;
    using base_t::parse_field_aux;
    using base_t::stream;

    //!\brief Befriend the base class to enable CRTP.
    friend base_t;
    //!\}

    //!\brief Print an error message with current line number in diagnostic.
    [[noreturn]] void error(auto const &... messages) const
    {
        std::string message = "[SeqAn3 BED format error in line " + detail::to_string(line) + "] ";
        ((message += detail::to_string(messages)), ...);

        throw format_error{message};
    }

    /*!\name Options
     * \{
     */
    //!\brief Whether to print warnings or not.
    bool print_warnings = true;
    //!\}

    /*!\name Raw record handling
     * \{
     */
    //!\brief The fields that this format supports [the base class accesses this type].
    using format_fields = decltype(ann_io::default_field_ids);
    //!\brief Type of the raw record.
    using raw_record_type =
      record<format_fields, seqan3::list_traits::repeat<format_fields::size, std::string_view>>;

    //!\brief Type of the low-level iterator.
    using lowlevel_iterator = detail::plaintext_input_iterator<plain_io::record_kind::line_and_fields>;

    //!\brief The raw record.
    raw_record_type   raw_record;
    //!\brief The header.
    // var_io::header    header;
    //!\brief Lowlevel stream iterator.
    lowlevel_iterator file_it;
    //!\brief Cache of the chromosome string.
    std::string       last_chrom;
    //!\brief Current line number in file.
    size_t            line           = 0;

    //!\brief Read the raw record [the base class invokes this function].
    void read_raw_record()
    {
        ++line;
        ++file_it;

        if (size_t field_num = file_it->fields.size(); field_num < 3)
            error("Expected at least 3 fields but got ", field_num);

        get<field::chrom>(raw_record)      = (*file_it).fields[0];
        get<field::chromStart>(raw_record) = (*file_it).fields[1];
        get<field::chromEnd>(raw_record)   = (*file_it).fields[2];
    }
    //!\}

    /*!\name Parsed record handling
     * \brief This is mostly done via the defaults in the base class.
     * \{
     */

    // implementation after class
    // template <typename t>
    //     requires detail::is_info_element_value_type<t> || detail::is_genotype_element_value_type<t>
    // static void init_element_value_type(ann_io::value_type_id const id, t & output);

    // implementation after class
    // struct parse_element_value_type_fn;

    // implementation after class
    // static size_t parse_element_value_type(var_io::value_type_id const               id,
    //                                        std::string_view const                    input_string,
    //                                        detail::is_info_element_value_type auto & output);

    //!\brief Parse the CHROM field. Note there is no index, as BED files don't store them.
    // void parse_field(vtag_t<field::chrom> const & /**/, auto & parsed_field)
    // {
    //     using parsed_field_t       = std::remove_cvref_t<decltype(parsed_field)>;
    //     std::string_view raw_field = get<field::chrom>(raw_record);
    //
    //     if (raw_field != last_chrom) last_chrom = raw_field;
    //     parsed_field = static_cast<parsed_field_t>(raw_field);
    // }

    /* chromStart and chromEnd are handled correctly by default */
    //!\}
public:
    /*!\name Constructors, destructor and assignment.
     * \{
     */
    format_input_handler()                                         = default; //!< Defaulted.
    format_input_handler(format_input_handler const &)             = delete;  //!< Deleted.
    format_input_handler(format_input_handler &&)                  = default; //!< Defaulted.
    ~format_input_handler()                                        = default; //!< Defaulted.
    format_input_handler & operator=(format_input_handler const &) = delete;  //!< Deleted.
    format_input_handler & operator=(format_input_handler &&)      = default; //!< Defaulted.

    /*!\brief Construct with an options object.
     * \param[in,out] str The input stream.
     * \param[in] options An object with options for the input handler.
     * \details
     *
     * The options argument is typically bio::var_io::reader_options, but any object with a subset of similarly named
     * members is also accepted. See bio::format_input_handler<vcf> for the supported options and defaults.
     */
    template <typename options_t>
    format_input_handler(std::istream & str, options_t const & options) : base_t{str}, file_it{str, false /*no_init!*/}

    {
        // extract options
        if constexpr (requires { (bool)options.print_warnings; })
        {
            print_warnings = options.print_warnings;
        }

        // std::string header_string;
        // while (file_it != std::default_sentinel && file_it.peak() == '#')
        // {
        //     ++file_it;
        //     ++line;
        //     header_string += file_it->line;
        //     header_string += "\n";
        // }
        // header = var_io::header{std::move(header_string)};
    }

    //!\brief Construct with only an input stream.
    format_input_handler(std::istream & str) : format_input_handler{str, int{}} {}
    //!\}

    //!\brief Return a reference to the header contained in the input handler.
    // var_io::header const & get_header() const { return header; }
};
}
