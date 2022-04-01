// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::format_input_handler<bio::sam>.
 * \author Svenja Mehringer <svenja.mehringer AT decode.is>
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/type_list/traits.hpp>

#include <bio/format/format_input_handler.hpp>
#include <bio/format/sam.hpp>
#include <bio/map_io/header.hpp>
#include <bio/map_io/misc.hpp>
#include <bio/map_io/sam_flag.hpp>
#include <bio/map_io/sam_tag_dictionary.hpp>
#include <bio/plain_io/reader.hpp>

namespace bio
{

/*!\brief Format input handler for the SAM format (bio::sam).
 * \ingroup format
 * \details
 *
 * ### Attention
 *
 * Most users should not perform I/O through input/output handlers but should instead use the respective
 * readers/writers. See the overview (TODO link) for more information.
 *
 * ### Options
 *
 * The following options are considered if the respective member variable is availabele in the object passed to
 * the constructor:
 *
 * | Member          | Type    | Default | Description                                                      |
 * |-----------------|---------|---------|------------------------------------------------------------------|
 * |`print_warnings` |`bool`   | `false` | Whether to print non-critical warnings to seqan3::debug_stream   |
 *
 */
template <>
class format_input_handler<sam> : public format_input_handler_base<format_input_handler<sam>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The type of the CRTP base class.
    using base_t = format_input_handler_base<format_input_handler<sam>>;
    using base_t::parse_field;
    using base_t::parse_field_aux;
    using base_t::stream;

    //!\brief Befriend the base class to enable CRTP.
    friend base_t;
    //!\}

    /*!\name Options
     * \{
     */
    //!\brief Whether to print warnings or not.
    bool print_warnings = true;
    //!\}

    //!\brief Throw a bio::format_error with an error message with current line number in diagnostic.
    [[noreturn]] void error(auto const &... messages) const
    {
        std::string message = "[B.I.O sam format error in line " + detail::to_string(line) + "] ";
        ((message += detail::to_string(messages)), ...);

        throw format_error{message};
    }

    //!\brief Print a B.I.O warning message with current line number in diagnostic.
    void warning(auto const &... messages) const
    {
        if (print_warnings)
        {
            seqan3::debug_stream << "[B.I.O sam format warning in line " << line << "] ";
            (seqan3::debug_stream << ... << messages);
            seqan3::debug_stream << std::endl;
        }
    }

    /*!\name Raw record handling
     * \{
     */
    //!\brief The fields that this format supports [the base class accesses this type].
    using format_fields = decltype(map_io::default_field_ids);
    //!\brief Type of the raw record.
    using raw_record_type =
      bio::record<format_fields, seqan3::list_traits::repeat<format_fields::size, std::string_view>>;

    //!\brief Type of the low-level iterator.
    using lowlevel_iterator = detail::plaintext_input_iterator<plain_io::record_kind::line_and_fields>;

    //!\brief The raw record.
    raw_record_type   raw_record;
    //!\brief The header.
    map_io::header    header;
    //!\brief Lowlevel stream iterator.
    lowlevel_iterator file_it;
    //!\brief Cache of the chromosome string.
    std::string       last_rname;
    //!\brief Cache of the chromosome index.
    int32_t           last_rname_idx = -1;
    //!\brief Current line number in file.
    size_t            line           = 0;

    //!\brief Read the raw record [the base class invokes this function].
    void read_raw_record()
    {
        ++line;
        ++file_it;

        // TODO assert number of fields

        get<field::qname>(raw_record) = (*file_it).fields[0];
        get<field::flag>(raw_record)  = (*file_it).fields[1];
        get<field::rname>(raw_record) = (*file_it).fields[2];
        get<field::pos>(raw_record)   = (*file_it).fields[3];
        get<field::mapq>(raw_record)  = (*file_it).fields[4];
        get<field::cigar>(raw_record) = (*file_it).fields[5];
        get<field::rnext>(raw_record) = (*file_it).fields[6];
        get<field::pnext>(raw_record) = (*file_it).fields[7];
        get<field::tlen>(raw_record)  = (*file_it).fields[8];
        get<field::seq>(raw_record)   = (*file_it).fields[9];
        get<field::qual>(raw_record)  = (*file_it).fields[10];

        // fields[10].end() that is guaranteed to be char*
        char const * end_qual        = (*file_it).fields[10].data() + (*file_it).fields[10].size() + 1 /*\t or \n*/;
        // line.end() that is guaranteed to be char*
        char const * end_line        = (*file_it).line.data() + (*file_it).line.size();
        // SAM tags go from end of qual til end of line (possibly empty)
        get<field::tags>(raw_record) = std::string_view{end_qual, static_cast<size_t>(end_line - end_qual)};
    }
    //!\}

    /*!\name Parsed record handling
     * \brief This is mostly done via the defaults in the base class.
     * \{
     */
    //!\brief Overload for parsing QNAME.
    template <typename parsed_field_t>
    void parse_field(vtag_t<field::qname> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field::qname>(raw_record);
        if (raw_field != "*")
            parse_field_aux(raw_field, parsed_field); // default parsing
    }

    //!\brief Overload for parsing FLAG.
    void parse_field(vtag_t<field::flag> const & tag, map_io::sam_flag & parsed_field)
    {
        uint16_t flag_integral{};
        parse_field(tag, flag_integral); // arithmetic parsing
        parsed_field = map_io::sam_flag{flag_integral};
    }

    //!\brief Parse the RNAME field. The reference names are stored in the header.
    template <typename parsed_field_t>
    void parse_field(vtag_t<field::rname> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field::rname>(raw_record);

        if (raw_field != "*")
        {
            size_t rname_pos = 0;

            if (auto it = header.rname_to_pos().find(raw_field); it == header.rname_to_pos().end())
            { // rname name was not in header, insert!
                header.push_back_rname(raw_field, raw_field.size(), "");

                rname_pos = header.rnames().size() - 1;

                warning("The reference sequence name \"", raw_field, "\" is not present in the header.");
            }
            else
            {
                rname_pos = it->second;
            }

            if constexpr (std::integral<std::remove_cvref_t<parsed_field_t>>)
                parsed_field = rname_pos;
            else
                parsed_field = header.rnames()[rname_pos];
        }
    }

    /* POS, MAPQ are handled correctly by default */

    //!\brief Overload for parsing CIGAR.
    void parse_field(vtag_t<field::cigar> const & /**/, std::vector<seqan3::cigar> & cigar_vector)
    {
        std::string_view raw_field = get<field::cigar>(raw_record);

        if (raw_field != "*")
        {
            uint32_t     cigar_count{};
            char const * ptr = raw_field.data();
            char const * end = ptr + raw_field.size();

            while (ptr < end)
            {
                auto const res = std::from_chars(ptr, end, cigar_count); // reads number up to next chatacter

                if (res.ec != std::errc{}) // parse number
                    error("Corrupted cigar string encountered");

                ptr = res.ptr + 1; // skip cigar operation character

                cigar_vector.emplace_back(cigar_count,
                                          seqan3::assign_char_strictly_to(*res.ptr, seqan3::cigar::operation{}));
            }
        }
    }

    //!\brief Parse the RNEXT field.
    template <typename parsed_field_t>
    void parse_field(vtag_t<field::rnext> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field::rnext>(raw_record);

        if (raw_field != "*")
        {
            if (raw_field == "=")
                raw_field = get<field::rname>(raw_record);

            size_t rname_pos = 0;

            if (auto it = header.rname_to_pos().find(raw_field); it == header.rname_to_pos().end())
            { // rname name was not in header, insert!
                header.push_back_rname(raw_field, raw_field.size(), "");

                rname_pos = header.rnames().size() - 1;

                warning("The mates reference sequence name \"", raw_field, "\" is not present in the header.");
            }
            else
            {
                rname_pos = it->second;
            }

            if constexpr (std::integral<std::remove_cvref_t<parsed_field_t>>)
                parsed_field = rname_pos;
            else
                parsed_field = header.rnames()[rname_pos];
        }
    }

    /* PNEXT, TLEN are handled correctly by default */

    //!\brief Overload for parsing SEQ.
    template <typename parsed_field_t>
    void parse_field(vtag_t<field::seq> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field::seq>(raw_record);

        if (raw_field != "*")
            parse_field_aux(raw_field, parsed_field); // reading into e.g. dna4 vector
    }

    //!\brief Overload for parsing QUAL.
    template <typename parsed_field_t>
    void parse_field(vtag_t<field::qual> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field::qual>(raw_record);

        if (raw_field != "*")
            parse_field_aux(raw_field, parsed_field); // reading into e.g. dna4 vector
    }

    //!\brief Overload for parsing the SAM tag dictionary.
    void parse_field(vtag_t<field::tags> const & /**/, map_io::sam_tag_dictionary & dictionary)
    {
        // we access the remaining fields from the iterator directly,
        // so we don't have to split the raw_record field again
        for (size_t i = 11; i < file_it->fields.size(); ++i)
            dictionary.parse_and_emplace(file_it->fields[i]);
    }

    //!\brief Overload for parsing the private data.
    void parse_field(vtag_t<field::_private> const & /**/, map_io::record_private_data & parsed_field)
    {
        parsed_field.header_ptr = &header;
    }
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
     * The options argument is typically bio::map_io::reader_options, but any object with a subset of similarly named
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

        std::string header_string;
        while (file_it != std::default_sentinel && file_it.peak() == '@')
        {
            ++file_it;
            ++line;
            header_string += file_it->line;
            header_string += "\n";
        }

        header = map_io::header{};
        header.read(header_string);
    }

    //!\brief Construct with only an input stream.
    format_input_handler(std::istream & str) : format_input_handler{str, int{}} {}
    //!\}

    //!\brief Return a reference to the header contained in the input handler.
    auto const & get_header() const { return header; }
};

} // namespace bio
