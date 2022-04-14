// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::map_io::tag_dictionary class and auxiliaries.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <map>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

#include <deque>
#include <seqan3/std/ranges>
#include <unordered_map>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>
#include <seqan3/utility/views/type_reduce.hpp>

#include <bio/exception.hpp>
#include <bio/misc.hpp>

namespace bio::map_io
{

/*!\brief The header used in mapping I/O.
 * \ingroup map_io
 * \details
 *
 * Each header line begins with the character `@` followed by one of the two-letter header record type codes
 * defined in this section. In the header, each line is tab-delimited and, apart from `@CO` lines, each data field
 * follows a format `TAG:VALUE` where TAG is a two-character string that defines the format and content of
 * VALUE. Thus header lines match `/^@(HD|SQ|RG|PG)(\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/` or are comment lines staring
 * with `@CO` followed by a tab and any character sequence.
 * Within each (non-`@CO`) header line, no field tag may appear more than once and the order in which the fields
 * appear is not significant.
 *
 * \sa https://samtools.github.io/hts-specs/SAMv1.pdf
 */
class header
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default construction.
    header()                           = default; //!< Defaulted.
    header(header const &)             = default; //!< Defaulted.
    header(header &&)                  = default; //!< Defaulted.
    ~header()                          = default; //!< Defaulted.
    header & operator=(header const &) = default; //!< Defaulted.
    header & operator=(header &&)      = default; //!< Defaulted.

    /*!\brief Construct from a range of reference ids.
     * \param[in] ref_ids The range over reference ids to redirect the pointer at.
     */
    template <typename ref_ids_type> // todo: restrict value type to be std::string_view constructible
    header(ref_ids_type & ref_ids) : reference_names_given_on_construction{true}
    {
        for (auto & id : ref_ids)
            reference_names.emplace_back(id);
    }

    /*!\brief Construct from a rvalue range of reference ids. Ids are copied over.
     * \param[in] ref_ids The range over reference ids to own.
     */
    template <typename ref_ids_type> // todo: restrict value type to be std::string_view constructible
    header(ref_ids_type && ref_ids) : reference_names_given_on_construction{true}
    {
        for (auto & id : ref_ids)
        {
            owned_reference_names.push_back(id);
            reference_names.emplace_back(owned_reference_names.back());
        }
    }
    //!\}

    //!\brief Defaulted three-way comparison.
    auto operator<=>(header const &) const = default;

    inline void read(std::string_view header_string);

private:
    //!\brief In case no reference information is given on construction, the reference names are stored here.
    std::deque<std::string> owned_reference_names;

    //!\brief The reference sequence names.
    std::vector<std::string_view> reference_names;

    //!\brief Additional information to the reference sequence (same ordering as `reference_names`).
    std::vector<std::tuple<int32_t, std::string>> reference_names_info{};

    //!\brief The mapping of reference name to position in the reference_names range and the rnames_info() range.
    std::unordered_map<std::string_view, int32_t> reference_name_to_pos{};

    //!\brief Whether reference sequence names were given to the header on construction.
    bool reference_names_given_on_construction{false};

    //!\brief Print a B.I.O warning message with current line number in diagnostic.
    /* [[noreturn]] compiler says this returns something...? */ void warning(auto const &... messages) const
    {
        // if (print_warnings)
        // {
        // seqan3::debug_stream << "[B.I.O sam format header warning in line " << line << "] "; // todo line numbers
        seqan3::debug_stream << "[B.I.O sam format header warning] ";
        (seqan3::debug_stream << ... << messages);
        seqan3::debug_stream << std::endl;
        // }
    }

public:
    /*!\name [HD] File-level meta data
     * \brief You can directly edit these member variables.
     * \{
     */
    std::string
      format_version{};       //!< [HD VN] The file format version. Note: this is overwritten by our formats on output.
    std::string sorting{};    //!< [HD SO] The sorting of the file. SAM: [unknown, unsorted, queryname, coordinate].
    std::string grouping{};   //!< [HD GO] The grouping of the file. SAM: [none, query, reference].
    std::string subsorting{}; //!< [HD SS] The sub-sorting of the file. SAM: [unknown, unsorted, queryname,
                              //!< coordinate]`(:[A-Za-z0-9_-]+)+`.
    //!\}

    /*!\name [SQ] Reference sequence dictionary
     * \brief You **CANNOT** directly edit these member variables. Please use the respective modifiers.
     * \{
     */

    /*!\brief [SQ SN] Reference sequence names
     *
     * \details
     *
     * This member function gives you access to the range of reference names.
     *
     * When reading a file, there are three scenarios:
     * 1) Reference id information is provided on construction. In this case, no copy is made but this function
     *    gives you a reference to the provided range. When reading the header or the records, their reference
     *    information will be checked against the given input.
     * 2) No reference information is provided on construction but the `@SQ` tags are present in the header.
     *    In this case, the reference id information is extracted from the header and this member function provides
     *    access to them. When reading the records, their reference id information will be checked against the header
     *    information.
     * 3) No reference information is provided on construction an no `@SQ` tags are present in the header.
     *    In this case, the reference information is parsed from the records field::ref_id and stored in the header.
     *    This member function then provides access to the unique list of reference names encountered in the records.
     */
    std::vector<std::string_view> const & rnames() { return reference_names; }

    /*!\brief [SQ LN,AH,AN,AS,M5,SP,UR] Reference sequence auxiliary information
     *
     * \details
     *
     * The reference information store the length (`@LN` tag) and
     * additional information of each reference sequence in the file. The record
     * must then store only the index of the reference.
     * The name and length information are required if the header is provided
     * and each reference sequence that is referred to in any of the records
     * must be present in the dictionary, otherwise a bio::format_error will
     * be thrown upon reading or writing a file.
     *
     * The additional information (2nd tuple entry) must model
     * the following formatting rules: The information is given in a tab separated
     * `TAG:VALUE` format, where TAG must be one of [AH, AN, AS, m5, SP, UR].
     * The following information and rules apply for each tag (taken from the SAM specs):
     *
     * * **AH:** Indicates that this sequence is an alternate locus. The value is the locus in the primary assembly for
     *           which this sequence is an alternative, in the format `chr:start-end`, `chr` (if known), or `*` (if
     *           unknown), where `chr` is a sequence in the primary assembly. Must not be present on sequences in the
     *           primary assembly.
     * * **AN:** Alternative reference sequence names. A comma-separated list of alternative names that tools may use
     *           when referring to this reference sequence. These alternative names are not used elsewhere within the
     *           SAM file; in  particular, they must not appear in alignment records’ RNAME or RNEXT fields. regular
     *           expression : `name (, name )*` where name is `[0-9A-Za-z][0-9A-Za-z*+.@ \|-]*`.
     * * **AS:** Genome assembly identifier.
     * * **M5:** MD5 checksum of the sequence.  See Section 1.3.1
     * * **SP:** Species.
     * * **UR:** URI of the sequence.  This value may start with one of the standard protocols, e.g http:  or ftp:. If
     *           it does not start with one of these protocols, it is assumed to be a file-system path
     */
    std::vector<std::tuple<int32_t, std::string>> const & rnames_info() { return reference_names_info; }

    //!\brief The mapping of reference name to position in the reference_names range and the rnames_info() range.
    std::unordered_map<std::string_view, int32_t> const & rname_to_pos() { return reference_name_to_pos; }

    //!\brief Append a new rname entry with it's length and extra information and update reference_name_to_pos.
    void push_back_rname(std::string_view const rname, int32_t const length, std::string_view const extra_info)
    {
        owned_reference_names.emplace_back(rname);
        reference_names.emplace_back(owned_reference_names.back());
        reference_names_info.emplace_back(length, extra_info);
        reference_name_to_pos[reference_names.back()] = reference_names.size() - 1;
    }
    //!\}

    /*!\name [RG] Read groups
     * \brief You can directly edit these member variables.
     * \{
     */
    /*!\brief The Read Group List
     *
     * \details
     *
     * The read group list stores the group id and
     * additional information of each read group in the file. The record
     * may store a RG tag information referencing one of the stored id's.
     * The id information is required if the \@RG header line is provided.
     *
     * The additional information (2nd tuple entry) for the SAM format must follow
     * the following formatting rules: The information is given in a tab separated
     * TAG:VALUE format, where TAG must be one of [AH, AN, AS, m5, SP, UR].
     * The following information and rules apply for each tag (taken from the SAM specs):
     *
     * * **BC:** Barcode sequence identifying the sample or library. This value is the expected barcode bases as read by
     *           the sequencing machine in the absence of errors. If there are several barcodes for the sample/library
     *           (e.g., one on each end of the template), the recommended implementation concatenates all the barcodes
     *           separating them with hyphens (`-`).
     * * **CN:** Name of sequencing center producing the read.
     * * **DS:** Description.  UTF-8 encoding may be used.
     * * **DT:** Date the run was produced (ISO8601 date or date/time).
     * * **FO:** Flow order. The array of nucleotide bases that correspond to the nucleotides used for each flow of each
     *           read. Multi-base flows are encoded in IUPAC format, and non-nucleotide flows by various other
     *           characters. Format : `/\*\|[ACMGRSVTWYHKDBN]+/`
     * * **KS:** The array of nucleotide bases that correspond to the key sequence of each read.
     * * **LB:** Library.
     * * **PG:** Programs used for processing the read group.
     * * **PI:** Predicted median insert size.
     * * **PL:** Platform/technology used to produce the reads. Valid values : CAPILLARY, LS454, ILLUMINA, SOLID,
     *           HELICOS, IONTORRENT, ONT, and PACBIO.
     * * **PM:** Platform model. Free-form text providing further details of the platform/technology used.
     * * **PU:** Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.
     * * **SM:** Sample. Use pool name where a pool is being sequenced.
     */
    std::vector<std::pair<std::string, std::string>> read_groups{};
    //!\}

    /*!\name [PG] Programm information
     * \brief You can directly edit these member variables.
     * \{
     */
    //!\brief Stores information of the program/tool that was used to create the file.
    struct program_info_t
    {
        std::string id{};                //!< A unique (file scope) id.
        std::string name{};              //!< The official name.
        std::string command_line_call{}; //!< The command line call that produces the file.
        std::string previous{};          //!< The id of the previous program if program calls were chained.
        std::string description{};       //!< A description of the program and/or program call.
        std::string version{};           //!< The program/tool version.
    };

    std::vector<program_info_t> program_infos{}; //!< The list of program information.
    //!\}

    /*!\name [CO] Comments
     * \brief You can directly edit these member variables.
     * \{
     */
    std::vector<std::string> comments{}; //!< The list of comments.
    //!\}
};

/*!\brief Reads the SAM header.
 * \param[in] header_string  The full header as a std::string_view.
 *
 * \throws bio::map_io::format_error if any unexpected character or format is encountered.
 *
 * \details
 *
 * Reading the header format is done according to the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
 *
 * The function throws a bio::map_io::format_error if any unknown tag was encountered. It will also fail if the format
 * is not in a correct state (e.g. required fields are not given), but throwing might occur downstream of the actual
 * error.
 */
void header::read(std::string_view header_string)
{
    auto stream_view = header_string | seqan3::views::single_pass_input;

    auto                it  = std::ranges::begin(stream_view);
    auto                end = std::ranges::end(stream_view);
    std::array<char, 2> raw_tag{};

    auto make_tag = [](uint8_t char1, uint8_t char2) constexpr
    {
        return static_cast<uint16_t>(char1) | (static_cast<uint16_t>(char2) << CHAR_BIT);
    };

    auto parse_and_make_tag = [&]()
    {
        raw_tag[0] = *it;
        ++it;
        raw_tag[1] = *it;
        ++it;
        return make_tag(raw_tag[0], raw_tag[1]);
    };

    auto skip_until = [&it](auto const & predicate)
    {
        while (!predicate(*it))
            ++it;
    };

    auto copy_tag_value_into = [&](auto & target)
    {
        skip_until(seqan3::is_char<':'>);
        ++it; // skip :
        while (!(seqan3::is_char<'\t'> || seqan3::is_char<'\n'>)(*it))
        {
            target.push_back(*it);
            ++it;
        }
    };

    auto copy_tag_into = [&raw_tag, &copy_tag_value_into](std::string & target)
    {
        // Some tags are not parsed individually. Instead, these are simply copied into a std::string.
        // Multiple tags must be separated by a `\t`, hence we prepend a tab to the string, except the first time.
        if (!target.empty())
            target.push_back('\t');
        target.push_back(raw_tag[0]);
        target.push_back(raw_tag[1]);
        target.push_back(':');
        copy_tag_value_into(target);
    };

    auto print_cerr_of_unspported_tag = [&it, &raw_tag](char const * const header_tag)
    { std::cerr << "Unsupported SAM header tag in @" << header_tag << ": " << raw_tag[0] << raw_tag[1] << '\n'; };

    while (it != end && seqan3::is_char<'@'>(*it))
    {
        ++it; // skip @

        switch (parse_and_make_tag())
        {
            case make_tag('H', 'D'): // HD (header) tag
                {
                    // All tags can appear in any order, VN is the only required tag
                    while (seqan3::is_char<'\t'>(*it))
                    {
                        ++it; // skip tab

                        switch (parse_and_make_tag())
                        {
                            case make_tag('V', 'N'): // parse required VN (version) tag
                                {
                                    copy_tag_value_into(this->format_version);
                                    break;
                                }
                            case make_tag('S', 'O'): // SO (sorting) tag
                                {
                                    copy_tag_value_into(this->sorting);
                                    break;
                                }
                            case make_tag('S', 'S'): // SS (sub-order) tag
                                {
                                    copy_tag_value_into(this->subsorting);
                                    break;
                                }
                            case make_tag('G', 'O'): // GO (grouping) tag
                                {
                                    copy_tag_value_into(this->grouping);
                                    break;
                                }
                            default: // unsupported header tag
                                {
                                    print_cerr_of_unspported_tag("HD");
                                    skip_until(seqan3::is_char<'\t'> || seqan3::is_char<'\n'>);
                                }
                        }
                    }
                    ++it; // skip newline

                    if (format_version.empty())
                        throw format_error{"The required VN tag in @HD is missing."};

                    break;
                }

            case make_tag('S', 'Q'): // SQ (sequence dictionary) tag
                {
                    std::string                      id;
                    std::tuple<int32_t, std::string> info{};
                    bool                             sequence_length_was_present{};

                    // All tags can appear in any order, SN and LN are required tags
                    while (seqan3::is_char<'\t'>(*it))
                    {
                        ++it; // skip tab

                        switch (parse_and_make_tag())
                        {
                            case make_tag('S', 'N'): // parse required SN (sequence name) tag
                                {
                                    copy_tag_value_into(id);
                                    break;
                                }
                            case make_tag('L', 'N'): // parse required LN (length) tag
                                {
                                    sequence_length_was_present = true;
                                    ++it; // skip :
                                    auto res = std::from_chars(&(*it), &(header_string.back()), get<0>(info));
                                    if (res.ec != std::errc{})
                                        throw format_error{"LN tag could not be parsed correctly."};
                                    skip_until(seqan3::is_char<'\t'> || seqan3::is_char<'\n'>);
                                    break;
                                }
                            default: // Any other tag
                                {
                                    copy_tag_into(std::get<1>(info));
                                }
                        }
                    }
                    ++it; // skip newline

                    if (id.empty())
                        throw format_error{"The required SN tag in @SQ is missing."};
                    if (!sequence_length_was_present)
                        throw format_error{"The required LN tag in @SQ is missing."};

                    if (!reference_names_given_on_construction)
                    { // Reference information was not given by the user but needs to be filled from the header
                        push_back_rname(id, std::get<0>(info), std::get<1>(info));
                    }
                    else
                    { // If reference information was given we check them against the once in the header
                        auto id_it = reference_name_to_pos.find(id);

                        if (id_it == reference_name_to_pos.end())
                        {
                            warning("The reference sequence name \"",
                                    id,
                                    "\" was present in the header but not in the "
                                    "user provided rnames.");

                            push_back_rname(id, std::get<0>(info), std::get<1>(info));
                        }
                        else
                        {
                            assert(id_it->second < static_cast<decltype(id_it->second)>(reference_names_info.size()));

                            if (std::get<0>(reference_names_info[id_it->second]) != std::get<0>(info))
                                warning("Provided and header-based reference length differ for rname :\"", id, "\".");

                            reference_names_info[id_it->second] = info; // update rname information
                        }
                    }
                    break;
                }

            case make_tag('R', 'G'): // RG (read group) tag
                {
                    std::pair<std::string, std::string> tmp{};

                    // All tags can appear in any order, SN and LN are required tags
                    while (seqan3::is_char<'\t'>(*it))
                    {
                        ++it; // skip tab

                        switch (parse_and_make_tag())
                        {
                            case make_tag('I', 'D'): // parse required ID tag
                                {
                                    copy_tag_value_into(get<0>(tmp));
                                    break;
                                }
                            default: // Any other tag
                                {
                                    copy_tag_into(get<1>(tmp));
                                }
                        }
                    }
                    ++it; // skip newline

                    if (get<0>(tmp).empty())
                        throw format_error{"The required ID tag in @RG is missing."};

                    read_groups.emplace_back(std::move(tmp));
                    break;
                }

            case make_tag('P', 'G'): // PG (program) tag
                {
                    program_info_t tmp{};

                    // All tags can appear in any order, ID is the only required tag
                    while (seqan3::is_char<'\t'>(*it))
                    {
                        ++it; // skip tab

                        switch (parse_and_make_tag())
                        {
                            case make_tag('I', 'D'): // read required ID tag
                                {
                                    copy_tag_value_into(tmp.id);
                                    break;
                                }
                            case make_tag('P', 'N'): // PN (program name) tag
                                {
                                    copy_tag_value_into(tmp.name);
                                    break;
                                }
                            case make_tag('P', 'P'): // PP (previous program) tag
                                {
                                    copy_tag_value_into(tmp.previous);
                                    break;
                                }
                            case make_tag('C', 'L'): // CL (command line) tag
                                {
                                    copy_tag_value_into(tmp.command_line_call);
                                    break;
                                }
                            case make_tag('D', 'S'): // DS (description) tag
                                {
                                    copy_tag_value_into(tmp.description);
                                    break;
                                }
                            case make_tag('V', 'N'): // VN (version) tag
                                {
                                    copy_tag_value_into(tmp.version);
                                    break;
                                }
                            default: // unsupported header tag
                                {
                                    print_cerr_of_unspported_tag("PG");
                                    skip_until(seqan3::is_char<'\t'> || seqan3::is_char<'\n'>);
                                }
                        }
                    }
                    ++it; // skip newline

                    if (tmp.id.empty())
                        throw format_error{"The required ID tag in @PG is missing."};

                    this->program_infos.push_back(std::move(tmp));
                    break;
                }

            case make_tag('C', 'O'): // CO (comment) tag
                {
                    ++it; // skip tab
                    std::string tmp{};
                    while (!seqan3::is_char<'\n'>(*it))
                    {
                        tmp.push_back(*it);
                        ++it;
                    }
                    ++it; // skip newline
                    comments.emplace_back(std::move(tmp));
                    break;
                }

            default:
                throw format_error{std::string{"Illegal SAM header tag starting with:"} + *it};
        }
    }
}

//!\brief A datastructure that contains private data of sam IO records.
//!\ingroup var_io
struct record_private_data
{
    //!\privatesection
    //!\brief Pointer to the header
    header const * header_ptr = nullptr;

    friend bool operator==(record_private_data const &, record_private_data const &) = default; //!< Defaulted.
};

} // namespace bio::map_io
