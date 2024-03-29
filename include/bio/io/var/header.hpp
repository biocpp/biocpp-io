// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::io::var::header and various auxiliary classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <algorithm>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <bio/ranges/container/dictionary.hpp>
#include <bio/ranges/container/small_string.hpp>

#include <bio/io/detail/charconv.hpp>
#include <bio/io/detail/views_eager_split.hpp>
#include <bio/io/exception.hpp>
#include <bio/io/misc.hpp>
#include <bio/io/var/misc.hpp>

namespace bio::io::var
{

/*!\brief Scoped (but weakly typed) enum for "Number" special values in bio::io::var::header INFO fields.
 * \ingroup var
 * \details
 *
 * ### Example
 *
 * TODO externalise
 * ```cpp
 * int32_t number = 3;
 * number = bio::io::header_number::A;
 * ```
 */
struct header_number
{
    //!\brief The anonymous enum.
    enum : int32_t
    {
        A   = -1, //!< One value per alternate allele.
        R   = -2, //!< One value for each possible allele (including ref) -> A + 1.
        G   = -3, //!< One value per Genotype.
        dot = -4, //!< Unknown, unspecified or unbounded.
    };

    //!\brief Convert a "header number" to a string.
    static inline std::string to_string(int32_t const n)
    {
        switch (n)
        {
            case A:
                return "A";
            case R:
                return "R";
            case G:
                return "G";
            case dot:
                return ".";
            default:
                return std::to_string(n);
        }
    }
};

/*!\brief The header used in variant I/O.
 * \ingroup var
 * \details
 *
 * TODO
 */
class header
{
public:
    // TODO Change the following into a ranges::dictionary, as well
    //!\brief The dictionary for the non-standard fields in a header entry.
    using other_fields_t = std::map<std::string, std::string, std::ranges::less>;

    /*!\name Header entry types
     * \{
     */
    //!\brief Type of the contig field header line.
    struct contig_t
    {
        int64_t        length = -1;    //!< Length of the contig (-1 if absent).
        other_fields_t other_fields{}; //!< Other entries.
        int32_t        idx = -1;       //!< The numeric ID.

        //!\brief Defaulted three-way comparisons.
        auto operator<=>(contig_t const &) const = default;
    };

    //!\brief Type of a INFO field header line.
    struct info_t
    {
        int32_t        number{};       //!< Number of values, see also bio::io::var::header_number.
        std::string    type{};         //!< Type of the field.
        type_enum      type_id{};      //!< Type of the field as bio::var::type_enum.
        std::string    description{};  //!< Description.
        other_fields_t other_fields{}; //!< Other entries.
        int32_t        idx = -1;       //!< The numeric ID.

        //!\brief Defaulted three-way comparisons.
        auto operator<=>(info_t const &) const = default;
    };

    //!\brief Type of a FILTER field header line.
    struct filter_t
    {
        std::string    description{};  //!< Description.
        other_fields_t other_fields{}; //!< Other entries.
        int32_t        idx = -1;       //!< The numeric ID.

        //!\brief Defaulted three-way comparisons.
        auto operator<=>(filter_t const &) const = default;
    };

    using format_t = info_t; //!< Type of a FORMAT field header line.

    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default construction.
    header() { add_pass_entry(); }
    header(header const &)             = default; //!< Defaulted.
    header(header &&)                  = default; //!< Defaulted.
    ~header()                          = default; //!< Defaulted.
    header & operator=(header const &) = default; //!< Defaulted.
    header & operator=(header &&)      = default; //!< Defaulted.

    //!\brief Construct from a header given as plaintext.
    explicit header(std::string_view plaintext_header)
    {
        add_pass_entry();

        if (plaintext_header.ends_with("\r\n"))
            plaintext_header = plaintext_header.substr(0, plaintext_header.size() - 2);
        else if (plaintext_header.ends_with("\n"))
            plaintext_header = plaintext_header.substr(0, plaintext_header.size() - 1);

        for (std::string_view const line : plaintext_header | io::detail::eager_split('\n'))
            parse_line(line);

        //TODO do we want to validate?
    }
    //!\}

    //!\brief Defaulted three-way comparison.
    auto operator<=>(header const &) const = default;

    /*!\name Header entries
     * \brief You can directly edit these member variables.
     * \{
     */
    std::string                               file_format = "VCFv4.3"; //!< The file format version.
    ranges::dictionary<std::string, filter_t> filters;                 //!< Header lines describing FILTER fields.
    ranges::dictionary<std::string, info_t>   infos;                   //!< Header lines describing INFO fields.
    ranges::dictionary<std::string, format_t> formats;                 //!< Header lines describing FORMAT fields.
    ranges::dictionary<std::string, contig_t> contigs;                 //!< Header lines describing contigs.
    std::vector<std::string>                  other_lines;             //!< Any other lines in the header.
    std::vector<std::string>                  column_labels;           //!< Standard column labels and sample names.
    //!\}

    /*!\name IDX maps
     * \{
     */
    //!\brief Hash map of IDX to string.
    using idx_to_string_map_t = std::unordered_map<int32_t, std::string>;

    //!\brief Global string to IDX mapping (filter, info, format).
    idx_to_string_map_t const & idx_to_string_map() const { return idx_to_string_map_; }
    //!\brief Global string to IDX mapping (contig).
    idx_to_string_map_t const & contig_idx_to_string_map() const { return contig_idx_to_string_map_; }

    //!\brief The largest IDX value used (filter, info, format).
    int32_t max_idx() const { return max_other_idx_; }
    //!\brief The largest IDX value used (contig).
    int32_t max_contig_idx() const { return max_contig_idx_; }

    /*!\brief Add missing IDX fields to header entries and ensure that everything has proper hash-entries.
     *
     * \details
     *
     * Does the following things:
     *
     *   1. ensure that "PASS" filter entry is present as first filter entry.
     *   2. assign a valid IDX value to all header entries that currently have IDX of -1.
     *   3. generate new idx-to-string maps.
     *
     * It does not:
     *
     *   * change existing IDX values of header entries (that are not -1).
     *   * re-use obsolete IDX values; call #idx_clear() before to achieve this.
     *
     * If you want to start with new IDX numbering, calls #idx_clear().
     */
    void idx_update()
    {
        idx_to_string_map_.clear();
        contig_idx_to_string_map_.clear();

        bool has_pass = false;
        for (auto && [id, filter] : filters)
        {
            if (id == "PASS")
            {
                has_pass              = true;
                filter.idx            = 0;
                idx_to_string_map_[0] = "PASS";
            }
            else
            {
                fix_idx(filter.idx, id);
            }
        }

        if (!has_pass)
            add_pass_entry();

        for (auto && [id, info] : infos)
            fix_idx(info.idx, id);

        for (auto && [id, format] : formats)
            fix_idx(format.idx, id);

        for (auto && [id, contig] : contigs)
            fix_contig_idx(contig.idx, id);
    }

    //!\brief Clear the IDX values from all header entries (sets them to -1).
    void idx_clear()
    {
        idx_to_string_map_.clear();
        contig_idx_to_string_map_.clear();

        max_contig_idx_ = -1;
        max_other_idx_  = 0;

        for (auto && [id, filter] : filters)
        {
            if (id == "PASS")
                filter.idx = 0;
            else
                filter.idx = -1;
        }

        for (info_t & info : infos | std::views::elements<1>)
            info.idx = -1;

        for (format_t & format : formats | std::views::elements<1>)
            format.idx = -1;

        for (contig_t & contig : contigs | std::views::elements<1>)
            contig.idx = -1;
    }

    /*!\brief Ensure the validity of the header entries.
     *
     * \details
     *
     * This finds the following kinds of errors:
     *
     *  * Non-unique ID or IDX values in header entries.
     *  * IDX-values that are -1.
     *  * Missing "PASS" entry in filters.
     *  * Strings in "other_lines" that actually should have been parsed.
     *  * Inconsistencies in the hash tables.
     *
     */
    void idx_validate() const
    {
        auto check_idx = [](int32_t const idx, std::string_view const id, auto const & map)
        {
            if (idx == -1)
            {
                throw format_error{"VCF Header entry with ID ",
                                   id,
                                   " has no IDX value set.\n"
                                   "Call idx_update() on the header."};
            }

            if (auto it = map.find(idx); it == map.end())
            {
                throw format_error{"VCF Header entry with ID ",
                                   id,
                                   " and IDX ",
                                   idx,
                                   " does not appear in reverse IDX map.\nCall idx_update() on the header."};
            }
            else if (std::string_view const stored_id = std::get<1>(*it); stored_id != id)
            {
                throw format_error{"VCF Header entry with ID ",
                                   id,
                                   " and IDX ",
                                   idx,
                                   " has different reverse ID lookup string: ",
                                   stored_id,
                                   "\nCall idx_clear() and "
                                   "idx_update() on the header."};
            }
        };

        bool has_pass = false;
        for (auto && [id, filter] : filters)
        {
            if (id == "PASS")
                has_pass = true;

            check_idx(filter.idx, id, idx_to_string_map());
        }

        if (!has_pass)
            throw format_error{"No VCF Header entry for the PASS-filter.\nCall idx_update() on the header."};

        for (auto && [id, info] : infos)
            check_idx(info.idx, id, idx_to_string_map());

        for (auto && [id, format] : formats)
            check_idx(format.idx, id, idx_to_string_map());

        for (auto && [id, contig] : contigs)
            check_idx(contig.idx, id, contig_idx_to_string_map());

        for (std::string_view const line : other_lines)
        {
            auto check = [line](std::string_view const keyword)
            {
                if (line.starts_with(keyword))
                {
                    throw format_error{"other_line in VCF Header entry is actually a ",
                                       keyword.substr(0, keyword.size() - 1), // remove "="
                                       "line.\nAdd data to the respective member instead."};
                }
            };

            check("FILTER=");
            check("INFO=");
            check("FORMAT=");
            check("contig=");
        }
    }
    //!\}

    /*!\name Convert to plaintext ("raw") header
     * \{
     */
    //!\brief Converts the header to plaintext (includes IDX entries).
    std::string to_plaintext() const { return to_plaintext_impl(true); }
    //!\brief Converts the header to plaintext (excludes IDX entries).
    std::string to_plaintext_without_idx() const { return to_plaintext_impl(false); }
    //!\}

private:
    //!\brief Whether the first line was read successfully.
    bool file_format_read = false;

    //!\brief Add implicit PASS filter.
    void add_pass_entry()
    {
        filters.emplace_back("PASS", filter_t{"\"All filters passed\"", {}, 0});
        idx_to_string_map_[0] = "PASS";
    }

    /*!\brief Set the correct IDX on an entry and updates the corresponding hash table.
     * \param[in,out] idx  The IDX value.
     * \param[in] id       The string ID.
     */
    void fix_idx(int32_t & idx, std::string const & id)
    {
        if (idx == -1)
        {
            // reverse lookup
            if (auto it =
                  std::ranges::find_if(idx_to_string_map_, [&](auto const & val) { return std::get<1>(val) == id; });
                it == idx_to_string_map_.end())
            {
                ++max_other_idx_;
                idx = max_other_idx_;
            }
            else
            {
                idx = std::get<0>(*it);
            }
        }

        if (auto it = idx_to_string_map_.find(idx); it == idx_to_string_map_.end())
            idx_to_string_map_[idx] = id;
        else if (std::string_view id_ = std::get<1>(*it); id_ != id)
            throw format_error{"Couldn't map IDX ", idx, " to ID ", id, ", because already mapped to ", id_, "."};
    }

    /*!\brief Set the correct IDX on a contig entry and updates the corresponding hash table.
     * \param[in,out] idx  The IDX value.
     * \param[in] id       The string ID.
     */
    void fix_contig_idx(int32_t & idx, std::string const & id)
    {
        if (idx == -1)
        {
            // reverse lookup
            if (auto it = std::ranges::find_if(contig_idx_to_string_map_,
                                               [&](auto const & val) { return std::get<1>(val) == id; });
                it == contig_idx_to_string_map_.end())
            {
                ++max_contig_idx_;
                idx = max_contig_idx_;
            }
            else
            {
                idx = std::get<0>(*it);
            }
        }

        if (auto it = contig_idx_to_string_map_.find(idx); it == contig_idx_to_string_map_.end())
            contig_idx_to_string_map_[idx] = id;
        else if (std::string_view id_ = std::get<1>(*it); id_ != id)
            throw format_error{"Couldn't map IDX ", idx, " to ID ", id, ", because already mapped to ", id_, "."};
    }

    /*!\name Advanced data fields
     * \brief You don't have to set these manually when creating a bio::io::var::header.
     * \{
     */
    idx_to_string_map_t idx_to_string_map_;        //!< Global IDX to string mapping (filter, info, format).
    idx_to_string_map_t contig_idx_to_string_map_; //!< Global IDX to string mapping (contig).

    int32_t max_other_idx_  = 0;  //!< The highest IDX value in use (defaults to 0, because PASS is used).
    int32_t max_contig_idx_ = -1; //!< The highest contig IDX value in use (defaults to -1, because none is used).
    //!\}

    /*!\name Functions for converting to text
     * \{
     */
    //!\brief Implementation function for creating the plaintext header.
    std::string to_plaintext_impl(bool const with_idx) const
    {
        std::string raw_data;

        constexpr auto quote_wrap = [](std::string in)
        {
            if (in.size() == 0)
                in = "\"\"";
            else if (in.front() != '\"')
                in.insert(in.begin(), '\"');

            if (in.size() == 1 || in.back() != '\"')
                in.push_back('\"');

            return in;
        };

        /* first of all, the structured data is transformed into plaintext */

        // file format
        ((raw_data += "##fileformat=") += file_format) += "\n";

        // filters
        for (auto && [id, filter] : filters)
        {
            (raw_data += "##FILTER=<ID=") += id;
            (raw_data += ",Description=") += quote_wrap(static_cast<std::string>(filter.description));

            for (auto [key, value] : filter.other_fields)
                (((raw_data += ",") += key) += "=") += value;

            if (with_idx)
                (raw_data += ",IDX=") += std::to_string(filter.idx); // TODO replace with std::to_chars

            raw_data += ">\n";
        }

        // TODO: think about if/which other_fields-value to quote_wrap
        //  infos
        for (auto && [id, info] : infos)
        {
            (raw_data += "##INFO=<ID=") += id;
            (raw_data += ",Number=") += header_number::to_string(info.number);
            (raw_data += ",Type=") += unparse_type(info.type, info.type_id);
            (raw_data += ",Description=") += quote_wrap(static_cast<std::string>(info.description));

            for (auto [key, value] : info.other_fields)
                (((raw_data += ",") += key) += "=") += value;

            if (with_idx)
                (raw_data += ",IDX=") += std::to_string(info.idx); // TODO replace with std::to_chars

            raw_data += ">\n";
        }

        // formats
        for (auto && [id, format] : formats)
        {
            (raw_data += "##FORMAT=<ID=") += id;
            (raw_data += ",Number=") += header_number::to_string(format.number);
            (raw_data += ",Type=") += unparse_type(format.type, format.type_id);
            (raw_data += ",Description=") += quote_wrap(static_cast<std::string>(format.description));

            for (auto [key, value] : format.other_fields)
                (((raw_data += ",") += key) += "=") += value;

            if (with_idx)
                (raw_data += ",IDX=") += std::to_string(format.idx); // TODO replace with std::to_chars

            raw_data += ">\n";
        }

        // contigs
        for (auto && [id, contig] : contigs)
        {
            (raw_data += "##contig=<ID=") += id;
            if (contig.length != -1)
                (raw_data += ",length=") += std::to_string(contig.length); // TODO replace with std::to_chars

            for (auto [key, value] : contig.other_fields)
                (((raw_data += ",") += key) += "=") += value;

            if (with_idx)
                (raw_data += ",IDX=") += std::to_string(contig.idx); // TODO replace with std::to_chars

            raw_data += ">\n";
        }

        // other lines
        for (auto const & line : other_lines)
            ((raw_data += "##") += line) += "\n";

        // column labals
        // TODO check that for 8/9 value are the required ones!
        raw_data += "#CHROM";
        for (size_t i = 1; i < column_labels.size(); ++i)
            (raw_data += "\t") += column_labels[i];
        raw_data += "\n";

        return raw_data;
    }

    //!\brief Turn bio::io::type_enum into string.
    static std::string_view unparse_type(std::string_view const type, type_enum const type_id)
    {
        std::string_view ret;

        switch (type_id)
        {
            case type_enum::int8:
            case type_enum::vector_of_int8:
            case type_enum::int16:
            case type_enum::vector_of_int16:
            case type_enum::int32:
            case type_enum::vector_of_int32:
                ret = "Integer";
                break;
            case type_enum::float32:
            case type_enum::vector_of_float32:
                ret = "Float";
                break;
            case type_enum::char8:
                ret = "Character";
                break;
            case type_enum::string:
                if (type == "Character") // multiple characters
                {
                    ret = "Character";
                    break;
                }
                [[fallthrough]];
            case type_enum::vector_of_string:
                ret = "String";
                break;
            case type_enum::flag:
                ret = "Flag";
                break;
            default:
                throw format_error{"Illegal type_id in INFO or FILTER header line."};
        }

        if (!type.empty() && type != ret)
            throw format_error{"Type string and type_id mismatch."};

        return ret;
    }
    //!\}

    /*!\name Functions for converting from text
     * \{
     */
    //!\brief Main function for interpreting a line (delegates to other functions).
    void parse_line(std::string_view const l)
    {
        if (file_format_read == false)
        {
            if (l.starts_with("##fileformat="))
            {
                file_format      = l.substr(13);
                file_format_read = true;
            }
            else
            {
                throw format_error{"File does not begin with \"##fileformat\"."};
            }
        }
        else if (l.starts_with("##fileformat="))
        {
            throw format_error{"File has two lines that begin with \"##fileformat\"."};
        }
        else if (l.starts_with("##INFO="))
        {
            parse_info_or_format_line(strip_angular_brackets(l.substr(7)), true);
        }
        else if (l.starts_with("##FILTER="))
        {
            parse_filter_line(strip_angular_brackets(l.substr(9)));
        }
        else if (l.starts_with("##FORMAT="))
        {
            parse_info_or_format_line(strip_angular_brackets(l.substr(9)), false);
        }
        else if (l.starts_with("##contig="))
        {
            parse_contig_line(strip_angular_brackets(l.substr(9)));
        }
        else if (l.starts_with("#CHROM"))
        {
            parse_column_labels_line(l.substr(1));
        }
        else if (l.starts_with("##"))
        {
            other_lines.push_back(static_cast<std::string>(l.substr(2)));
        }
        else
        {
            throw format_error{"Plaintext header contains lines that don't start with \"##\" or \"#CHROM\"."};
        }
    }

    //!\brief Parse an INFO or FORMAT line.
    void parse_info_or_format_line(std::string_view const l, bool const is_info)
    {
        info_t new_entry;
        new_entry.other_fields = to_dictionary(l);

        /* ID */
        auto id = new_entry.other_fields.extract("ID");
        if (id.empty())
            throw format_error{"INFO or FORMAT line does not contain ID field."};

        std::string id_str = std::move(id.mapped());

        /* Number */
        auto number = new_entry.other_fields.extract("Number");
        if (number.empty())
            throw format_error{"INFO or FORMAT line does not contain Number field."};
        else
            new_entry.number = parse_number(number.mapped());

        /* Type */
        auto type = new_entry.other_fields.extract("Type");
        if (type.empty())
        {
            throw format_error{"INFO or FORMAT line does not contain Type field."};
        }
        else
        {
            new_entry.type    = type.mapped();
            new_entry.type_id = parse_type(type.mapped(), new_entry.number);
        }
        if (auto it = new_entry.other_fields.find("IntegerBits"); it != new_entry.other_fields.end())
        {
            std::string_view number = strip_quotes(it->second);

            if (number == "8")
            {
                switch (new_entry.type_id)
                {
                    case type_enum::int8:
                    case type_enum::int16:
                    case type_enum::int32:
                        new_entry.type_id = type_enum::int8;
                        break;
                    case type_enum::vector_of_int8:
                    case type_enum::vector_of_int16:
                    case type_enum::vector_of_int32:
                        new_entry.type_id = type_enum::vector_of_int8;
                        break;
                    default:
                        break;
                }
            }
            else if (number == "16")
            {
                switch (new_entry.type_id)
                {
                    case type_enum::int8:
                    case type_enum::int16:
                    case type_enum::int32:
                        new_entry.type_id = type_enum::int16;
                        break;
                    case type_enum::vector_of_int8:
                    case type_enum::vector_of_int16:
                    case type_enum::vector_of_int32:
                        new_entry.type_id = type_enum::vector_of_int16;
                        break;
                    default:
                        break;
                }
            }
            // integer fields are set to 32 by default, so no need to handle this case here
            // if the number is something else or the string is not a number, we also don't
            // do anything and just assume 32Bit Integer
        }

        /* Description */
        auto description = new_entry.other_fields.extract("Description");
        if (description.empty())
            throw format_error{"INFO or FORMAT line does not contain Description field."};
        else
            new_entry.description = description.mapped();

        /* IDX */
        auto idx = new_entry.other_fields.extract("IDX");
        if (!idx.empty())
            io::detail::string_to_number(idx.mapped(), new_entry.idx);
        max_other_idx_ = std::max(max_other_idx_, new_entry.idx);
        fix_idx(new_entry.idx, id_str);

        if (is_info)
        {
            if (infos.contains(id_str))
                throw format_error{"Duplicate INFO ID \"", id_str, "\" in HEADER."};

            infos.emplace_back(std::move(id_str), std::move(new_entry));
        }
        else
        {
            if (formats.contains(id_str))
                throw format_error{"Duplicate FORMAT ID \"", id_str, "\" in HEADER."};

            formats.emplace_back(std::move(id_str), std::move(new_entry));
        }
    }

    //!\brief Parse a FILTER line.
    void parse_filter_line(std::string_view const l)
    {
        filter_t new_entry;
        new_entry.other_fields = to_dictionary(l);

        /* ID */
        auto id = new_entry.other_fields.extract("ID");
        if (id.empty())
            throw format_error{"FILTER line does not contain ID field."};

        std::string id_str = std::move(id.mapped());

        /* Description */
        auto description = new_entry.other_fields.extract("Description");
        if (description.empty())
            throw format_error{"FILTER line does not contain Description field."};
        else
            new_entry.description = description.mapped();

        /* IDX */
        auto idx = new_entry.other_fields.extract("IDX");
        if (!idx.empty())
            io::detail::string_to_number(idx.mapped(), new_entry.idx);
        max_other_idx_ = std::max(max_other_idx_, new_entry.idx);
        fix_idx(new_entry.idx, id_str);

        // PASS line was added by us before and is now swapped with user-provided
        if (filters.size() > 0 && get<0>(filters.front()) == "PASS" && id_str == "PASS")
        {
            get<1>(filters[0]) = new_entry;
        }
        else
        {
            if (filters.contains(id_str))
                throw format_error{"Duplicate FILTER ID \"", id_str, "\" in HEADER."};

            filters.emplace_back(std::move(id_str), std::move(new_entry));
        }
    }

    //!\brief Parse a CONTIG line.
    void parse_contig_line(std::string_view const l)
    {
        contig_t new_entry;
        new_entry.other_fields = to_dictionary(l);

        /* ID */
        auto id = new_entry.other_fields.extract("ID");
        if (id.empty())
            throw format_error{"FILTER line does not contain ID field."};

        std::string id_str = std::move(id.mapped());

        /* Length */
        auto length = new_entry.other_fields.extract("length");
        if (!length.empty())
            io::detail::string_to_number(length.mapped(), new_entry.length);

        /* IDX */
        auto idx = new_entry.other_fields.extract("IDX");
        if (idx.empty())
            new_entry.idx = ++max_contig_idx_;
        else
            io::detail::string_to_number(idx.mapped(), new_entry.idx);
        max_contig_idx_ = std::max(max_contig_idx_, new_entry.idx);
        fix_contig_idx(new_entry.idx, id_str);

        if (contigs.contains(id_str))
            throw format_error{"Duplicate CONTIG ID \"", id_str, "\" in HEADER."};

        contigs.emplace_back(std::move(id_str), std::move(new_entry));
    }

    //!\brief Parse the line with column labels / sample names.
    void parse_column_labels_line(std::string_view const l)
    {
        for (std::string_view field : l | io::detail::eager_split('\t'))
            column_labels.push_back(static_cast<std::string>(field));
    }

    //!\brief Return a substring from the argument that does not contain enclosing angular brackets.
    static inline std::string_view strip_angular_brackets(std::string_view const in)
    {
        if (in.size() < 2 || in.front() != '<' || in.back() != '>')
            throw format_error{"Structured line does not contain \"<\" and \">\" at right places."};
        return in.substr(1, in.size() - 2);
    }

    //!\brief Return a substring from the argument that does not contain enclosing quotes (if present).
    static inline std::string_view strip_quotes(std::string_view const in)
    {
        return (in.size() < 2 || in.front() != '"' || in.back() != '"') ? in : in.substr(1, in.size() - 2);
    }

    //!\brief Turn a string into a bio::io::var::header_number.
    static int32_t parse_number(std::string_view const in)
    {
        switch (in[0])
        {
            case 'A':
                return header_number::A;
            case 'R':
                return header_number::R;
            case 'G':
                return header_number::G;
            case '.':
                return header_number::dot;
            default:
                {
                    int32_t ret = 0;
                    io::detail::string_to_number(in, ret);
                    return ret;
                }
        }
        return header_number::dot;
    }

    /*!\brief Turn a string into bio::io::var::type_enum.
     * \param[in] in        The input string.
     * \param[in] number    The accompanying number value from the header.
     * \return The dynamic type id.
     */
    static type_enum parse_type(std::string_view const in, int32_t const number)
    {
        type_enum ret{};

        if (in == "Flag")
        {
            ret = type_enum::flag;
            if (number != 0)
                throw format_error{"Flags must always have number 0 in header."};
            return ret;
        }

        if (number == 0)
            throw format_error{"Only flags may have number 0 in header."};

        if (in == "Integer")
        {
            if (number == 1)
                ret = type_enum::int32;
            else
                ret = type_enum::vector_of_int32;
        }
        else if (in == "Float")
        {
            if (number == 1)
                ret = type_enum::float32;
            else
                ret = type_enum::vector_of_float32;
        }
        else if (in == "Character")
        {
            if (number == 1)
                ret = type_enum::char8;
            else
                ret = type_enum::string;
        }
        else if (in == "String")
        {
            if (number == 1)
                ret = type_enum::string;
            else
                ret = type_enum::vector_of_string;
        }
        else
        {
            throw format_error{"Cannot convert the following string to a type identifier: ", in};
        }

        return ret;
    }

    //!\brief Turns a string containing value-pairs into a dictionary.
    static other_fields_t to_dictionary(std::string_view const value_pairs)
    {
        other_fields_t ret;

        for (std::string_view const pair : value_pairs | io::detail::eager_split(',', true))
        {
            auto pair_split = pair | io::detail::eager_split('=');
            auto it1        = pair_split.begin();
            auto it2        = std::ranges::next(it1);
            auto it3        = std::ranges::next(it2); // TODO whats going on here?

            if (it1 == std::default_sentinel || it2 == std::default_sentinel) //|| it3 != std::default_sentinel)
            {
                throw format_error{"Could not parse the following string into a dictionary: ", pair};
            }

            ret[static_cast<std::string>(*it1)] = static_cast<std::string>(*it2);
        }

        return ret;
    }
    //!\}
};

// clang-format off
/*!\brief A compile-time mapping of INFO strings to bio::io::var::type_enum.
 *!\ingroup var
 *
 * This mapping allows accessing bio::io::var::info_variant_shallow and bio::io::var::info_variant_deep by e.g.
 * `get<"AA">(variant)` instead of `get<static_cast<size_t>(type_enum::string)>(variant)`.
 *
 * The predefined mappings are as defined in "Table 1" of the
 * [VCF specifiction](https://samtools.github.io/hts-specs/VCFv4.3.pdf), and also shown in bio::io::var::reserved_infos.
 *
 * ### Customisation point
 *
 * This variable template is a customisation point, which means that you can specialise it for your own strings.
 *
 * TODO add example snippet
 */
template <ranges::small_string str>
inline constexpr meta::ignore_t info_key2type_enum{};

//!\cond
template <> inline constexpr type_enum info_key2type_enum<"AA"       > = type_enum::string;
template <> inline constexpr type_enum info_key2type_enum<"AC"       > = type_enum::vector_of_int32;
template <> inline constexpr type_enum info_key2type_enum<"AD"       > = type_enum::vector_of_int32;
template <> inline constexpr type_enum info_key2type_enum<"ADF"      > = type_enum::vector_of_int32;
template <> inline constexpr type_enum info_key2type_enum<"ADR"      > = type_enum::vector_of_int32;
template <> inline constexpr type_enum info_key2type_enum<"AF"       > = type_enum::vector_of_float32;
template <> inline constexpr type_enum info_key2type_enum<"AN"       > = type_enum::int32;
template <> inline constexpr type_enum info_key2type_enum<"BQ"       > = type_enum::float32;
template <> inline constexpr type_enum info_key2type_enum<"CIGAR"    > = type_enum::vector_of_string;
template <> inline constexpr type_enum info_key2type_enum<"DB"       > = type_enum::flag;
template <> inline constexpr type_enum info_key2type_enum<"DP"       > = type_enum::int32;
template <> inline constexpr type_enum info_key2type_enum<"END"      > = type_enum::int32;
template <> inline constexpr type_enum info_key2type_enum<"H2"       > = type_enum::flag;
template <> inline constexpr type_enum info_key2type_enum<"H3"       > = type_enum::flag;
template <> inline constexpr type_enum info_key2type_enum<"MQ"       > = type_enum::float32;
template <> inline constexpr type_enum info_key2type_enum<"MQ0"      > = type_enum::int32;
template <> inline constexpr type_enum info_key2type_enum<"NS"       > = type_enum::int32;
template <> inline constexpr type_enum info_key2type_enum<"SB"       > = type_enum::vector_of_int32;
template <> inline constexpr type_enum info_key2type_enum<"SOMATIC"  > = type_enum::flag;
template <> inline constexpr type_enum info_key2type_enum<"VALIDATED"> = type_enum::flag;
template <> inline constexpr type_enum info_key2type_enum<"1000G"    > = type_enum::flag;
//!\endcond

//!\brief A table of reserved INFO entries.
//!\ingroup var
inline ranges::dictionary<std::string_view, header::info_t> const reserved_infos
{
    {"AA",       {               1, "String",  type_enum::string,            "\"Ancestral allele\""}},
    {"AC",       {header_number::A, "Integer", type_enum::vector_of_int32,   "\"Allele count in genotypes, for each ALT allele, in the same order as listed\""}},
    {"AD",       {header_number::R, "Integer", type_enum::vector_of_int32,   "\"Total read depth for each allele\""}},
    {"ADF",      {header_number::R, "Integer", type_enum::vector_of_int32,   "\"Read depth for each allele on the forward strand\""}},
    {"ADR",      {header_number::R, "Integer", type_enum::vector_of_int32,   "\"Read depth for each allele on the reverse strand\""}},
    {"AF",       {header_number::A, "Float",   type_enum::vector_of_float32, "\"Allele frequency for each ALT allele in the same order as listed\""}},
    {"AN",       {               1, "Integer", type_enum::int32,             "\"Total number of alleles in called genotypes\""}},
    {"BQ",       {               1, "Float",   type_enum::float32,           "\"RMS base quality\""}},
    {"CIGAR",    {header_number::A, "String",  type_enum::vector_of_string,  "\"Cigar string describing how to align an alternate allele to the reference allele\""}},
    {"DB",       {               0, "Flag",    type_enum::flag,              "\"dbSNP membership\""}},
    {"DP",       {               1, "Integer", type_enum::int32,             "\"Combined depth across samples\""}},
    {"END",      {               1, "Integer", type_enum::int32,             "\"End position on CHROM (used with symbolic alleles; see below)\""}},
    {"H2",       {               0, "Flag",    type_enum::flag,              "\"HapMap2 membership\""}},
    {"H3",       {               0, "Flag",    type_enum::flag,              "\"HapMap3 membership\""}},
    {"MQ",       {               1, "Float",   type_enum::float32,           "\"RMS mapping quality\""}},
    {"MQ0",      {               1, "Integer", type_enum::int32,             "\"Number of MAPQ == 0 reads\""}},
    {"NS",       {               1, "Integer", type_enum::int32,             "\"Number of samples with data\""}},
    {"SB",       {               4, "Integer", type_enum::vector_of_int32,   "\"Strand bias\""}},
    {"SOMATIC",  {               0, "Flag",    type_enum::flag,              "\"Somatic mutation (for cancer genomics)\""}},
    {"VALIDATED",{               0, "Flag",    type_enum::flag,              "\"Validated by follow-up experiment\""}},
    {"1000G",    {               0, "Flag",    type_enum::flag,              "\"1000 Genomes membership\""}}
};

/*!\brief A compile-time mapping of FORMAT strings to bio::io::var::type_enum.
 *!\ingroup var
 *
 * This mapping allows accessing bio::io::var::genotype_variant_shallow and bio::io::var::genotype_variant_deep by e.g.
 * `get<"GT">(variant)` instead of `get<static_cast<size_t>(type_enum::string)>(variant)`.
 *
 * The predefined mappings are as defined in "Table 2" of the
 * [VCF specifiction](https://samtools.github.io/hts-specs/VCFv4.3.pdf), and also shown in
 * bio::io::var::reserved_formats.
 *
 * ### Customisation point
 *
 * This variable template is a customisation point, which means that you can specialise it for your own strings.
 *
 * TODO add example snippet
 */
template <ranges::small_string str>
inline constexpr meta::ignore_t format_key2type_enum{};

//!\cond
template <> inline constexpr type_enum format_key2type_enum<"AD">  = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"ADF"> = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"ADR"> = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"DP">  = type_enum::int32;
template <> inline constexpr type_enum format_key2type_enum<"EC">  = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"FT">  = type_enum::string;
template <> inline constexpr type_enum format_key2type_enum<"GL">  = type_enum::vector_of_float32;
template <> inline constexpr type_enum format_key2type_enum<"GP">  = type_enum::vector_of_float32;
template <> inline constexpr type_enum format_key2type_enum<"GQ">  = type_enum::int32;
template <> inline constexpr type_enum format_key2type_enum<"GT">  = type_enum::string;
template <> inline constexpr type_enum format_key2type_enum<"HQ">  = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"LAA"> = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"LAD"> = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"LGT"> = type_enum::vector_of_string;
template <> inline constexpr type_enum format_key2type_enum<"LPL"> = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"MQ">  = type_enum::int32;
template <> inline constexpr type_enum format_key2type_enum<"PL">  = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"PP">  = type_enum::vector_of_int32;
template <> inline constexpr type_enum format_key2type_enum<"PQ">  = type_enum::int32;
template <> inline constexpr type_enum format_key2type_enum<"PS">  = type_enum::int32;
//!\endcond

//!\brief A table of reserved FORMAT entries.
//!\ingroup var
inline ranges::dictionary<std::string_view, header::format_t> const reserved_formats
{
    {"AD",  {  header_number::R, "Integer", type_enum::vector_of_int32,   "\"Read depth for each allele\""}},
    {"ADF", {  header_number::R, "Integer", type_enum::vector_of_int32,   "\"Read depth for each allele on the forward strand\""}},
    {"ADR", {  header_number::R, "Integer", type_enum::vector_of_int32,   "\"Read depth for each allele on the reverse strand\""}},
    {"DP",  {                 1, "Integer", type_enum::int32,             "\"Read depth\""}},
    {"EC",  {  header_number::A, "Integer", type_enum::vector_of_int32,   "\"Expected alternate allele counts\""}},
    {"FT",  {                 1, "String",  type_enum::string,            "\"Filter indicating if this genotype was “called”\""}},
    {"GL",  {  header_number::G, "Float",   type_enum::vector_of_float32, "\"Genotype likelihoods\""}},
    {"GP",  {  header_number::G, "Float",   type_enum::vector_of_float32, "\"Genotype posterior probabilities\""}},
    {"GQ",  {                 1, "Integer", type_enum::int32,             "\"Conditional genotype quality\""}},
    {"GT",  {                 1, "String",  type_enum::string,            "\"Genotype\""}},
    {"HQ",  {                 2, "Integer", type_enum::vector_of_int32,   "\"Haplotype quality\""}},
    {"LAA", {header_number::dot, "Integer", type_enum::vector_of_int32,   "\"Strictly increasing, 1-based indices into ALT, indicating which alternate alleles are relevant (local) for the current sample\""}},
    {"LAD", {header_number::dot, "Integer", type_enum::vector_of_int32,   "\"Read depth for the reference and each of the local alternate alleles listed in LAA\""}},
    {"LGT", {header_number::dot, "String",  type_enum::vector_of_string,  "\"Genotype against the local alleles\""}},
    {"LPL", {header_number::dot, "Integer", type_enum::vector_of_int32,   "\"Phred-scaled genotype likelihoods rounded to the closest integer for genotypes that involve the reference and the local alternative alleles listed in LAA\""}},
    {"MQ",  {                 1, "Integer", type_enum::int32,             "\"RMS mapping quality\""}},
    {"PL",  {  header_number::G, "Integer", type_enum::vector_of_int32,   "\"Phred-scaled genotype likelihoods rounded to the closest integer\""}},
    {"PP",  {  header_number::G, "Integer", type_enum::vector_of_int32,   "\"Phred-scaled genotype posterior probabilities rounded to the closest integer\""}},
    {"PQ",  {                 1, "Integer", type_enum::int32,             "\"Phasing quality\""}},
    {"PS",  {                 1, "Integer", type_enum::int32,             "\"Phase set\""}}
};
// clang-format on

} // namespace bio::io::var

namespace bio::io::var::detail
{

//!\brief Data structure that represents the beginning of a BCF file.
struct bcf_header
{
    std::array<char, 3> magic{};         //!< The magic bytes.
    uint8_t             major_version{}; //!< The major version.
    uint8_t             minor_version{}; //!< The minor version.
    uint32_t            l_text;          //!< Length of the text field.
    std::string         text;            //!< The text (VCF) header.
};

} // namespace bio::io::var::detail
