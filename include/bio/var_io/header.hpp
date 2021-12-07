// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::var_io::header and various auxiliary classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <map>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

#include <bio/detail/charconv.hpp>
#include <bio/detail/views_eager_split.hpp>
#include <bio/exception.hpp>
#include <bio/misc.hpp>
#include <bio/var_io/dynamic_type.hpp>

namespace bio::var_io
{
/*!\brief Scoped (but weakly typed) enum for "Number" special values in bio::var_io::header INFO fields.
 * \ingroup var_io
 * \details
 *
 * ### Example
 *
 * TODO externalise
 * ```cpp
 * int32_t number = 3;
 * number = bio::header_number::A;
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
 * \ingroup var_io
 * \details
 *
 * TODO
 */
class header
{
public:
    //!\brief The dictionary for the non-standard fields in a header entry.
    using other_fields_t = std::map<std::string, std::string, std::ranges::less>;

    /*!\name Header entry types
     * \{
     */
    //!\brief Type of the contig field header line.
    struct contig_t
    {
        std::string    id;             //!< The ID.
        int64_t        length = -1;    //!< Length of the contig (-1 if absent).
        other_fields_t other_fields{}; //!< Other entries.
        int32_t        idx = -1;       //!< The numeric ID.

        //!\brief Defaulted three-way comparisons.
        auto operator<=>(contig_t const &) const = default;
    };

    //!\brief Type of a INFO field header line.
    struct info_t
    {
        std::string     id;             //!< The ID.
        int32_t         number{};       //!< Number of values, see also bio::var_io::header_number.
        dynamic_type_id type{};         //!< Type of the field.
        std::string     description{};  //!< Description.
        other_fields_t  other_fields{}; //!< Other entries.
        int32_t         idx = -1;       //!< The numeric ID.

        //!\brief Defaulted three-way comparisons.
        auto operator<=>(info_t const &) const = default;
    };

    //!\brief Type of a FILTER field header line.
    struct filter_t
    {
        std::string    id;             //!< The ID.
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
    header() { init(); }
    header(header const &) = default;             //!< Defaulted.
    header(header &&)      = default;             //!< Defaulted.
    ~header()              = default;             //!< Defaulted.
    header & operator=(header const &) = default; //!< Defaulted.
    header & operator=(header &&) = default;      //!< Defaulted.

    //!\brief Construct from a header given as plaintext.
    explicit header(std::string_view plaintext_header)
    {
        init();

        if (plaintext_header.ends_with("\r\n"))
            plaintext_header = plaintext_header.substr(0, plaintext_header.size() - 2);
        else if (plaintext_header.ends_with("\n"))
            plaintext_header = plaintext_header.substr(0, plaintext_header.size() - 1);

        for (std::string_view const line : plaintext_header | detail::eager_split('\n'))
            parse_line(line);
    }
    //!\}

    //!\brief Defaulted three-way comparison.
    auto operator<=>(header const &) const = default;

    /*!\name Header entries
     * \brief You can directly edit these member variables.
     * \{
     */
    std::string              file_format = "VCFv4.3"; //!< The file format version.
    std::vector<filter_t>    filters;                 //!< Header lines describing FILTER fields.
    std::vector<info_t>      infos;                   //!< Header lines describing INFO fields.
    std::vector<format_t>    formats;                 //!< Header lines describing FORMAT fields.
    std::vector<contig_t>    contigs;                 //!< Header lines describing contigs.
    std::vector<std::string> other_lines;             //!< Any other lines in the header.
    std::vector<std::string> column_labels;           //!< Standard column labels and sample names.
    //!\}

    /*!\name Hash maps
     * \{
     */

    //!\brief Hash map of a string-view to a position (in vector).
    using string_to_pos_t = std::unordered_map<std::string_view, size_t>;
    //!\brief Hash map of IDX to position (in vector).
    using idx_to_pos_t    = std::unordered_map<int32_t, size_t>;
    //!\brief Hash map of a string to a an IDX (this hash-map owns the strings).
    using string_to_idx_t = std::unordered_map<std::string, int32_t>;

    //!\brief ID-string to position in #filters.
    string_to_pos_t const & string_to_filter_pos() const { return string_to_filter_pos_; }
    //!\brief IDX to position in #filters.
    idx_to_pos_t const &    idx_to_filter_pos() const { return idx_to_filter_pos_; }
    //!\brief ID-string to position in #infos.
    string_to_pos_t const & string_to_info_pos() const { return string_to_info_pos_; }
    //!\brief IDX to position in #infos.
    idx_to_pos_t const &    idx_to_info_pos() const { return idx_to_info_pos_; }
    //!\brief ID-string to position in #formats.
    string_to_pos_t const & string_to_format_pos() const { return string_to_format_pos_; }
    //!\brief IDX to position in #formats.
    idx_to_pos_t const &    idx_to_format_pos() const { return idx_to_format_pos_; }
    //!\brief ID-string to position in #contigs.
    string_to_pos_t const & string_to_contig_pos() const { return string_to_contig_pos_; }
    //!\brief IDX to position in #contigs.
    idx_to_pos_t const &    idx_to_contig_pos() const { return idx_to_contig_pos_; }

    //!\brief Global string to IDX mapping (filter, info, format).
    string_to_idx_t const & string_to_idx() const { return string_to_idx_; }
    //!\brief Global string to IDX mapping (contig).
    string_to_idx_t const & contig_string_to_idx() const { return contig_string_to_idx_; }
    //!\}

    /*!\name Update, reset and inspect
     * \{
     */
    /*!\brief Add missing IDX fields to header entries and ensure that everything has proper hash-entries.
     *
     * \details
     *
     * Does the following things:
     *
     *   1. ensure that "PASS" filter entry is present as first filter entry.
     *   2. assign a valid IDX value to all header entries that currently have IDX of -1.
     *   3. generate correct hash-table mapping for all header entries.
     *
     * It does not:
     *
     *   * remove obsolete hash table entries; call #reset_hash() before to achieve this.
     *   * change existing IDX values of header entries (that are not -1).
     *   * re-use obsolete IDX values; call #reset_idx() before to achieve this.
     *
     * If you only add new entries to the header or change members other than `.id` and `.idx`, it is sufficient to
     * call this function.
     *
     * If you delete header entries, change their order, or change the `.id` / `.idx` members, you need to call
     * #reset_hash() before calling #add_missing(). If you want to start with new IDX numbering, calls #reset_idx().
     */
    void add_missing()
    {
        bool has_pass = false;
        for (size_t count = 0; count < filters.size(); ++count)
        {
            if (filters[count].id == "PASS")
            {
                has_pass                      = true;
                filters[count].idx            = 0;
                string_to_idx_["PASS"]        = 0;
                string_to_filter_pos_["PASS"] = count;
                idx_to_filter_pos_[0]         = count;
            }
            else
            {
                add_idx_and_hash_entries<entry_kind::filter>(count);
            }
        }

        if (!has_pass)
        {
            filters.push_back({"PASS", "\"All filters passed\"", {}, 0});
            string_to_idx_["PASS"]        = 0;
            string_to_filter_pos_["PASS"] = filters.size() - 1;
            idx_to_filter_pos_[0]         = filters.size() - 1;
        }

        for (size_t count = 0; count < infos.size(); ++count)
            add_idx_and_hash_entries<entry_kind::info>(count);

        for (size_t count = 0; count < formats.size(); ++count)
            add_idx_and_hash_entries<entry_kind::format>(count);

        for (size_t count = 0; count < contigs.size(); ++count)
            add_idx_and_hash_entries<entry_kind::contig>(count);
    }

    //!\brief Clear the IDX values from all header entries (sets them to -1); implicitly calls #reset_hash().
    void reset_idx()
    {
        max_contig_idx = -1;
        max_other_idx  = 0;

        for (filter_t & filter : filters)
        {
            if (filter.id == "PASS")
                filter.idx = 0;
            else
                filter.idx = -1;
        }

        for (info_t & info : infos)
            info.idx = -1;

        for (format_t & format : formats)
            format.idx = -1;

        for (contig_t & contig : contigs)
            contig.idx = -1;

        reset_hash();
    }

    //!\brief Clear all hash-maps.
    void reset_hash()
    {
        string_to_filter_pos_.clear();
        idx_to_filter_pos_.clear();
        string_to_info_pos_.clear();
        idx_to_info_pos_.clear();
        string_to_format_pos_.clear();
        idx_to_format_pos_.clear();
        string_to_contig_pos_.clear();
        idx_to_contig_pos_.clear();
        string_to_idx_.clear();
        contig_string_to_idx_.clear();
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
     *
     * This does **not** find:
     *  * Any inconsistencies in the hash tables; if in doubt, just call #reset_hash() followed by #add_missing().
     */
    void validate()
    {
        // TODO
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
    //!\brief An enum denoting the different kinds of header entries.
    enum class entry_kind
    {
        contig,
        filter,
        format,
        info
    };

    //!\brief Whether the first line was read successfully.
    bool file_format_read = false;

    //!\brief Add implicit PASS filter.
    void init()
    {
        filters.push_back({"PASS", "\"All filters passed\"", {}, 0});
        string_to_idx_["PASS"]        = 0;
        string_to_filter_pos_["PASS"] = 0;
        idx_to_filter_pos_[0]         = 0;
    }

    /*!\brief Set the correct IDX on an entry and updates the corresponding hash tables.
     * \tparam k           Pass #filters, #infos, #formats or #contigs.
     * \param[in] entry_no The position to update in the entries.
     *
     */
    template <entry_kind k>
    void add_idx_and_hash_entries(size_t const entry_no)
    {
        auto & entries = [this]() -> auto &
        {
            if constexpr (k == entry_kind::contig)
                return contigs;
            else if constexpr (k == entry_kind::filter)
                return filters;
            else if constexpr (k == entry_kind::format)
                return formats;
            else
                return infos;
        }
        ();

        string_to_idx_t & string_to_idx = k == entry_kind::contig ? contig_string_to_idx_ : string_to_idx_;
        int32_t &         max_idx       = k == entry_kind::contig ? max_contig_idx : max_other_idx;

        string_to_pos_t & string_to_pos = k == entry_kind::contig   ? string_to_contig_pos_
                                          : k == entry_kind::filter ? string_to_filter_pos_
                                          : k == entry_kind::format ? string_to_format_pos_
                                                                    : string_to_info_pos_;
        idx_to_pos_t &    idx_to_pos    = k == entry_kind::contig   ? idx_to_contig_pos_
                                          : k == entry_kind::filter ? idx_to_filter_pos_
                                          : k == entry_kind::format ? idx_to_format_pos_
                                                                    : idx_to_info_pos_;
        auto &            entry         = entries[entry_no];

        std::string_view stable_string_ref{};

        if (auto it = string_to_idx.find(entry.id); it != string_to_idx.end())
        {
            stable_string_ref = it->first;
            int32_t idx       = it->second;

            if (entry.idx == -1)
                entry.idx = idx;
            else if (idx != entry.idx)
                throw std::runtime_error{"Mismatching IDX values in header entry and hash-table. Call reset_hash()."};
        }
        else
        {
            if (entry.idx == -1)
                entry.idx = ++max_idx;

            auto [it2, insert_successful] = string_to_idx.emplace(entry.id, entry.idx);
            assert(insert_successful);
            stable_string_ref = it2->first;
        }

        /* string_to_pos has views as keys, so we cannot pass entry.id here, because
         * subsequent growing of e.g. std::vector<info_t> will invalidate that string.
         * This is because std::string has SSO which means string might be on stack.
         *
         * The strings in string_to_idx are stable with regards to growth of elements.
         */
        string_to_pos[stable_string_ref] = entry_no;
        idx_to_pos[entry.idx]            = entry_no;
    }

    /*!\name Advanced data fields
     * \brief You don't have to set these manually when creating a bio::var_io::header.
     * \{
     */
    string_to_pos_t string_to_filter_pos_; //!< ID-string to position in #filters.
    idx_to_pos_t    idx_to_filter_pos_;    //!< IDX to position in #filters.
    string_to_pos_t string_to_info_pos_;   //!< ID-string to position in #infos.
    idx_to_pos_t    idx_to_info_pos_;      //!< IDX to position in #infos.
    string_to_pos_t string_to_format_pos_; //!< ID-string to position in #formats.
    idx_to_pos_t    idx_to_format_pos_;    //!< IDX to position in #formats.
    string_to_pos_t string_to_contig_pos_; //!< ID-string to position in #contigs.
    idx_to_pos_t    idx_to_contig_pos_;    //!< IDX to position in #contigs.

    // TODO possibly store strings for these dictionaries in extra concatenated_container for cache-efficiency
    string_to_idx_t string_to_idx_;        //!< Global string to IDX mapping (filter, info, format).
    string_to_idx_t contig_string_to_idx_; //!< Global string to IDX mapping (contig).

    int32_t max_other_idx  = 0;  //!< The highest IDX value in use (defaults to 0, because PASS is used).
    int32_t max_contig_idx = -1; //!< The highest contig IDX value in use (defaults to -1, because none is used).
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
        for (auto const & filter : filters)
        {
            (raw_data += "##FILTER=<ID=") += filter.id;
            (raw_data += ",Description=") += quote_wrap(static_cast<std::string>(filter.description));

            for (auto [key, value] : filter.other_fields)
                (((raw_data += ",") += key) += "=") += value;

            if (with_idx)
                (raw_data += ",IDX=") += std::to_string(filter.idx); // TODO replace with std::to_chars

            raw_data += ">\n";
        }

        // TODO: think about if/which other_fields-value to quote_wrap
        //  infos
        for (auto const & info : infos)
        {
            (raw_data += "##INFO=<ID=") += info.id;
            (raw_data += ",Number=") += header_number::to_string(info.number);
            (raw_data += ",Type=") += unparse_type(info.type);
            (raw_data += ",Description=") += quote_wrap(static_cast<std::string>(info.description));

            for (auto [key, value] : info.other_fields)
                (((raw_data += ",") += key) += "=") += value;

            if (with_idx)
                (raw_data += ",IDX=") += std::to_string(info.idx); // TODO replace with std::to_chars

            raw_data += ">\n";
        }

        // formats
        for (auto const & format : formats)
        {
            (raw_data += "##FORMAT=<ID=") += format.id;
            (raw_data += ",Number=") += header_number::to_string(format.number);
            (raw_data += ",Type=") += unparse_type(format.type);
            (raw_data += ",Description=") += quote_wrap(static_cast<std::string>(format.description));

            for (auto [key, value] : format.other_fields)
                (((raw_data += ",") += key) += "=") += value;

            if (with_idx)
                (raw_data += ",IDX=") += std::to_string(format.idx); // TODO replace with std::to_chars

            raw_data += ">\n";
        }

        // contigs
        for (auto const & contig : contigs)
        {
            (raw_data += "##contig=<ID=") += contig.id;
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

    //!\brief Turn bio::dynamic_type_id into string.
    static std::string unparse_type(dynamic_type_id const id)
    {
        // TODO replace with string_view
        switch (id)
        {
            case dynamic_type_id::int32:
            case dynamic_type_id::vector_of_int32:
                return "Integer";
            case dynamic_type_id::float32:
            case dynamic_type_id::vector_of_float32:
                return "Float";
            case dynamic_type_id::char8:
            case dynamic_type_id::vector_of_char8:
                return "Character";
            case dynamic_type_id::string:
            case dynamic_type_id::vector_of_string:
                return "String";
            case dynamic_type_id::flag:
                return "Flag";
            default:
                throw format_error{"Illegal type in INFO or FILTER header line."};
        }
        return "";
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
        else
            new_entry.id = id.mapped();

        /* Number */
        auto number = new_entry.other_fields.extract("Number");
        if (number.empty())
            throw format_error{"INFO or FORMAT line does not contain Number field."};
        else
            new_entry.number = parse_number(number.mapped());

        /* Type */
        auto type = new_entry.other_fields.extract("Type");
        if (type.empty())
            throw format_error{"INFO or FORMAT line does not contain Type field."};
        else
            new_entry.type = parse_type(type.mapped(), new_entry.number);

        /* Description */
        auto description = new_entry.other_fields.extract("Description");
        if (description.empty())
            throw format_error{"INFO or FORMAT line does not contain Description field."};
        else
            new_entry.description = description.mapped();

        /* IDX */
        auto idx = new_entry.other_fields.extract("IDX");
        if (!idx.empty())
            detail::string_to_number(idx.mapped(), new_entry.idx);
        max_other_idx = std::max(max_other_idx, new_entry.idx);

        if (is_info)
        {
            if (string_to_info_pos_.contains(new_entry.id))
                throw format_error{std::string{"Duplicate INFO ID \""} + std::string{new_entry.id} + "\" in HEADER."};

            infos.push_back(std::move(new_entry));
            add_idx_and_hash_entries<entry_kind::info>(infos.size() - 1);
        }
        else
        {
            if (string_to_format_pos_.contains(new_entry.id))
                throw format_error{std::string{"Duplicate FORMAT ID \""} + std::string{new_entry.id} + "\" in HEADER."};

            formats.push_back(std::move(new_entry));
            add_idx_and_hash_entries<entry_kind::format>(formats.size() - 1);
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
        else
            new_entry.id = id.mapped();

        /* Description */
        auto description = new_entry.other_fields.extract("Description");
        if (description.empty())
            throw format_error{"FILTER line does not contain Description field."};
        else
            new_entry.description = description.mapped();

        /* IDX */
        auto idx = new_entry.other_fields.extract("IDX");
        if (!idx.empty())
            detail::string_to_number(idx.mapped(), new_entry.idx);
        max_other_idx = std::max(max_other_idx, new_entry.idx);

        // PASS line was added by us before and is now swapped with user-provided
        if (filters.size() > 0 && filters.front().id == "PASS" && new_entry.id == "PASS")
        {
            std::swap(filters[0], new_entry);
        }
        else
        {
            if (string_to_filter_pos_.contains(new_entry.id))
                throw format_error{std::string{"Duplicate FILTER ID \""} + std::string{new_entry.id} + "\" in HEADER."};

            filters.push_back(std::move(new_entry));
        }

        add_idx_and_hash_entries<entry_kind::filter>(filters.size() - 1);
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
        else
            new_entry.id = id.mapped();

        /* Length */
        auto length = new_entry.other_fields.extract("length");
        if (!length.empty())
            detail::string_to_number(length.mapped(), new_entry.length);

        /* IDX */
        auto idx = new_entry.other_fields.extract("IDX");
        if (idx.empty())
            new_entry.idx = ++max_contig_idx;
        else
            detail::string_to_number(idx.mapped(), new_entry.idx);
        max_contig_idx = std::max(max_contig_idx, new_entry.idx);

        if (string_to_contig_pos_.contains(new_entry.id))
            throw format_error{std::string{"Duplicate CONTIG ID \""} + std::string{new_entry.id} + "\" in HEADER."};

        contigs.push_back(std::move(new_entry));

        add_idx_and_hash_entries<entry_kind::contig>(contigs.size() - 1);
    }

    //!\brief Parse the line with column labels / sample names.
    void parse_column_labels_line(std::string_view const l)
    {
        for (std::string_view field : l | detail::eager_split('\t'))
            column_labels.push_back(static_cast<std::string>(field));
    }

    //!\brief Return a substring from the argument that does not contain enclosing angular brackets.
    static inline std::string_view strip_angular_brackets(std::string_view const in)
    {
        if (in.size() < 2 || in.front() != '<' || in.back() != '>')
            throw format_error{"Structured line does not contain \"<\" and \">\" at right places."};
        return in.substr(1, in.size() - 2);
    }

    //!\brief Turn a string into a bio::var_io::header_number.
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
                    detail::string_to_number(in, ret);
                    return ret;
                }
        }
        return header_number::dot;
    }

    /*!\brief Turn a string into bio::var_io::dynamic_type_id.
     * \param[in] in        The input string.
     * \param[in] number    The accompanying number value from the header.
     * \return The dynamic type id.
     */
    static dynamic_type_id parse_type(std::string_view const in, int32_t const number)
    {
        dynamic_type_id ret{};

        if (in == "Flag")
        {
            ret = dynamic_type_id::flag;
            if (number != 0)
                throw format_error{std::string{"Flags must always have number 0 in header."}};
            return ret;
        }

        if (number == 0)
            throw format_error{std::string{"Only flags may have number 0 in header."}};

        if (in == "Integer")
        {
            if (number == 1)
                ret = dynamic_type_id::int32;
            else
                ret = dynamic_type_id::vector_of_int32;
        }
        else if (in == "Float")
        {
            if (number == 1)
                ret = dynamic_type_id::float32;
            else
                ret = dynamic_type_id::vector_of_float32;
        }
        else if (in == "Character")
        {
            if (number == 1)
                ret = dynamic_type_id::char8;
            else
                ret = dynamic_type_id::vector_of_char8;
        }
        else if (in == "String")
        {
            if (number == 1)
                ret = dynamic_type_id::string;
            else
                ret = dynamic_type_id::vector_of_string;
        }
        else
        {
            throw format_error{std::string{"Cannot convert the following string to a type identifier: "} +
                               std::string{in}};
        }

        return ret;
    }

    //!\brief Turns a string containing value-pairs into a dictionary.
    static other_fields_t to_dictionary(std::string_view const value_pairs)
    {
        other_fields_t ret;

        for (std::string_view const pair : value_pairs | detail::eager_split(',', true))
        {
            auto pair_split = pair | detail::eager_split('=');
            auto it1        = pair_split.begin();
            auto it2        = std::ranges::next(it1);
            auto it3        = std::ranges::next(it2); // TODO whats going on here?

            if (it1 == std::default_sentinel || it2 == std::default_sentinel) //|| it3 != std::default_sentinel)
            {
                throw format_error{std::string{"Could not parse the following string into a dictionary: "} +
                                   std::string{pair}};
            }

            ret[static_cast<std::string>(*it1)] = static_cast<std::string>(*it2);
        }

        return ret;
    }
    //!\}
};

// clang-format off
//!\brief A table of reserved INFO entries.
inline std::unordered_map<std::string_view, header::info_t> const reserved_infos =
{
    {"AA",       {"AA",                   1, dynamic_type_id::string,            "\"Ancestral allele\""}},
    {"AC",       {"AC",    header_number::A, dynamic_type_id::vector_of_int32,   "\"Allele count in genotypes, for each ALT allele, in the same order as listed\""}},
    {"AD",       {"AD",    header_number::R, dynamic_type_id::vector_of_int32,   "\"Total read depth for each allele\""}},
    {"ADF",      {"ADF",   header_number::R, dynamic_type_id::vector_of_int32,   "\"Read depth for each allele on the forward strand\""}},
    {"ADR",      {"ADR",   header_number::R, dynamic_type_id::vector_of_int32,   "\"Read depth for each allele on the reverse strand\""}},
    {"AF",       {"AF",    header_number::A, dynamic_type_id::vector_of_float32, "\"Allele frequency for each ALT allele in the same order as listed\""}},
    {"AN",       {"AN",                   1, dynamic_type_id::int32,             "\"Total number of alleles in called genotypes\""}},
    {"BQ",       {"BQ",                   1, dynamic_type_id::float32,           "\"RMS base quality\""}},
    {"CIGAR",    {"CIGAR", header_number::A, dynamic_type_id::vector_of_string,  "\"Cigar string describing how to align an alternate allele to the reference allele\""}},
    {"DB",       {"DB",                   0, dynamic_type_id::flag,              "\"dbSNP membership\""}},
    {"DP",       {"DP",                   1, dynamic_type_id::int32,             "\"Combined depth across samples\""}},
    {"END",      {"END",                  1, dynamic_type_id::int32,             "\"End position on CHROM (used with symbolic alleles; see below)\""}},
    {"H2",       {"H2",                   0, dynamic_type_id::flag,              "\"HapMap2 membership\""}},
    {"H3",       {"H3",                   0, dynamic_type_id::flag,              "\"HapMap3 membership\""}},
    {"MQ",       {"MQ",                   1, dynamic_type_id::float32,           "\"RMS mapping quality\""}},
    {"MQ0",      {"MQ0",                  1, dynamic_type_id::int32,             "\"Number of MAPQ == 0 reads\""}},
    {"NS",       {"NS",                   1, dynamic_type_id::int32,             "\"Number of samples with data\""}},
    {"SB",       {"SB",                   4, dynamic_type_id::vector_of_int32,   "\"Strand bias\""}},
    {"SOMATIC",  {"SOMATIC",              0, dynamic_type_id::flag,              "\"Somatic mutation (for cancer genomics)\""}},
    {"VALIDATED",{"VALIDATED",            0, dynamic_type_id::flag,              "\"Validated by follow-up experiment\""}},
    {"1000G",    {"1000G",                0, dynamic_type_id::flag,              "\"1000 Genomes membership\""}}
};
// clang-format on

// clang-format off
//!\brief A table of reserved FORMAT entries.
inline std::unordered_map<std::string_view, header::format_t> const reserved_formats =
{
    {"AD",  {"AD",  header_number::R, dynamic_type_id::vector_of_int32,   "\"Read depth for each allele\""}},
    {"ADF", {"ADF", header_number::R, dynamic_type_id::vector_of_int32,   "\"Read depth for each allele on the forward strand\""}},
    {"ADR", {"ADR", header_number::R, dynamic_type_id::vector_of_int32,   "\"Read depth for each allele on the reverse strand\""}},
    {"DP",  {"DP",                 1, dynamic_type_id::int32,             "\"Read depth\""}},
    {"EC",  {"EC",  header_number::A, dynamic_type_id::vector_of_int32,   "\"Expected alternate allele counts\""}},
    {"FT",  {"FT",                 1, dynamic_type_id::string,            "\"Filter indicating if this genotype was “called”\""}},
    {"GL",  {"GL",  header_number::G, dynamic_type_id::vector_of_float32, "\"Genotype likelihoods\""}},
    {"GP",  {"GP",  header_number::G, dynamic_type_id::vector_of_float32, "\"Genotype posterior probabilities\""}},
    {"GQ",  {"GQ",                 1, dynamic_type_id::int32,             "\"Conditional genotype quality\""}},
    {"GT",  {"GT",                 1, dynamic_type_id::string,            "\"Genotype\""}},
    {"HQ",  {"HQ",                 2, dynamic_type_id::vector_of_int32,   "\"Haplotype quality\""}},
    {"MQ",  {"MQ",                 1, dynamic_type_id::int32,             "\"RMS mapping quality\""}},
    {"PL",  {"PL",  header_number::G, dynamic_type_id::vector_of_int32,   "\"Phred-scaled genotype likelihoods rounded to the closest integer\""}},
    {"PP",  {"PP",  header_number::G, dynamic_type_id::vector_of_int32,   "\"Phred-scaled genotype posterior probabilities rounded to the closest integer\""}},
    {"PQ",  {"PQ",                 1, dynamic_type_id::int32,             "\"Phasing quality\""}},
    {"PS",  {"PS",                 1, dynamic_type_id::int32,             "\"Phase set\""}}
};
// clang-format on

//!\brief A datastructure that contains private data of variant IO records.
//!\ingroup var_io
struct record_private_data
{
    //!\privatesection
    //!\brief Pointer to the header
    header const * header_ptr = nullptr;

    //!\brief Defaulted three-way comparison.
    friend bool operator==(record_private_data const &, record_private_data const &) = default;
    // TODO pointer to bcf-record
};

} // namespace bio::var_io

namespace bio::detail
{

//!\brief Data structure that represents the beginning of a BCF style (bit-compatible).
struct bcf_header
{
    std::array<char, 3> magic{};         //!< The magic bytes.
    uint8_t             major_version{}; //!< The major version.
    uint8_t             minor_version{}; //!< The minor version.
    uint32_t            l_text;          //!< Length of the text field.
    std::string         text;            //!< The text (VCF) header.
};

} // namespace bio::detail
