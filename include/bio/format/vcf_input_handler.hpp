// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::format_input_handler<bio::vcf>.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <ranges>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/views/char_strictly_to.hpp>
#include <seqan3/core/debug_stream.hpp> //TODO evaluate if there is a better solution
#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/views/to.hpp>

#include <bio/detail/magic_get.hpp>
#include <bio/format/format_input_handler.hpp>
#include <bio/format/vcf.hpp>
#include <bio/plain_io/reader.hpp>
#include <bio/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/var_io/dynamic_type.hpp>
#include <bio/var_io/header.hpp>
#include <bio/var_io/misc.hpp>
#include <bio/var_io/reader_options.hpp>

namespace bio
{

/*!\brief Format input handler for the VCF format (bio::vcf).
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
 * ### Performance
 *
 * TODO after genotype redesign
 */
template <>
class format_input_handler<vcf> : public format_input_handler_base<format_input_handler<vcf>>
{
private:
    /*!\name CRTP related entities
     * \{
     */
    //!\brief The type of the CRTP base class.
    using base_t = format_input_handler_base<format_input_handler<vcf>>;
    using base_t::parse_field;
    using base_t::parse_field_aux;
    using base_t::stream;

    //!\brief Befriend the base class to enable CRTP.
    friend base_t;
    //!\}

    //!\brief Print an error message with current line number in diagnostic.
    [[noreturn]] void error(auto const &... messages) const
    {
        std::string message = "[SeqAn3 VCF format error in line " + detail::to_string(line) + "] ";
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
    using format_fields = decltype(var_io::default_field_ids);
    //!\brief Type of the raw record.
    using raw_record_type =
      record<format_fields,
             seqan3::list_traits::concat<seqan3::list_traits::repeat<format_fields::size - 1, std::string_view>,
                                         seqan3::type_list<var_io::record_private_data>>>;

    //!\brief Type of the low-level iterator.
    using lowlevel_iterator = detail::plaintext_input_iterator<plain_io::record_kind::line_and_fields>;

    //!\brief The raw record.
    raw_record_type   raw_record;
    //!\brief The header.
    var_io::header    header;
    //!\brief Lowlevel stream iterator.
    lowlevel_iterator file_it;
    //!\brief Cache of the chromosome string.
    std::string       last_chrom;
    //!\brief Cache of the chromosome index.
    int32_t           last_chrom_idx = -1;
    //!\brief Current line number in file.
    size_t            line           = 0;

    //!\brief Read the raw record [the base class invokes this function].
    void read_raw_record()
    {
        ++line;
        ++file_it;

        if (size_t field_num = file_it->fields.size(); field_num < 8)
            error("Expected at least 8 fields but got ", field_num);

        get<field::chrom>(raw_record)  = (*file_it).fields[0];
        get<field::pos>(raw_record)    = (*file_it).fields[1];
        get<field::id>(raw_record)     = (*file_it).fields[2];
        get<field::ref>(raw_record)    = (*file_it).fields[3];
        get<field::alt>(raw_record)    = (*file_it).fields[4];
        get<field::qual>(raw_record)   = (*file_it).fields[5];
        get<field::filter>(raw_record) = (*file_it).fields[6];
        get<field::info>(raw_record)   = (*file_it).fields[7];

        // fields[7].end() that is guaranteed to be char*
        char const * end_qual             = (*file_it).fields[7].data() + (*file_it).fields[7].size();
        // line.end() that is guaranteed to be char*
        char const * end_line             = (*file_it).line.data() + (*file_it).line.size();
        // genotypes go from end of qual til end of line (possibly empty)
        get<field::genotypes>(raw_record) = std::string_view{end_qual, static_cast<size_t>(end_line - end_qual)};
    }
    //!\}

    /*!\name Parsed record handling
     * \brief This is mostly done via the defaults in the base class.
     * \{
     */

    //!\brief Visitor definition for #parse_dynamic_type.
    struct parse_dynamic_type_fn; // implementation after class

    /*!\brief Create an bio::var_io::dynamic_type from a string and a known bio::var_io::dynamic_type_id.
     * \param[in]  id           ID of the type that shall be read.
     * \param[in]  input_string The string data to read from.
     * \param[out] output       The object to store the result in.
     * \returns The number of elements stored in the output in case ID is one of the "vector_of_"-types; 1 otherwise.
     */
    static size_t parse_dynamic_type(var_io::dynamic_type_id const  id,
                                     std::string_view const         input_string,
                                     detail::is_dynamic_type auto & output); // implementation after class

    //!\brief Parse the CHROM field. Reading chrom as number means getting the index (not converting string to number).
    void parse_field(vtag_t<field::chrom> const & /**/, auto & parsed_field)
    {
        using parsed_field_t       = std::remove_cvref_t<decltype(parsed_field)>;
        std::string_view raw_field = get<field::chrom>(raw_record);

        if (raw_field != last_chrom)
        {
            if (auto it = header.string_to_contig_pos().find(raw_field);
                it == header.string_to_contig_pos().end()) // contig name was not in header, insert!
            {
                header.contigs.push_back({.id = static_cast<std::string>(raw_field)});
                header.add_missing();

                last_chrom_idx = header.contigs[header.contigs.size() - 1].idx;

                if (print_warnings)
                {
                    seqan3::debug_stream << "[bio::var_io::warning] The contig name \"" << raw_field
                                         << "\" found on line " << line << " is not present in the header.\n";
                }
            }
            else
            {
                last_chrom_idx = header.contigs[it->second].idx;
            }

            last_chrom = raw_field;
        }

        if constexpr (std::integral<parsed_field_t>)
            parsed_field = static_cast<parsed_field_t>(last_chrom_idx);
        else
            parsed_field = static_cast<parsed_field_t>(raw_field);
    }

    /* POS, ID, REF are handled correctly by default */

    //!\brief Overload for parsing ALT.
    template <detail::back_insertable parsed_field_t>
        requires std::ranges::range<std::ranges::range_reference_t<parsed_field_t>>
    void parse_field(vtag_t<field::alt> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field::alt>(raw_record);

        if (raw_field != ".")
        {
            for (std::string_view subfield : raw_field | detail::eager_split(','))
            {
                parsed_field.emplace_back();
                // delegate parsing of element to base
                parse_field_aux(subfield, parsed_field.back());
            }
        }
    }

    //!\brief Overload for parsing QUAL.
    void parse_field(vtag_t<field::qual> const & /**/, seqan3::arithmetic auto & parsed_field)
    {
        std::string_view raw_field = get<field::qual>(raw_record);

        if (raw_field == ".")
        {
            parsed_field = var_io::missing_value<std::remove_cvref_t<decltype(parsed_field)>>;
        }
        else
        {
            parse_field_aux(raw_field, parsed_field); // arithmetic parsing
        }
    }

    //!\brief Overload for parsing FILTER.
    template <detail::back_insertable parsed_field_t>
    void parse_field(vtag_t<field::filter> const & /**/, parsed_field_t & parsed_field)
    {
        using elem_t               = std::ranges::range_value_t<parsed_field_t>;
        std::string_view raw_field = get<field::filter>(raw_record);

        if (raw_field == ".")
            return;

        for (std::string_view subfield : raw_field | detail::eager_split(';'))
        {
            int32_t idx = -1;
            if (auto it = header.string_to_filter_pos().find(subfield);
                it == header.string_to_filter_pos().end()) // filter name was not in header, insert!
            {
                header.filters.push_back(
                  {.id = static_cast<std::string>(subfield), .description = "\"Automatically added by SeqAn3.\""});

                header.add_missing(); // update IDX and hash-tables

                idx = header.filters.back().idx;

                if (print_warnings)
                {
                    seqan3::debug_stream << "[bio::var_io::warning] The filter name \"" << subfield
                                         << "\" found on line " << line << " was not present in the header.\n";
                }
            }
            else
            {
                idx = header.filters[it->second].idx;
            }

            if constexpr (std::integral<elem_t>)
                parsed_field.push_back(static_cast<elem_t>(idx));
            else
                parsed_field.push_back(static_cast<elem_t>(subfield));
        }
    }

    //!\brief Helper function that adds new INFO entry to the header.
    void add_info_to_header(std::string_view const info_name, std::string_view const info_value)
    {
        if (print_warnings)
        {
            seqan3::debug_stream << "[bio::var_io::warning] The INFO name \"" << info_name << "\" found on line "
                                 << line << " was not present in the header.\n";
        }

        if (var_io::reserved_infos.contains(info_name))
        {
            header.infos.push_back(var_io::reserved_infos.at(info_name));
        }
        else
        {
            var_io::header::info_t info;
            info.id          = info_name;
            info.description = "\"Automatically added by B.I.O..\"";

            if (info_value.empty()) // no "=" → flag
            {
                info.type   = var_io::dynamic_type_id::flag;
                info.number = 0;
            }
            else if (info_value.find(',') != std::string_view::npos) // found comma → assume vector-of-strings
            {
                info.type   = var_io::dynamic_type_id::vector_of_string;
                info.number = var_io::header_number::dot;
            }
            else // assume string as type
            {
                info.type   = var_io::dynamic_type_id::string;
                info.number = 1;
            }

            // create a new header with new info and replace current one
            header.infos.push_back(std::move(info));
        }

        header.add_missing();
    }

    //!\brief Overload for parsing INFO.
    template <detail::back_insertable parsed_field_t>
        requires detail::info_element_reader_concept<std::ranges::range_reference_t<parsed_field_t>>
    void parse_field(vtag_t<field::info> const & /**/, parsed_field_t & parsed_field)
    {
        using key_t   = detail::first_elem_t<std::ranges::range_reference_t<parsed_field_t>>;
        using value_t = detail::second_elem_t<std::ranges::range_reference_t<parsed_field_t>>;
        // TODO this function handles value_t being string or string_view but the concept currently disallows this

        std::string_view raw_field = get<field::info>(raw_record);

        if (raw_field == ".")
            return;

        for (std::string_view subfield : raw_field | detail::eager_split(';'))
        {
            std::ranges::range_value_t<parsed_field_t> ret{};
            auto & [parsed_key, parsed_value] = ret;

            auto key_value_view = subfield | detail::eager_split('=');
            auto key_it         = key_value_view.begin();
            auto val_it         = std::ranges::next(key_it);
            auto post_it        = std::ranges::next(val_it);

            if (key_it == key_value_view.end() || (val_it != key_value_view.end() && post_it != key_value_view.end()))
                error("Could not read INFO fields from the following string: ", subfield);

            std::string_view key = *key_it;
            std::string_view val = val_it == key_value_view.end() ? std::string_view{} : *val_it;

            /* PARSE KEY */
            size_t info_pos = -1;
            if (auto it = header.string_to_info_pos().find(key);
                it == header.string_to_info_pos().end()) // info name was not in header, insert!
            {
                add_info_to_header(key, val);
                info_pos = header.infos.size() - 1;
            }
            else
            {
                info_pos = it->second;
            }

            if constexpr (std::same_as<key_t, int32_t>)
                parsed_key = static_cast<key_t>(header.infos[info_pos].idx);
            else // key_t is std::string or std::string_view
                parsed_key = static_cast<key_t>(key);

            /* PARSE VALUE */
            if (val.empty()) // no "=" → flag
            {
                if constexpr (detail::is_dynamic_type<value_t>)
                {
                    if (header.infos[info_pos].type != var_io::dynamic_type_id::flag ||
                        header.infos[info_pos].number != 0)
                    {
                        error("INFO field \"", key, "\" is not a flag and should come with a value -- but does not.");
                    }

                    parsed_value = true; // set variant to boolean/flag state
                }
                // no-op for std::string
            }
            else // any other type than flag
            {
                if constexpr (detail::is_dynamic_type<value_t>)
                {
                    int32_t num_val = parse_dynamic_type(header.infos[info_pos].type, val, parsed_value);
                    if (int32_t exp_val = header.infos[info_pos].number;
                        print_warnings && num_val != exp_val && exp_val >= 0)
                    {
                        seqan3::debug_stream << "[bio::var_io::warning] Expected to find " << exp_val
                                             << " values for the INFO field " << key << " but found: " << num_val
                                             << "\n";
                    }
                }
                else // val_t == std::string or std::string_view
                {
                    parsed_value = value_t{val};
                }
            }

            parsed_field.push_back(std::move(ret));
        }
    }

    //!\brief Helper function that adds new FORMAT entry to the header.
    void add_format_to_header(std::string_view const format_name)
    {
        if (print_warnings)
        {
            seqan3::debug_stream << "[bio::var_io::warning] The FORMAT name \"" << format_name << "\" found on line "
                                 << line << " was not present in the header.\n";
        }

        if (var_io::reserved_formats.contains(format_name))
        {
            header.formats.push_back(var_io::reserved_formats.at(format_name));
        }
        else
        {
            var_io::header::format_t format;

            format.id          = format_name;
            format.number      = 1;
            format.type        = var_io::dynamic_type_id::string;
            format.description = "\"Automatically added by B.I.O..\"";

            // create a new header with new format and replace current one
            header.formats.push_back(std::move(format));
        }

        header.add_missing();
    }

    //!\brief Overload for parsing bcf style GENOTYPES.
    template <detail::back_insertable field_t>
        requires detail::genotype_bcf_style_reader_concept<std::ranges::range_reference_t<field_t>>
    void parse_field(vtag_t<field::genotypes> const & /**/, field_t & parsed_field)
    {
        using genotype_field_t = std::ranges::range_reference_t<field_t>;

        size_t column_number          = file_it->fields.size();
        size_t expected_column_number = header.column_labels.size();

        if (column_number != expected_column_number)
            error("Expected ", expected_column_number, " columns in line but found ", column_number, ".");

        if (column_number <= 8) // there are no genotypes
            return;

        std::string_view format_names = file_it->fields[8];

        /* parse keys */
        size_t formats = 0;
        for (std::string_view format_name : format_names | detail::eager_split(':'))
        {
            size_t format_pos = -1;
            if (auto it = header.string_to_format_pos().find(format_name);
                it == header.string_to_format_pos().end()) // format name was not in header, insert!
            {
                add_format_to_header(format_name);
                format_pos = header.formats.size() - 1;
            }
            else
            {
                format_pos = it->second;
            }

            parsed_field.push_back({{}, {}});
            auto & [current_id, current_value] = parsed_field.back();

            if constexpr (std::same_as<int32_t, detail::first_elem_t<genotype_field_t>>)
                current_id = header.formats[format_pos].idx;
            else
                detail::string_copy(format_name, current_id);

            auto const & format = header.formats[format_pos];

            detail::init_dynamic_type(format.type, current_value);
            auto reserve = [s = column_number - 8](auto & vec) { vec.reserve(s); };
            std::visit(reserve, current_value);

            ++formats;
        }

        /* parse values/samples */
        for (size_t i = 9; i < column_number; ++i)
        {
            std::string_view sample = file_it->fields[i];

            auto fields_view = sample | detail::eager_split(':');
            auto fields_it   = std::ranges::begin(fields_view);

            for (size_t j = 0; j < formats; ++j)
            {
                std::string_view field{};

                if (fields_it != std::default_sentinel)
                {
                    field = *fields_it;
                    ++fields_it;

                    auto parse_and_append = [field](auto & variant)
                    {
                        if constexpr (std::same_as<std::ranges::range_value_t<decltype(variant)>, bool>)
                            variant.push_back(true);
                        else
                            variant.emplace_back();

                        parse_dynamic_type_fn{field}(variant.back());
                    };

                    auto & [current_id, current_value] = parsed_field[j];

                    std::visit(parse_and_append, current_value);
                }
                else // this handles trailing dropped fields
                {
                    break;
                }
            }
        }
    }

    //!\brief Overload for parsing vcf style GENOTYPES.
    void parse_field(vtag_t<field::genotypes> const & /**/,
                     detail::genotypes_vcf_style_reader_concept auto & parsed_field)
    {
        size_t column_number          = file_it->fields.size();
        size_t expected_column_number = header.column_labels.size();

        if (column_number != expected_column_number)
            error("Expected ", expected_column_number, " columns in line but found ", column_number, ".");

        if (column_number <= 8) // there are no genotypes
            return;

        auto & [parsed_format, parsed_samples] = parsed_field;

        using string_t  = std::remove_cvref_t<decltype(parsed_format[0])>;
        using variant_t = std::remove_cvref_t<decltype(parsed_samples[0][0])>;

        /* parse formats */
        std::string_view    format_names = file_it->fields[8];
        std::vector<size_t> format_map; // ATTENTION ALWAYS DYNAMICALLY ALLOCATES HERE
        for (std::string_view format_name : format_names | detail::eager_split(':'))
        {
            size_t format_pos = -1;
            if (auto it = header.string_to_format_pos().find(format_name);
                it == header.string_to_format_pos().end()) // format name was not in header, insert!
            {
                add_format_to_header(format_name);
                format_pos = header.formats.size() - 1;
            }
            else
            {
                format_pos = it->second;
            }

            parsed_format.push_back(static_cast<string_t>(format_name));
            format_map.push_back(format_pos);
        }

        /* parse samples */
        parsed_samples.resize(column_number - 9);

        size_t sample_num = 0;
        for (std::string_view sample : file_it->fields | std::views::drop(9))
        {
            auto & parsed_sample = parsed_samples[sample_num++];

            size_t field_num = 0;
            for (std::string_view field : sample | detail::eager_split(':'))
            {
                variant_t var;

                var_io::dynamic_type_id id = header.formats[format_map[field_num++]].type;

                parse_dynamic_type(id, field, var);

                parsed_sample.push_back(std::move(var));
            }
        }
    }

    //!\brief Overload for parsing the private data.
    void parse_field(vtag_t<field::_private> const & /**/, var_io::record_private_data & parsed_field)
    {
        parsed_field = var_io::record_private_data{.header_ptr = &header};
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

        std::string header_string;
        while (file_it != std::default_sentinel && file_it.peak() == '#')
        {
            ++file_it;
            ++line;
            header_string += file_it->line;
            header_string += "\n";
        }
        header = var_io::header{std::move(header_string)};
    }

    //!\brief Construct with only an input stream.
    format_input_handler(std::istream & str) : format_input_handler{str, int{}} {}
    //!\}

    //!\brief Return a reference to the header contained in the input handler.
    var_io::header const & get_header() const { return header; }
};

//!\brief Visitor definition for format_input_handler<vcf>::parse_dynamic_type.
struct format_input_handler<vcf>::parse_dynamic_type_fn
{
    //!\brief The input data.
    std::string_view                  input;
    //!\brief Missing value constant.
    static constexpr std::string_view missing = ".";

    //!\brief Parse into bool.
    constexpr size_t operator()(bool & output) const
    {
        output = true;
        return 0;
    }

    //!\brief Parse into bool.
    inline size_t operator()(std::vector<bool>::reference output) const
    {
        output = true;
        return 0;
    }

    //!\brief Parse into character.
    constexpr size_t operator()(char & output) const
    {
        assert(input.size() == 1);
        output = input[0];
        return 1;
    }

    //!\brief Parse into number.
    template <seqan3::arithmetic arith_t>
    inline size_t operator()(arith_t & output) const
    {
        if (input == missing)
            output = var_io::missing_value<arith_t>;
        else
            detail::string_to_number(input, output);
        return 1;
    }

    //!\brief Parse into std::string.
    inline size_t operator()(std::string & output) const
    {
        if (input != missing)
            output = std::string{input};
        return 1;
    }

    //!\brief Parse into std::string_view.
    inline size_t operator()(std::string_view & output) const
    {
        output = input;
        return 1;
    }

    //!\brief Parse into vector of something.
    template <typename elem_t>
    inline size_t operator()(std::vector<elem_t> & vec) const
    {
        if (input != missing)
        {
            for (std::string_view const s : input | detail::eager_split(','))
            {
                vec.emplace_back();
                parse_dynamic_type_fn{s}(vec.back());
            }
        }

        return vec.size();
    }
};

template <detail::is_dynamic_type output_t>
inline size_t format_input_handler<vcf>::parse_dynamic_type(var_io::dynamic_type_id const id,
                                                            std::string_view const        input_string,
                                                            output_t &                    output)
{
    detail::init_dynamic_type(id, output);
    return std::visit(parse_dynamic_type_fn{input_string}, output);
}

} // namespace bio
