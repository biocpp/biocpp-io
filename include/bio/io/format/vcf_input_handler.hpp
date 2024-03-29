// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the bio::io::format_input_handler<bio::io::vcf>.
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

#include <bio/meta/tag/vtag.hpp>
#include <bio/meta/type_list/traits.hpp>
#include <bio/ranges/concept.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/format/format_input_handler.hpp>
#include <bio/io/format/vcf.hpp>
#include <bio/io/stream/detail/fast_streambuf_iterator.hpp>
#include <bio/io/txt/reader.hpp>
#include <bio/io/var/header.hpp>
#include <bio/io/var/misc.hpp>
#include <bio/io/var/reader_options.hpp>

namespace bio::io
{

/*!\brief Format input handler for the VCF format (bio::io::vcf).
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
 * |`print_warnings` |`bool`   | `false` | Whether to print non-critical warnings to std::cerr              |
 *
 * ### Performance
 *
 * TODO after genotype redesign
 */
template <>
class format_input_handler<vcf> :
  public format_input_handler_base<format_input_handler<vcf>>,
  public var::format_handler_mixin
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
        throw format_error{"[VCF format error in record ", line, "] ", messages...};
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
    using format_fields   = decltype(var::detail::field_ids);
    //!\brief Type of the raw record.
    using raw_record_type = io::detail::tuple_record<
      format_fields,
      meta::list_traits::concat<meta::list_traits::repeat<format_fields::size - 1, std::string_view>,
                                meta::type_list<var::record_private_data>>>;

    //!\brief Type of the low-level iterator.
    using lowlevel_iterator = txt::detail::input_iterator<txt::record_kind::line_and_fields>;

    //!\brief The raw record.
    raw_record_type   raw_record;
    //!\brief The header.
    var::header       header;
    //!\brief Lowlevel stream iterator.
    lowlevel_iterator file_it;
    //!\brief Cache of the chromosome string.
    std::string       last_chrom;
    //!\brief Cache of the chromosome index.
    int32_t           last_chrom_idx = -1;
    //!\brief Cache of the number of alleles.
    size_t            n_alts         = -1;
    //!\brief Current line number in file.
    size_t            line           = 0; //TODO this should be "record_no" and not "line"

    //!\brief Read the raw record [the base class invokes this function].
    void read_raw_record()
    {
        ++line;
        ++file_it;

        n_alts = 1; // this is overwritten by the ALT parse function (unless ALT is not being parsed)

        if (size_t field_num = file_it->fields.size(); field_num < 8)
            error("Expected at least 8 fields but got ", field_num);

        get<detail::field::chrom>(raw_record)  = (*file_it).fields[0];
        get<detail::field::pos>(raw_record)    = (*file_it).fields[1];
        get<detail::field::id>(raw_record)     = (*file_it).fields[2];
        get<detail::field::ref>(raw_record)    = (*file_it).fields[3];
        get<detail::field::alt>(raw_record)    = (*file_it).fields[4];
        get<detail::field::qual>(raw_record)   = (*file_it).fields[5];
        get<detail::field::filter>(raw_record) = (*file_it).fields[6];
        get<detail::field::info>(raw_record)   = (*file_it).fields[7];

        // fields[7].end() that is guaranteed to be char*
        char const * end_qual = (*file_it).fields[7].data() + (*file_it).fields[7].size();
        // line.end() that is guaranteed to be char*
        char const * end_line = (*file_it).line.data() + (*file_it).line.size();
        // genotypes go from end of qual til end of line (possibly empty)
        get<detail::field::genotypes>(raw_record) =
          std::string_view{end_qual, static_cast<size_t>(end_line - end_qual)};
    }
    //!\}

    /*!\name Parsed record handling
     * \brief This is mostly done via the defaults in the base class.
     * \{
     */

    // implementation after class
    template <typename t>
        //!\cond REQ
        requires(var::detail::is_info_variant<t> || var::detail::is_genotype_variant<t>)
    //!\endcond
    static void init_element_value_type(var::type_enum const id, t & output);

    // implementation after class
    struct parse_element_value_type_fn;

    // implementation after class
    static size_t parse_element_value_type(var::type_enum const                id,
                                           std::string_view const              input_string,
                                           var::detail::is_info_variant auto & output);

    //!\brief Parse the CHROM field. Reading chrom as number means getting the index (not converting string to number).
    void parse_field(meta::vtag_t<detail::field::chrom> const & /**/, auto & parsed_field)
    {
        using parsed_field_t       = std::remove_cvref_t<decltype(parsed_field)>;
        std::string_view raw_field = get<detail::field::chrom>(raw_record);

        if (raw_field != last_chrom)
        {
            // contig name was not in header, insert!
            if (auto it = header.contigs.find(raw_field); it == header.contigs.end())
            {
                header.contigs.emplace_back(static_cast<std::string>(raw_field), var::header::contig_t{});
                header.idx_update();

                last_chrom_idx = get<1>(header.contigs.back()).idx;

                if (print_warnings)
                {
                    std::cerr << "[bio::io::var::warning] The contig name \"" << raw_field << "\" found on line "
                              << line << " is not present in the header.\n";
                }
            }
            else
            {
                last_chrom_idx = get<1>(*it).idx;
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
    template <ranges::back_insertable parsed_field_t>
        requires std::ranges::range<std::ranges::range_reference_t<parsed_field_t>>
    void parse_field(meta::vtag_t<detail::field::alt> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<detail::field::alt>(raw_record);

        n_alts = 0;

        if (raw_field != ".")
        {
            for (std::string_view subfield : raw_field | detail::eager_split(','))
            {
                ++n_alts; // cache the number of alts

                parsed_field.emplace_back();
                // delegate parsing of element to base
                parse_field_aux(subfield, parsed_field.back());
            }
        }
    }

    //!\brief Overload for parsing QUAL.
    void parse_field(meta::vtag_t<detail::field::qual> const & /**/, meta::arithmetic auto & parsed_field)
    {
        std::string_view raw_field = get<detail::field::qual>(raw_record);

        if (raw_field == ".")
        {
            parsed_field = var::missing_value<std::remove_cvref_t<decltype(parsed_field)>>;
        }
        else
        {
            parse_field_aux(raw_field, parsed_field); // arithmetic parsing
        }
    }

    //!\brief Overload for parsing FILTER.
    template <ranges::back_insertable parsed_field_t>
    void parse_field(meta::vtag_t<detail::field::filter> const & /**/, parsed_field_t & parsed_field)
    {
        using elem_t               = std::ranges::range_value_t<parsed_field_t>;
        std::string_view raw_field = get<detail::field::filter>(raw_record);

        if (raw_field == ".")
            return;

        for (std::string_view subfield : raw_field | detail::eager_split(';'))
        {
            int32_t idx = -1;
            // filter name was not in header, insert!
            if (auto it = header.filters.find(subfield); it == header.filters.end())
            {
                header.filters.emplace_back(static_cast<std::string>(subfield),
                                            var::header::filter_t{.description = "\"Automatically added by SeqAn3.\""});

                header.idx_update(); // update IDX and hash-tables

                idx = get<1>(header.filters.back()).idx;

                if (print_warnings)
                {
                    std::cerr << "[bio::io::var::warning] The filter name \"" << subfield << "\" found on line " << line
                              << " was not present in the header.\n";
                }
            }
            else
            {
                idx = get<1>(*it).idx;
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
            std::cerr << "[bio::io::var::warning] The INFO name \"" << info_name << "\" found on line " << line
                      << " was not present in the header.\n";
        }

        if (var::reserved_infos.contains(info_name))
        {
            header.infos.emplace_back(info_name, var::reserved_infos.at(info_name));
        }
        else
        {
            var::header::info_t info;
            info.description = "\"Automatically added by BioC++.\"";

            if (info_value.empty()) // no "=" → flag
            {
                info.type    = "Flag";
                info.type_id = var::type_enum::flag;
                info.number  = 0;
            }
            else if (info_value.find(',') != std::string_view::npos) // found comma → assume vector-of-strings
            {
                info.type    = "String";
                info.type_id = var::type_enum::vector_of_string;
                info.number  = var::header_number::dot;
            }
            else // assume string as type
            {
                info.type    = "String";
                info.type_id = var::type_enum::string;
                info.number  = 1;
            }

            // create a new header with new info and replace current one
            header.infos.emplace_back(info_name, std::move(info));
        }

        header.idx_update();
    }

    //!\brief Overload for parsing INFO.
    template <ranges::back_insertable parsed_field_t>
        requires var::detail::info_element_reader_concept<std::ranges::range_value_t<parsed_field_t>>
    void parse_field(meta::vtag_t<detail::field::info> const & /**/, parsed_field_t & parsed_field)
    {
        using key_t = std::remove_cvref_t<
          std::tuple_element_t<0, std::remove_reference_t<std::ranges::range_value_t<parsed_field_t>>>>;
        using value_t = std::remove_cvref_t<
          std::tuple_element_t<1, std::remove_reference_t<std::ranges::range_value_t<parsed_field_t>>>>;
        // TODO this function handles value_t being string or string_view but the concept currently disallows this

        std::string_view raw_field = get<detail::field::info>(raw_record);

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
            // info name was not in header, insert!
            if (auto it = header.infos.find(key); it == header.infos.end())
            {
                add_info_to_header(key, val);
                info_pos = header.infos.size() - 1;
            }
            else
            {
                info_pos = it - header.infos.begin();
            }

            var::header::info_t & info = get<1>(header.infos[info_pos]);

            if constexpr (std::same_as<key_t, int32_t>)
                parsed_key = static_cast<key_t>(info.idx);
            else // key_t is std::string or std::string_view
                parsed_key = static_cast<key_t>(key);

            /* PARSE VALUE */
            if (val.empty()) // no "=" → flag
            {
                if constexpr (var::detail::is_info_variant<value_t>)
                {
                    if (info.type_id != var::type_enum::flag || info.number != 0)
                    {
                        error("INFO field \"", key, "\" is not a flag and should come with a value -- but does not.");
                    }

                    parsed_value = true; // set variant to boolean/flag state
                }
                // no-op for std::string
            }
            else // any other type than flag
            {
                if constexpr (var::detail::is_info_variant<value_t>)
                {
                    int32_t num_val = parse_element_value_type(info.type_id, val, parsed_value);
                    if (int32_t exp_val = info.number; print_warnings && num_val != exp_val && exp_val >= 0)
                    {
                        std::cerr << "[bio::io::var::warning] Expected to find " << exp_val
                                  << " values for the INFO field " << key << " but found: " << num_val << "\n";
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
            std::cerr << "[bio::io::var::warning] The FORMAT name \"" << format_name << "\" found on line " << line
                      << " was not present in the header.\n";
        }

        if (var::reserved_formats.contains(format_name))
        {
            header.formats.emplace_back(format_name, var::reserved_formats.at(format_name));
        }
        else
        {
            var::header::format_t format;

            format.number      = 1;
            format.type        = "String";
            format.type_id     = var::type_enum::string;
            format.description = "\"Automatically added by BioC++.\"";

            // create a new header with new format and replace current one
            header.formats.emplace_back(format_name, std::move(format));
        }

        header.idx_update();
    }

    //!\brief Overload for parsing GENOTYPES.
    template <ranges::back_insertable field_t>
        //!\cond REQ
        requires var::detail::genotype_reader_concept<std::ranges::range_value_t<field_t>>
    //!\endcond
    void parse_field(meta::vtag_t<detail::field::genotypes> const & /**/, field_t & parsed_field);

    //!\brief Overload for parsing the private data.
    void parse_field(meta::vtag_t<detail::field::_private> const & /**/, var::record_private_data & parsed_field)
    {
        parsed_field.header_ptr  = &header;
        parsed_field.raw_record  = nullptr;
        parsed_field.record_core = nullptr;
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
     * The options argument is typically bio::io::var::reader_options, but any object with a subset of similarly
     * named members is also accepted. See bio::io::format_input_handler<vcf> for the supported options and defaults.
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
        header = var::header{std::move(header_string)};
    }

    //!\brief Construct with only an input stream.
    format_input_handler(std::istream & str) : format_input_handler{str, int{}} {}
    //!\}

    //!\brief Return a reference to the header contained in the input handler.
    var::header const & get_header() const { return header; }

    //!\brief This resets the stream iterator after region-seek.
    void reset_stream() { file_it = lowlevel_iterator{*stream, false}; }
};

// ----------------------------------------------------------------------------
// out-of-line definitions of some members
// ----------------------------------------------------------------------------

/*!\brief Initialise an object of dynamic type to a given ID.
 * \tparam     t        Type of the output
 * \param[in]  id       The ID.
 * \param[out] output   The object being initialised.
 */
template <typename t>
    //!\cond REQ
    requires(var::detail::is_info_variant<t> || var::detail::is_genotype_variant<t>)
//!\endcond
inline void format_input_handler<vcf>::init_element_value_type(var::type_enum const id, t & output)
{
    switch (id)
    {
        case var::type_enum::char8:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::char8);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::int8:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::int8);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::int16:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::int16);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::int32:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::int32);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::float32:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::float32);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::string:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::string);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::vector_of_int8:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::vector_of_int8);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::vector_of_int16:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::vector_of_int16);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::vector_of_int32:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::vector_of_int32);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::vector_of_float32:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::vector_of_float32);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::vector_of_string:
            {
                constexpr size_t id = static_cast<size_t>(var::type_enum::vector_of_string);
                output.template emplace<id>();
                return;
            }
        case var::type_enum::flag:
            {
                if constexpr (var::detail::is_genotype_variant<t>)
                {
                    throw unreachable_code{std::source_location::current()};
                }
                else
                {
                    constexpr size_t id = static_cast<size_t>(var::type_enum::flag);
                    output.template emplace<id>();
                }
                return;
            }
    }
}

//!\brief Visitor definition for format_input_handler<vcf>::parse_element_value_type.
struct format_input_handler<vcf>::parse_element_value_type_fn
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
    template <meta::arithmetic arith_t>
    inline size_t operator()(arith_t & output) const
    {
        if (input == missing)
            output = var::missing_value<arith_t>;
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
                parse_element_value_type_fn{s}(vec.back());
            }
        }

        return vec.size();
    }
};

/*!\brief Parse text input into a bio::io::var::info_element_value_type /
 * bio::io::var::genotype_variant. \param[in]  id           ID of the type that shall be read. \param[in]
 * input_string The string data to read from. \param[out] output       The object to store the result into. \returns The
 * number of elements stored in the output in case ID is one of the "vector_of_"-types; 1 otherwise.
 */
template <var::detail::is_info_variant output_t>
inline size_t format_input_handler<vcf>::parse_element_value_type(var::type_enum const   id,
                                                                  std::string_view const input_string,
                                                                  output_t &             output)
{
    init_element_value_type(id, output);
    return std::visit(parse_element_value_type_fn{input_string}, output);
}

//!\brief Overload for reading the GENOTYPE field.
template <ranges::back_insertable field_t>
    requires var::detail::genotype_reader_concept<std::ranges::range_value_t<field_t>>
inline void format_input_handler<vcf>::parse_field(meta::vtag_t<detail::field::genotypes> const & /**/,
                                                   field_t & parsed_field)
{
    using genotype_field_t = std::remove_cvref_t<std::ranges::range_value_t<field_t>>;

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
        // format name was not in header, insert!
        if (auto it = header.formats.find(format_name); it == header.formats.end())
        {
            add_format_to_header(format_name);
            format_pos = header.formats.size() - 1;
        }
        else
        {
            format_pos = it - header.formats.begin();
        }

        std::ranges::range_value_t<field_t> new_element;
        auto & [current_id, current_value] = new_element;

        if constexpr (std::same_as<int32_t &, std::tuple_element_t<0, genotype_field_t> &>)
            current_id = get<1>(header.formats[format_pos]).idx;
        else
            detail::string_copy(format_name, current_id);

        auto const & format = get<1>(header.formats[format_pos]);

        init_element_value_type(format.type_id, current_value);
        auto reserve = meta::overloaded(
          [&]<typename t>(ranges::concatenated_sequences<t> & seqs)
          {
              size_t n_samples       = column_number - 8;
              size_t concat_capacity = 0;

              seqs.reserve(n_samples);

              switch (format.number)
              {
                  case bio::io::var::header_number::A:
                      concat_capacity = n_samples * n_alts;
                      break;
                  case bio::io::var::header_number::R:
                      concat_capacity = n_samples * (n_alts + 1);
                      break;
                  case bio::io::var::header_number::G:
                      concat_capacity = n_samples * (var::detail::vcf_gt_formula(n_alts, n_alts) + 1);
                      break;
                  case bio::io::var::header_number::dot:
                      // assume 1 value per sample if nothing else is known
                      concat_capacity = n_samples;
                      break;
                  case 0:
                      throw unreachable_code{std::source_location::current()};
                      break;
                  default:
                      concat_capacity = n_samples * format.number;
                      break;
              }

              seqs.concat_reserve(concat_capacity);
          },
          [&]<typename t>(std::vector<t> & vec) { vec.reserve(column_number - 8); });
        std::visit(reserve, current_value);

        ++formats;

        parsed_field.push_back(std::move(new_element));
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

                auto parse_and_append = meta::overloaded(
                  [field]<typename t>(ranges::concatenated_sequences<t> & seqs)
                  {
                      seqs.push_back();
                      if (field != ".")
                      {
                          for (std::string_view const s : field | detail::eager_split(','))
                          {
                              std::ranges::range_value_t<t> back_back;
                              parse_element_value_type_fn{s}(back_back);
                              seqs.push_back_inner(back_back);
                          }
                      }
                  },
                  [field]<typename t>(std::vector<t> & seqs)
                  {
                      seqs.emplace_back();
                      parse_element_value_type_fn{field}(seqs.back());
                  });

                auto && [current_id, current_value] = parsed_field[j];

                std::visit(parse_and_append, current_value);
            }
            else // this handles trailing dropped fields
            {
                break;
            }
        }
    }
}

} // namespace bio::io
