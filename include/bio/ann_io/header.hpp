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

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <bio/detail/charconv.hpp>
#include <bio/detail/views_eager_split.hpp>
#include <bio/exception.hpp>
#include <bio/misc.hpp>
#include <bio/ann_io/misc.hpp>

namespace bio::ann_io
{
class header
{
public:
    std::vector<std::pair<std::string, std::string>> browser_values{};
    std::vector<std::pair<std::string, std::string>> track_values{};

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

    //!\brief Construct from a header given as plaintext.
    explicit header(std::string_view plaintext_header)
    {
        if (plaintext_header.ends_with("\r\n"))
            plaintext_header = plaintext_header.substr(0, plaintext_header.size() - 2);
        else if (plaintext_header.ends_with("\n"))
            plaintext_header = plaintext_header.substr(0, plaintext_header.size() - 1);

        for (std::string_view const line : plaintext_header | detail::eager_split('\n'))
            parse_line(line);
    }

    /*!\name Convert to plaintext ("raw") header
     * \{
     */
    //!\brief Converts the header to plaintext (includes IDX entries).
    std::string to_plaintext() const { return to_plaintext_impl(); }
    //!\}

private:
    void parse_line(std::string_view const l)
    {
        if (l.starts_with("browser"))
        {
            auto pair_split = l.substr(8) | detail::eager_split(' ');
            auto it1        = pair_split.begin();
            auto it2        = std::ranges::next(it1);
            auto it3        = std::ranges::next(it2); // TODO whats going on here?

            if (it1 == std::default_sentinel || it2 == std::default_sentinel) //|| it3 != std::default_sentinel)
            {
                throw format_error{std::string{"Could not parse the following string into a dictionary: "} +
                                   std::string{l}};
            }

            browser_values.emplace_back(static_cast<std::string>(*it1), static_cast<std::string>(*it2));
        }
        else if (l.starts_with("track"))
        {
            for (std::string_view const pair : l.substr(6) | detail::eager_split(' ', true))
            {
                auto pair_split = pair | detail::eager_split('=');
                auto it1        = pair_split.begin();
                auto it2        = std::ranges::next(it1);

                if (it1 == std::default_sentinel || it2 == std::default_sentinel) //|| it3 != std::default_sentinel)
                {
                    throw format_error{std::string{"Could not parse the following string into a dictionary: "} +
                                       std::string{pair}};
                }

                track_values.emplace_back(static_cast<std::string>(*it1), static_cast<std::string>(strip_quotes(*it2)));
            }
        }
    }

    //!\brief Return a substring from the argument that does not contain enclosing quotes (if present).
    static inline std::string_view strip_quotes(std::string_view const in)
    {
        return (in.size() < 2 || in.front() != '"' || in.back() != '"') ? in : in.substr(1, in.size() - 2);
    }

    /*!\name Functions for converting to text
     * \{
     */
    //!\brief Implementation function for creating the plaintext header.
    std::string to_plaintext_impl() const
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

        constexpr auto is_number = [](const std::string& s)
        {
            std::string::const_iterator it = s.begin();
            while (it != s.end() && std::isdigit(*it)) ++it;
            return !s.empty() && it == s.end();
        };

        /* First print out browser settings */

        for (auto const & e : browser_values)
        {
            (((((raw_data += "browser ") += e.first) += ' ') += e.second) += '\n');
        }
        raw_data += "track ";
        for (auto const & e : track_values)
        {
            ((((raw_data += e.first) += '=') += is_number(e.second) ? e.second : quote_wrap(e.second)) += ' ');
        }

        return raw_data.substr(0, raw_data.size() - 1);
    }
    //
    // //!\brief Turn bio::value_type_id into string.
    // static std::string unparse_type(value_type_id const id)
    // {
    //     // TODO replace with string_view
    //     switch (id)
    //     {
    //         case value_type_id::int8:
    //         case value_type_id::vector_of_int8:
    //         case value_type_id::int16:
    //         case value_type_id::vector_of_int16:
    //         case value_type_id::int32:
    //         case value_type_id::vector_of_int32:
    //             return "Integer";
    //         case value_type_id::float32:
    //         case value_type_id::vector_of_float32:
    //             return "Float";
    //         case value_type_id::char8:
    //         case value_type_id::vector_of_char8:
    //             return "Character";
    //         case value_type_id::string:
    //         case value_type_id::vector_of_string:
    //             return "String";
    //         case value_type_id::flag:
    //             return "Flag";
    //         default:
    //             throw format_error{"Illegal type in INFO or FILTER header line."};
    //     }
    //     return "";
    // }
    //!\}
};
}
