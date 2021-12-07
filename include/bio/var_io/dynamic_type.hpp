// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the "dynamic typing" and auxiliaries for variant IO.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <seqan3/core/debug_stream/debug_stream_type.hpp>

#include <bio/detail/concept.hpp>
#include <bio/misc.hpp>

namespace bio::var_io
{

//!\brief Enumerator to ease "dynamic typing" in variant IO.
//!\ingroup var_io
enum class dynamic_type_id
{
    char8,
    int32,
    float32,
    string,
    vector_of_char8,
    vector_of_int8,
    vector_of_int16,
    vector_of_int32,
    vector_of_float32,
    vector_of_string,
    flag
};

/*!\brief Variant to handle "dynamic typing" in variant IO.
 * \ingroup var_io
 * \details
 *
 * TODO
 */
template <ownership own = ownership::shallow>
using dynamic_type =
  std::variant<char,
               int32_t,
               float,
               std::conditional_t<own == ownership::shallow, std::string_view, std::string>,
               std::vector<char>,
               std::vector<int8_t>,
               std::vector<int16_t>,
               std::vector<int32_t>,
               std::vector<float>,
               std::vector<std::conditional_t<own == ownership::shallow, std::string_view, std::string>>,
               bool>;

/*!\brief Variant to handle "dynamic typing" in variant IO.
 * \ingroup var_io
 * \details
 *
 * This type is similar to bio::var_io::dynamic_type except that it encodes a range of the respective types.
 * It is used for the genotype field in VCF/BCF files.
 *
 * It does not contain an entry for bio::var_io::dynamic_type_id::flag, because flags cannot appear in
 * the genotype field.
 */
template <ownership own = ownership::shallow>
using dynamic_vector_type =
  std::variant<std::vector<char>,
               std::vector<int32_t>,
               std::vector<float>,
               std::vector<std::conditional_t<own == ownership::shallow, std::string_view, std::string>>,
               std::vector<std::vector<char>>,
               std::vector<std::vector<int8_t>>,
               std::vector<std::vector<int16_t>>,
               std::vector<std::vector<int32_t>>,
               std::vector<std::vector<float>>,
               std::vector<std::vector<std::conditional_t<own == ownership::shallow, std::string_view, std::string>>>
               /* no flag here */>;

} // namespace bio::var_io

namespace seqan3
{

//!\brief TODO implement me properly
template <typename char_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, bio::var_io::dynamic_type_id const & id)
{
    // TODO print nice string
    s << (size_t)id;
    return s;
}

} // namespace seqan3

namespace bio::detail
{
//!\brief Auxilliary concept that encompasses bio::var_io::dynamic_type.
template <typename t>
concept is_dynamic_type = one_of<t, var_io::dynamic_type<ownership::shallow>, var_io::dynamic_type<ownership::deep>>;

//!\brief Auxilliary concept that encompasses bio::var_io::dynamic_vector_type.
template <typename t>
concept is_dynamic_vector_type =
  one_of<t, var_io::dynamic_vector_type<ownership::shallow>, var_io::dynamic_vector_type<ownership::deep>>;

/*!\brief Initialise an object of dynamic type to a given ID.
 * \tparam     t        Type of the output
 * \param[in]  id       The ID.
 * \param[out] output   The object being initialised.
 */
template <typename t>
    requires is_dynamic_type<t> || is_dynamic_vector_type<t>
inline void init_dynamic_type(var_io::dynamic_type_id const id, t & output)
{
    switch (id)
    {
        case var_io::dynamic_type_id::char8:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::char8);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::int32:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::int32);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::float32:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::float32);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::string:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::string);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::vector_of_char8:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::vector_of_char8);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::vector_of_int8:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::vector_of_int8);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::vector_of_int16:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::vector_of_int16);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::vector_of_int32:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::vector_of_int32);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::vector_of_float32:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::vector_of_float32);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::vector_of_string:
            {
                constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::vector_of_string);
                output.template emplace<id>();
                return;
            }
        case var_io::dynamic_type_id::flag:
            {
                if constexpr (is_dynamic_vector_type<t>)
                {
                    throw std::logic_error{"bio::var_io::dynamic_vector_type cannot be initialised to flag state."};
                }
                else
                {
                    constexpr size_t id = static_cast<size_t>(var_io::dynamic_type_id::flag);
                    output.template emplace<id>();
                }
                return;
            }
    }
}

} // namespace bio::detail
