// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <bio/test/tmp_filename.hpp>

#include <bio/io/txt/reader.hpp>

inline constexpr std::string_view input_no_header =
  R"raw(foo bar
bax
bat baz
)raw";

inline constexpr std::string_view input_with_extraline =
  R"raw(header
foo bar
bax
bat baz
)raw";

inline constexpr std::string_view input_with_header =
  R"raw(# header 1
# header 2
foo bar
bax
bat baz
)raw";

inline std::vector<std::string_view> lines_comp = {
  "foo bar",
  "bax",
  "bat baz",
};

inline std::vector<std::vector<std::string_view>> fields_comp = {
  {"foo", "bar" },
  {"bax"      },
  {"bat", "baz"}
};

void do_compare_linewise(auto & reader)
{
    auto it = reader.begin();
    ASSERT_TRUE(it != reader.end());
    EXPECT_EQ(*it, lines_comp[0]);

    ++it;
    ASSERT_TRUE(it != reader.end());
    EXPECT_EQ(*it, lines_comp[1]);

    ++it;
    ASSERT_TRUE(it != reader.end());
    EXPECT_EQ(*it, lines_comp[2]);

    ++it;
    ASSERT_TRUE(it == reader.end());
}

TEST(reader, line_wise_stream)
{
    std::istringstream str{static_cast<std::string>(input_no_header)};

    bio::io::txt::reader reader{str};

    do_compare_linewise(reader);
}

TEST(reader, line_wise_stream_header_first_line)
{
    std::istringstream str{static_cast<std::string>(input_with_extraline)};

    bio::io::txt::reader reader{str, bio::io::txt::header_kind::first_line};

    EXPECT_EQ(reader.header(), "header");

    do_compare_linewise(reader);
}

TEST(reader, line_wise_stream_header_starts_with)
{
    std::istringstream str{static_cast<std::string>(input_with_header)};

    bio::io::txt::reader reader{str, bio::io::txt::header_kind::starts_with{'#'}};

    EXPECT_EQ(reader.header(), "# header 1\n# header 2");

    do_compare_linewise(reader);
}

TEST(reader, line_wise_file)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << input_no_header;
    }

    bio::io::txt::reader reader{filename.get_path()};

    do_compare_linewise(reader);
}

TEST(reader, line_wise_file_header_first_line)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << input_with_extraline;
    }

    bio::io::txt::reader reader{filename.get_path(), bio::io::txt::header_kind::first_line};

    EXPECT_EQ(reader.header(), "header");

    do_compare_linewise(reader);
}

TEST(reader, line_wise_file_header_starts_with)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << input_with_header;
    }

    bio::io::txt::reader reader{filename.get_path(), bio::io::txt::header_kind::starts_with{'#'}};

    EXPECT_EQ(reader.header(), "# header 1\n# header 2");

    do_compare_linewise(reader);
}

//--------------------------------- fields wise ----------------------------

void do_compare_fields(auto & reader)
{
    auto it = reader.begin();
    ASSERT_TRUE(it != reader.end());
    EXPECT_EQ(it->line, lines_comp[0]);
    ASSERT_EQ(it->fields.size(), 2ull);
    EXPECT_EQ(it->fields[0], fields_comp[0][0]);
    EXPECT_EQ(it->fields[1], fields_comp[0][1]);

    ++it;
    ASSERT_TRUE(it != reader.end());
    EXPECT_EQ(it->line, lines_comp[1]);
    ASSERT_EQ(it->fields.size(), 1ull);
    EXPECT_EQ(it->fields[0], fields_comp[1][0]);

    ++it;
    ASSERT_TRUE(it != reader.end());
    EXPECT_EQ(it->line, lines_comp[2]);
    ASSERT_EQ(it->fields.size(), 2ull);
    EXPECT_EQ(it->fields[0], fields_comp[2][0]);
    EXPECT_EQ(it->fields[1], fields_comp[2][1]);

    ++it;
    ASSERT_TRUE(it == reader.end());
}

TEST(reader, field_wise_stream)
{
    std::istringstream str{static_cast<std::string>(input_no_header)};

    bio::io::txt::reader reader{str, ' '};

    do_compare_fields(reader);
}

TEST(reader, field_wise_stream_move)
{
    std::istringstream str{static_cast<std::string>(input_no_header)};

    bio::io::txt::reader reader{std::move(str), ' '};

    do_compare_fields(reader);
}

TEST(reader, field_wise_stream_header_first_line)
{
    std::istringstream str{static_cast<std::string>(input_with_extraline)};

    bio::io::txt::reader reader{str, ' ', bio::io::txt::header_kind::first_line};

    EXPECT_EQ(reader.header(), "header");

    do_compare_fields(reader);
}

TEST(reader, field_wise_stream_header_starts_with)
{
    std::istringstream str{static_cast<std::string>(input_with_header)};

    bio::io::txt::reader reader{str, ' ', bio::io::txt::header_kind::starts_with{'#'}};

    EXPECT_EQ(reader.header(), "# header 1\n# header 2");

    do_compare_fields(reader);
}

TEST(reader, field_wise_file)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << input_no_header;
    }

    bio::io::txt::reader reader{filename.get_path(), ' '};

    do_compare_fields(reader);
}

TEST(reader, field_wise_file_header_first_line)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << input_with_extraline;
    }

    bio::io::txt::reader reader{filename.get_path(), ' ', bio::io::txt::header_kind::first_line};

    EXPECT_EQ(reader.header(), "header");

    do_compare_fields(reader);
}

TEST(reader, field_wise_file_header_starts_with)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << input_with_header;
    }

    bio::io::txt::reader reader{filename.get_path(), ' ', bio::io::txt::header_kind::starts_with{'#'}};

    EXPECT_EQ(reader.header(), "# header 1\n# header 2");

    do_compare_fields(reader);
}

//--------------------------------- empty file ----------------------------

TEST(reader, empty_file)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};
    }

    bio::io::txt::reader reader{filename.get_path()};

    EXPECT_TRUE(reader.begin() == reader.end());
}

TEST(reader, empty_file_first_line)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << "header\n";
    }

    bio::io::txt::reader reader{filename.get_path(), ' ', bio::io::txt::header_kind::first_line};

    auto it = reader.begin();
    ASSERT_TRUE(it == reader.end());
    EXPECT_EQ(reader.header(), "header");
}

TEST(reader, empty_file_starts_with)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << "# header 1\n# header 2\n";
    }

    bio::io::txt::reader reader{filename.get_path(), ' ', bio::io::txt::header_kind::starts_with{'#'}};

    auto it = reader.begin();
    EXPECT_TRUE(it == reader.end());
    EXPECT_EQ(reader.header(), "# header 1\n# header 2");
}

//--------------------------------- fancy tests ----------------------------

// The same as line_wise_file but with buffer size of 3byte so we get many overflows
// (even more than one per line)
TEST(reader, overflow)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << input_no_header;
    }

    bio::io::txt::reader reader{filename.get_path(),
                                bio::io::txt::header_kind::none,
                                bio::io::transparent_istream_options{.buffer1_size = 3}};

    do_compare_linewise(reader);
}

// no EOL character at end of line / file
TEST(reader, no_eol)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << "header";
    }

    bio::io::txt::reader reader{filename.get_path(), ' ', bio::io::txt::header_kind::first_line};

    auto it = reader.begin();
    ASSERT_TRUE(it == reader.end());
    EXPECT_EQ(reader.header(), "header");
}
