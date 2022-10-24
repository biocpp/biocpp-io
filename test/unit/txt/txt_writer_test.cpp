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

#include <bio/io/txt/writer.hpp>

inline constexpr std::string_view compare =
  R"(foo bar
bax
bat baz 3.4 7
)";

TEST(writer, line_wise_push_back_stream)
{
    std::ostringstream str;

    bio::io::txt::writer writer{str};

    writer.push_back("foo bar");
    writer.push_back("bax");
    writer.push_back("bat baz 3.4 7");

    str.flush();
    EXPECT_EQ(str.str(), compare);
}

TEST(writer, line_wise_emplace_back_stream)
{
    std::ostringstream str;

    bio::io::txt::writer writer{str};

    writer.emplace_back("foo", " ", "bar");
    writer.emplace_back("bax");
    writer.emplace_back("bat", " ", "baz", " ", 3.4, " ", 7);

    str.flush();
    EXPECT_EQ(str.str(), compare);
}

TEST(writer, line_wise_push_back_file)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        bio::io::txt::writer writer{filename.get_path()};
        writer.push_back("foo bar");
        writer.push_back("bax");
        writer.push_back("bat baz 3.4 7");
    }

    std::ifstream     ifs{filename.get_path()};
    std::stringstream buffer;
    buffer << ifs.rdbuf();

    buffer.flush();
    EXPECT_EQ(buffer.str(), compare);
}

TEST(writer, line_wise_emplace_back_file)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        bio::io::txt::writer writer{filename.get_path()};
        writer.emplace_back("foo", " ", "bar");
        writer.emplace_back("bax");
        writer.emplace_back("bat", " ", "baz", " ", 3.4, " ", 7);
    }

    std::ifstream     ifs{filename.get_path()};
    std::stringstream buffer;
    buffer << ifs.rdbuf();

    buffer.flush();
    EXPECT_EQ(buffer.str(), compare);
}

//--------------------------------- fields wise ----------------------------

using rng_t = std::vector<std::string_view>;

TEST(writer, field_wise_push_back_stream)
{
    std::ostringstream str;

    bio::io::txt::writer writer{str, ' '};

    writer.push_back(rng_t{"foo", "bar"});
    writer.push_back(rng_t{"bax"});
    writer.push_back(rng_t{"bat", "baz", "3.4", "7"});

    str.flush();
    EXPECT_EQ(str.str(), compare);
}

TEST(writer, field_wise_emplace_back_stream)
{
    std::ostringstream str;

    bio::io::txt::writer writer{str, ' '};

    writer.emplace_back("foo", "bar"); // space is inserted automatically
    writer.emplace_back("bax");
    writer.emplace_back("bat", "baz", 3.4, 7); // space is inserted automatically

    str.flush();
    EXPECT_EQ(str.str(), compare);
}

TEST(writer, field_wise_push_back_file)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        bio::io::txt::writer writer{filename.get_path(), ' '};
        writer.push_back(rng_t{"foo", "bar"});
        writer.push_back(rng_t{"bax"});
        writer.push_back(rng_t{"bat", "baz", "3.4", "7"});
    }

    std::ifstream     ifs{filename.get_path()};
    std::stringstream buffer;
    buffer << ifs.rdbuf();

    buffer.flush();
    EXPECT_EQ(buffer.str(), compare);
}

TEST(writer, field_wise_emplace_back_file)
{
    bio::test::tmp_filename filename{"txt_test"};

    {
        bio::io::txt::writer writer{filename.get_path(), ' '};
        writer.emplace_back("foo", "bar"); // space is inserted automatically
        writer.emplace_back("bax");
        writer.emplace_back("bat", "baz", 3.4, 7); // space is inserted automatically
    }

    std::ifstream     ifs{filename.get_path()};
    std::stringstream buffer;
    buffer << ifs.rdbuf();

    buffer.flush();
    EXPECT_EQ(buffer.str(), compare);
}

// TODO write by record
