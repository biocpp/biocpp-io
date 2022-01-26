// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <bio/format/bcf_output_handler.hpp>
#include <bio/var_io/reader.hpp>

#include "bcf_data.hpp"
#include "vcf_data.hpp"

enum class style
{
    def,
    vcf,
    bcf
};

template <style s, bio::ownership own>
void field_types()
{
    std::ostringstream ostr{};

    {
        bio::format_output_handler<bio::bcf> handler{ostr, bio::var_io::writer_options{}};

        bio::var_io::header hdr{example_from_spec_bcf_header};
        hdr.add_missing();
        handler.set_header(std::move(hdr));

        auto recs = []()
        {
            if constexpr (s == style::def)
                return example_records_default_style<own>();
            else if constexpr (s == style::vcf)
                return example_records_vcf_style<own>();
            else
                return example_records_bcf_style<own>();
        }();

        for (auto & rec : recs)
            handler.write_record(rec);
    }

    EXPECT_EQ(ostr.str(), example_from_spec_bcf_unbgzf_our);
}

TEST(bcf_output, default_style_shallow)
{
    field_types<style::def, bio::ownership::shallow>();
}

TEST(bcf_output, default_style_deep)
{
    field_types<style::def, bio::ownership::deep>();
}

// TEST(bcf_output, vcf_style_shallow)
// {
//     field_types<style::vcf, bio::ownership::shallow>();
// }
//
// TEST(bcf_output, vcf_style_deep)
// {
//     field_types<style::vcf, bio::ownership::deep>();
// }

TEST(bcf_output, bcf_style_shallow)
{
    field_types<style::bcf, bio::ownership::shallow>();
}

TEST(bcf_output, bcf_style_deep)
{
    field_types<style::bcf, bio::ownership::deep>();
}

TEST(bcf_output, novariant)
{
    std::ostringstream ostr{};

    {
        bio::format_output_handler<bio::bcf> handler{ostr, bio::var_io::writer_options{}};

        bio::var_io::header hdr{example_from_spec_header};
        hdr.add_missing();
        handler.set_header(std::move(hdr));

        auto records = example_records_novariant(); // < records is a tuple here

        std::apply([&](auto &... recs) { (handler.write_record(recs), ...); }, records);
    }

    EXPECT_EQ(ostr.str(), example_from_spec_bcf_unbgzf_our);
}

// TEST(bcf_output, novariant_vcf_style_genotypes)
// {
//     std::ostringstream ostr{};
//
//     {
//         bio::format_output_handler<bio::bcf> handler{ostr, bio::var_io::writer_options{}};
//
//         bio::var_io::header hdr{example_from_spec_header};
//         hdr.add_missing();
//         handler.set_header(std::move(hdr));
//
//         auto records = example_records_novariant_vcf_style_genotypes(); // < records is a tuple here
//
//         std::apply([&](auto &... recs) { (handler.write_record(recs), ...); }, records);
//     }
//
//     EXPECT_EQ(ostr.str(), example_from_spec_header_regenerated_no_IDX + example_from_spec_records);
// }

