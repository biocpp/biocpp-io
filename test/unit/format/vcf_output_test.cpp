// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <bio/format/vcf_output_handler.hpp>
#include <bio/var_io/reader.hpp>

#include "vcf_data.hpp"

template <bio::ownership own>
void field_types_bcf_style()
{
    using namespace std::literals;

    std::ostringstream ostr{};

    {
        bio::format_output_handler<bio::vcf> handler{ostr, bio::var_io::writer_options{}};

        bio::var_io::header hdr{example_from_spec_header};
        hdr.add_missing();
        handler.set_header(std::move(hdr));

        auto recs = example_records_bcf_style<own>();

        for (auto & rec : recs)
            handler.write_record(rec);
    }

    EXPECT_EQ(ostr.str(), example_from_spec_header_regenerated_no_IDX + example_from_spec_records);
}

TEST(vcf_output, bcf_style_shallow)
{
    field_types_bcf_style<bio::ownership::shallow>();
}

TEST(vcf_output, bcf_style_deep)
{
    field_types_bcf_style<bio::ownership::deep>();
}

template <bio::ownership own>
void field_types_vcf_style()
{
    using namespace std::literals;

    std::ostringstream ostr{};

    {
        bio::format_output_handler<bio::vcf> handler{ostr, bio::var_io::writer_options{}};

        bio::var_io::header hdr{example_from_spec_header};
        handler.set_header(std::move(hdr));

        auto recs = example_records_vcf_style<own>();

        for (auto & rec : recs)
            handler.write_record(rec);
    }

    EXPECT_EQ(ostr.str(), example_from_spec_header_regenerated_no_IDX + example_from_spec_records);
}

TEST(vcf_output, vcf_style_shallow)
{
    field_types_vcf_style<bio::ownership::shallow>();
}

TEST(vcf_output, vcf_style_deep)
{
    field_types_vcf_style<bio::ownership::deep>();
}
