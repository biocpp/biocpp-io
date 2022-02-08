// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <ranges>
#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/views/char_strictly_to.hpp>
#include <seqan3/utility/views/to.hpp>

#include <bio/detail/magic_get.hpp>
#include <bio/var_io/reader.hpp>

//=============================================================================
// Official example
//=============================================================================

// https://samtools.github.io/hts-specs/VCFv4.3.pdf
inline std::string const example_from_spec_records =
R"(20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4
20	1230237	.	T	.	47	PASS	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2
20	1234567	microsat1	GTC	G,GTCT	50	PASS	NS=3;DP=9;AA=G	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3
)";

inline std::string const example_from_spec_header =
R"(##fileformat=VCFv4.3
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
)";

inline std::string const example_from_spec = example_from_spec_header + example_from_spec_records;

inline std::string const example_from_spec_header_regenerated =
R"(##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
##FILTER=<ID=q10,Description="Quality below 10",IDX=7>
##FILTER=<ID=s50,Description="Less than 50% of samples have data",IDX=8>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data",IDX=1>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth",IDX=2>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency",IDX=3>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele",IDX=4>
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129",IDX=5>
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership",IDX=6>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=9>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality",IDX=10>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth",IDX=2>
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality",IDX=11>
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x,IDX=0>
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##phasing=partial
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
)";

inline std::string const example_from_spec_header_regenerated_no_IDX =
R"(##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##phasing=partial
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
)";

inline std::string const example_from_spec_bgzipped{
"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x36\x03\x8d"
"\x93\x5b\x73\xda\x3a\x14\x85\x9f\x95\x5f\xe1\x49\xe6\xbc\xf9\x18\x5d\x2c\x19"
"\xc4\x71\x67\x9c\xb8\x06\x66\x12\xc2\xc5\xa7\x7d\x16\x20\x83\x67\x7c\x8b\x25"
"\xda\x32\xd3\x1f\x5f\xd9\x86\x34\x38\x33\x6d\x5e\xb0\x25\x6b\x7d\x6b\xef\xb5"
"\xc5\xdd\x5d\x92\x66\x32\x29\xeb\x5c\x68\xff\xcb\x43\xf4\xcd\x75\xc8\xcd\x5d"
"\xbb\x19\x0a\x2d\x7d\x0c\xe1\x08\x0e\x21\x35\x7b\xaa\x3c\xd6\x5b\xe9\xe7\xa7"
"\x59\x5e\x1d\xb5\xd0\x69\x59\x2c\xea\x72\x5f\x8b\xfc\x0b\x71\x90\x39\x50\xcb"
"\x44\xd6\xb2\x30\x67\x1a\x39\x1f\x0c\x06\x4a\xbe\x0c\x5e\x77\xd5\x00\x41\x08"
"\x27\xb2\x28\x73\xa9\x16\x69\x56\xea\x7f\xe7\x0f\xf7\x33\xc2\x9c\x44\x28\x2d"
"\x0c\x60\x5b\x16\x3a\xdd\xfb\xff\xcd\x42\xe3\x6b\x67\xb2\xd8\xeb\x83\xcf\xb0"
"\x4b\xe8\x88\xb9\xb6\x50\x4a\xe6\x9b\xec\xe4\xdf\x13\x66\xe7\x3b\xea\x27\x08"
"\xb3\xed\x2e\x19\x0a\x26\xe1\xd6\x4b\x88\x37\xda\x31\x34\x4c\x12\xc6\x36\x72"
"\x83\x77\xc2\x56\x95\xdc\xa6\x52\xf9\xb7\xd3\x32\x2f\x2d\x25\xaa\x54\x16\xea"
"\xd6\xd6\xe2\x47\x69\x8a\x38\xf9\x3f\x3e\x19\xd3\xea\x20\x54\x5a\xec\xfd\x4a"
"\xd4\x3a\x15\x99\xd9\x99\xcd\xa3\xe7\xb6\x88\xf9\xda\x9e\x1f\xf3\x8d\xac\x7d"
"\x64\xc7\xa7\x4a\xfa\xb3\x42\xcb\xbd\xac\xed\x50\xaa\x6d\x9d\x56\x4d\x04\xfe"
"\x6d\x77\xc4\x2a\x13\x6b\x2d\xf2\x2a\x93\xca\xfa\x9a\xea\x83\x65\xe2\x13\xb7"
"\x9f\xde\xe2\xc2\xc5\x47\x70\x71\xa9\x45\x66\x85\xb2\xd2\x87\x6b\x79\x10\x5d"
"\xe4\x41\x27\x8f\xb2\x52\xe8\x6b\x71\x90\x65\x32\x93\x56\x54\xcb\x97\xa3\xc9"
"\xfc\xd4\x23\x04\xbd\x02\xd6\xba\x36\xad\xf7\x10\xcd\xa8\x74\x6d\x6a\xe8\x60"
"\xbd\x1e\xee\x2f\x08\x78\x29\x42\xf4\x00\xbb\xcd\x7a\xbe\xb0\x72\xd9\x9c\x52"
"\x87\xb4\xb2\xad\xcd\x31\xcd\x76\x16\xc2\xa3\x6b\xd6\x14\xff\x95\x35\x15\xd5"
"\x93\xa8\xf0\x1b\x5a\x8b\x88\x66\x8f\xf1\xe7\x55\x0b\x79\x41\xf0\x5a\xb2\x3c"
"\x8a\x2c\xd5\x27\x6b\x23\xb3\xf2\xbb\x85\x60\x5f\xa0\x68\x4f\xf0\x28\x95\xb2"
"\xf4\x41\x14\x16\x85\xff\x34\x63\x54\xe7\x31\x1e\xc4\x37\x69\xed\x2e\x63\x8c"
"\x9e\x57\x4f\x41\xdc\x22\x26\xf1\x07\x72\x6c\x2e\xba\x36\x1f\xdf\x89\x97\x1f"
"\xb9\x05\x17\xb5\x75\x6e\xa7\x4f\xf9\xd8\x5d\x5a\x49\xb1\x7b\x73\x95\xde\xe8"
"\xa7\xaf\x55\xe0\x3f\xe8\x4d\xfc\xd9\xbb\x32\x1e\xa6\xab\xe7\x27\xb0\x78\x5e"
"\x83\x59\x08\x56\x9f\x23\x10\x3c\xc6\x60\xf9\x7f\xf0\x08\xba\x94\x41\x33\x60"
"\xd0\x79\x81\x79\x60\xfe\xf3\x10\x9d\x9f\xf8\xfc\x24\x37\x18\x02\xe4\x12\x0f"
"\x82\x5a\x31\x48\x5d\x4c\x3d\x30\x01\x01\xc0\x23\xb0\x08\xd6\x6b\x30\x5f\xfb"
"\x64\x1c\x2e\x7c\xe4\x8e\x83\xc8\x87\x0e\x1d\x87\xf7\xe3\x29\x06\x93\x98\x4f"
"\x96\x3c\x5c\xf0\xe9\x12\xc0\x9f\x90\xbb\x43\x8e\x38\x45\x36\x45\x00\x75\xcb"
"\xe1\x65\x39\x40\xdc\x25\x9c\x72\xc7\x76\x5a\x37\x8f\x10\x08\x1c\x10\x1b\x17"
"\x02\xcc\xad\xf9\xed\x81\x3a\x0f\x88\xbc\xf7\xfc\x11\x37\x8c\xa1\x4d\xa1\x59"
"\xa2\xe6\x9d\x33\x6a\x13\x00\x07\xe6\x9b\x59\xb7\x64\x84\x20\x1b\xb1\xb6\x13"
"\x17\x12\x4a\x8d\xc3\xc4\x8e\x01\xf3\x5e\x7b\xc1\xad\x0f\xec\x7c\x08\x21\x36"
"\x74\x18\xf3\xc6\x41\xe0\xc7\xa6\xb1\x2b\x57\xf4\x13\x73\x8c\x38\xe3\x98\xd8"
"\xd8\x03\xd8\xb8\x62\x0e\x39\x1a\xda\x18\xe0\x01\xe6\x84\x72\xb7\x75\xc5\x04"
"\x62\xe2\xb5\x1d\x39\xc0\xf5\x7a\xb9\x91\x16\xfe\xae\x1f\xea\x72\x8f\x53\x66"
"\x33\x78\x89\xcf\x3d\xe7\xd5\x74\xc4\x8c\xd7\x99\xed\x52\x53\x7e\x9e\x6e\xeb"
"\x52\x09\x8d\x0c\xe7\xc1\xf4\x64\x7e\x63\x60\xa2\xb8\xb2\x1a\x35\x4e\x93\x57"
"\x27\x03\x42\x6d\x91\xe6\x05\x73\xe4\x71\xdc\x8d\x02\x72\x72\xf3\x0b\xb2\x85"
"\x75\xde\x6d\x06\x00\x00\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42"
"\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00", 851};

//=============================================================================
// minimal fields
//=============================================================================

inline std::string const minimal_field_rows =
R"(20	14370	.	G	.	.	.	.	.	.	.	.
20	17330	.	T	.	.	.	.	.	.	.	.
20	1110696	.	A	.	.	.	.	.	.	.	.
20	1230237	.	T	.	.	.	.	.	.	.	.
20	1234567	.	GTC	.	.	.	.	.	.	.	.
)";

//=============================================================================
// Header misses identifiers
//=============================================================================

inline std::string const incomplete_header_before =
R"(##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data",IDX=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=2>
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##phasing=partial
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
)";

inline std::string const incomplete_header_after =
R"(##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
##FILTER=<ID=q10,Description="Automatically added by SeqAn3.",IDX=9>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data",IDX=1>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples",IDX=3>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency for each ALT allele in the same order as listed",IDX=4>
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership",IDX=5>
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership",IDX=6>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele",IDX=10>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=2>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality",IDX=7>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth",IDX=3>
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality",IDX=8>
##contig=<ID=20,IDX=0>
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##phasing=partial
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
)";


//=============================================================================
// records
//=============================================================================

/* auxiliary stuff */
using tf_view = decltype(std::string_view{} | seqan3::views::char_strictly_to<seqan3::dna5>);

bool operator==(tf_view const & lhs, tf_view const & rhs)
{
    return std::ranges::equal(lhs, rhs);
}

namespace seqan3
{

template <typename char_t, typename byte_type>
//!\cond
    requires std::same_as<std::remove_cvref_t<byte_type>, bio::var_io::record_private_data>
//!\endcond
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, byte_type &&)
{
    return s;
}


template <typename char_t, typename agg_t>
    requires bio::detail::aggregate_of_two<std::remove_cvref_t<agg_t>>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, agg_t && agg)
{
    s << '[' << bio::detail::get_first(agg) << ", " << bio::detail::get_second(agg) << ']';
    return s;
}


} // namespace bio

template <bio::ownership own>
auto make_ref(std::string_view const str)
{
    if constexpr (own == bio::ownership::shallow)
        return tf_view{{str, {}}, {}};
    else
        return str | seqan3::views::char_strictly_to<seqan3::dna5> | seqan3::views::to<std::vector>;
}

template <bio::ownership own, typename int_t = int32_t>
auto example_records_default_style()
{
    using record_t = bio::record<decltype(bio::var_io::default_field_ids),
                                 decltype(bio::var_io::field_types<own>)>;

    bio::var_io::record_private_data priv{};
    constexpr int_t mv = bio::var_io::missing_value<int_t>;
    using ivec          = std::vector<int_t>;
    using ivecvec       = std::vector<std::vector<int_t>>;
    using fvec          = std::vector<float>;
    using svec          = std::conditional_t<own == bio::ownership::shallow,
                                             std::vector<std::string_view>,
                                             std::vector<std::string>>;

    // clang-format off
    std::vector<record_t> recs{
    {"20", 14370,   "rs6054257", make_ref<own>("G"),   {"A"},        29, {"PASS"}, {{"NS",(int_t)3}, {"DP", (int_t)14}, {"AF", fvec{0.5f}}, {"DB", true}, {"H2", true}        }, { {"GT", svec{"0|0", "1|0", "1/1"}}, {"GQ", ivec{48, 48, 43}}, {"DP", ivec{1, 8, 5}}, {"HQ", ivecvec{{51,51}, {51,51}, {mv,mv} }}}, priv},
    {"20", 17330,   ".",         make_ref<own>("T"),   {"A"},        3,  {"q10"},  {{"NS",(int_t)3}, {"DP", (int_t)11}, {"AF", fvec{0.017f}}                                  }, { {"GT", svec{"0|0", "0|1", "0/0"}}, {"GQ", ivec{49,  3, 41}}, {"DP", ivec{3, 5, 3}}, {"HQ", ivecvec{{58,50}, {65, 3}          }}}, priv},
    {"20", 1110696, "rs6040355", make_ref<own>("A"),   {"G","T"},    67, {"PASS"}, {{"NS",(int_t)2}, {"DP", (int_t)10}, {"AF", fvec{0.333f,0.667f}}, {"AA", "T"}, {"DB", true}}, { {"GT", svec{"1|2", "2|1", "2/2"}}, {"GQ", ivec{21,  2, 35}}, {"DP", ivec{6, 0, 4}}, {"HQ", ivecvec{{23,27}, {18, 2}          }}}, priv},
    {"20", 1230237, ".",         make_ref<own>("T"),   {},           47, {"PASS"}, {{"NS",(int_t)3}, {"DP", (int_t)13}, {"AA", "T"}                                           }, { {"GT", svec{"0|0", "0|0", "0/0"}}, {"GQ", ivec{54, 48, 61}}, {"DP", ivec{7, 4, 2}}, {"HQ", ivecvec{{56,60}, {51,51}          }}}, priv},
    {"20", 1234567, "microsat1", make_ref<own>("GTC"), {"G","GTCT"}, 50, {"PASS"}, {{"NS",(int_t)3}, {"DP", (int_t)9 }, {"AA", "G"}                                           }, { {"GT", svec{"0/1", "0/2", "1/1"}}, {"GQ", ivec{35, 17, 40}}, {"DP", ivec{4, 2, 3}}                                             }, priv},
    };
    // clang-format on

    return recs;

}

template <bio::ownership own, typename int_t = int32_t>
auto example_records_vcf_style()
{
    using record_t = bio::record<decltype(bio::var_io::default_field_ids),
                                 decltype(bio::var_io::field_types_vcf_style<own>)>;

    bio::var_io::record_private_data priv{};
    constexpr int_t mv = bio::var_io::missing_value<int_t>;
    using ivec = std::vector<int_t>;
    using fvec = std::vector<float>;

    // clang-format off
    std::vector<record_t> recs{
    {"20", 14370,   "rs6054257", make_ref<own>("G"),   {"A"},        29, {"PASS"}, {{"NS",(int_t)3}, {"DP", (int_t)14}, {"AF", fvec{0.5f}}, {"DB", true}, {"H2", true}        }, { {"GT", "GQ", "DP", "HQ"}, {{"0|0", (int_t)48,(int_t)1,ivec{51,51}}, {"1|0",(int_t)48,(int_t)8,ivec{51,51}}, {"1/1",(int_t)43,(int_t)5,ivec{mv,mv}}}}, priv},
    {"20", 17330,   ".",         make_ref<own>("T"),   {"A"},        3,  {"q10"},  {{"NS",(int_t)3}, {"DP", (int_t)11}, {"AF", fvec{0.017f}}                                  }, { {"GT", "GQ", "DP", "HQ"}, {{"0|0", (int_t)49,(int_t)3,ivec{58,50}}, {"0|1",(int_t) 3,(int_t)5,ivec{65, 3}}, {"0/0",(int_t)41,(int_t)3            }}}, priv},
    {"20", 1110696, "rs6040355", make_ref<own>("A"),   {"G","T"},    67, {"PASS"}, {{"NS",(int_t)2}, {"DP", (int_t)10}, {"AF", fvec{0.333f,0.667f}}, {"AA", "T"}, {"DB", true}}, { {"GT", "GQ", "DP", "HQ"}, {{"1|2", (int_t)21,(int_t)6,ivec{23,27}}, {"2|1",(int_t) 2,(int_t)0,ivec{18, 2}}, {"2/2",(int_t)35,(int_t)4            }}}, priv},
    {"20", 1230237, ".",         make_ref<own>("T"),   {},           47, {"PASS"}, {{"NS",(int_t)3}, {"DP", (int_t)13}, {"AA", "T"}                                           }, { {"GT", "GQ", "DP", "HQ"}, {{"0|0", (int_t)54,(int_t)7,ivec{56,60}}, {"0|0",(int_t)48,(int_t)4,ivec{51,51}}, {"0/0",(int_t)61,(int_t)2            }}}, priv},
    {"20", 1234567, "microsat1", make_ref<own>("GTC"), {"G","GTCT"}, 50, {"PASS"}, {{"NS",(int_t)3}, {"DP", (int_t)9 }, {"AA", "G"}                                           }, { {"GT", "GQ", "DP"      }, {{"0/1", (int_t)35,(int_t)4            }, {"0/2",(int_t)17,(int_t)2            }, {"1/1",(int_t)40,(int_t)3            }}}, priv}
    };
    // clang-format on

    return recs;
}

template <bio::ownership own, typename int_t = int32_t>
auto example_records_bcf_style()
{
    using record_t = bio::record<decltype(bio::var_io::default_field_ids),
                                 decltype(bio::var_io::field_types_bcf_style<own>)>;

    bio::var_io::record_private_data priv{};
    constexpr int_t mv  = bio::var_io::missing_value<int_t>;
    using ivec          = std::vector<int_t>;
    using ivecvec       = std::vector<std::vector<int_t>>;
    using fvec          = std::vector<float>;
    using svec          = std::conditional_t<own == bio::ownership::shallow,
                                             std::vector<std::string_view>,
                                             std::vector<std::string>>;

    // clang-format off
    std::vector<record_t> recs{
    {0, 14370,   "rs6054257", make_ref<own>("G"),   {"A"},        29, {0}, {{1,(int_t)3}, {2, (int_t)14}, {3, fvec{0.5f}}, {5, true}, {6, true}        }, { {9, svec{"0|0", "1|0", "1/1"}}, {10, ivec{48, 48, 43}}, {2, ivec{1, 8, 5}}, {11, ivecvec{{51,51}, {51,51}, {mv,mv} }}}, priv},
    {0, 17330,   ".",         make_ref<own>("T"),   {"A"},        3,  {7}, {{1,(int_t)3}, {2, (int_t)11}, {3, fvec{0.017f}}                            }, { {9, svec{"0|0", "0|1", "0/0"}}, {10, ivec{49,  3, 41}}, {2, ivec{3, 5, 3}}, {11, ivecvec{{58,50}, {65, 3}          }}}, priv},
    {0, 1110696, "rs6040355", make_ref<own>("A"),   {"G","T"},    67, {0}, {{1,(int_t)2}, {2, (int_t)10}, {3, fvec{0.333f,0.667f}}, {4, "T"}, {5, true}}, { {9, svec{"1|2", "2|1", "2/2"}}, {10, ivec{21,  2, 35}}, {2, ivec{6, 0, 4}}, {11, ivecvec{{23,27}, {18, 2}          }}}, priv},
    {0, 1230237, ".",         make_ref<own>("T"),   {},           47, {0}, {{1,(int_t)3}, {2, (int_t)13}, {4, "T"}                                     }, { {9, svec{"0|0", "0|0", "0/0"}}, {10, ivec{54, 48, 61}}, {2, ivec{7, 4, 2}}, {11, ivecvec{{56,60}, {51,51}          }}}, priv},
    {0, 1234567, "microsat1", make_ref<own>("GTC"), {"G","GTCT"}, 50, {0}, {{1,(int_t)3}, {2, (int_t)9 }, {4, "G"}                                     }, { {9, svec{"0/1", "0/2", "1/1"}}, {10, ivec{35, 17, 40}}, {2, ivec{4, 2, 3}}                                           }, priv},
    };
    // clang-format on

    return recs;
}

auto example_records_novariant()
{
    using namespace std::string_view_literals;

    bio::var_io::record_private_data priv{};
    constexpr int32_t mv = bio::var_io::missing_value<int32_t>;
    using ivec          = std::vector<int32_t>;
    using ivecvec       = std::vector<std::vector<int32_t>>;
    using fvec          = std::vector<float>;
    using svec          = std::vector<std::string_view>;

    auto rec0 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 14370,
                                 "rs6054257",
                                 "G",
                                 svec{"A"},
                                 29.0,
                                 svec{"PASS"},
                                 std::tuple{std::pair{"NS", 3},
                                            std::pair{"DP", 14},
                                            std::pair{"AF", fvec{0.5f}},
                                            std::pair{"DB", true},
                                            std::pair{"H2", true}},
                                 std::tuple{std::pair{"GT", svec{"0|0", "1|0", "1/1"}},
                                            std::pair{"GQ", ivec{48, 48, 43}},
                                            std::pair{"DP", ivec{1, 8, 5}},
                                            std::pair{"HQ", ivecvec{{51,51}, {51,51}, {mv,mv} }}},
                                 priv);

    auto rec1 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 17330,
                                 ".",
                                 "T",
                                 svec{"A"},
                                 3.0,
                                 svec{"q10"},
                                 std::tuple{std::pair{"NS", 3},
                                            std::pair{"DP", 11},
                                            std::pair{"AF", fvec{0.017f}}},
                                 std::tuple{std::pair{"GT", svec{"0|0", "0|1", "0/0"}},
                                            std::pair{"GQ", ivec{49,  3, 41}},
                                            std::pair{"DP", ivec{3, 5, 3}},
                                            std::pair{"HQ", ivecvec{{58,50}, {65, 3}}}},
                                 priv);

    auto rec2 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 1110696,
                                 "rs6040355",
                                 "A",
                                 svec{"G","T"},
                                 67,
                                 svec{"PASS"},
                                 std::tuple{std::pair{"NS",2},
                                            std::pair{"DP", 10},
                                            std::pair{"AF", fvec{0.333f,0.667f}},
                                            std::pair{"AA", "T"},
                                            std::pair{"DB", true}},
                                 std::tuple{std::pair{"GT", svec{"1|2", "2|1", "2/2"}},
                                            std::pair{"GQ", ivec{21,  2, 35}},
                                            std::pair{"DP", ivec{6, 0, 4}},
                                            std::pair{"HQ", ivecvec{{23,27}, {18, 2}}}},
                                 priv);

    auto rec3 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 1230237,
                                 ".",
                                 "T",
                                 svec{},
                                 47,
                                 svec{"PASS"},
                                 std::tuple{std::pair{"NS",3},
                                            std::pair{"DP", 13},
                                            std::pair{"AA", "T"} },
                                 std::tuple{std::pair{"GT", svec{"0|0", "0|0", "0/0"}},
                                            std::pair{"GQ", ivec{54, 48, 61}},
                                            std::pair{"DP", ivec{7, 4, 2}},
                                            std::pair{"HQ", ivecvec{{56,60}, {51,51}}}},
                                 priv);
    auto rec4 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 1234567,
                                 "microsat1",
                                 "GTC",
                                 svec{"G","GTCT"},
                                 50,
                                 svec{"PASS"},
                                 std::tuple{std::pair{"NS",3},
                                            std::pair{"DP", 9 },
                                            std::pair{"AA", "G"}},
                                 std::tuple{std::pair{"GT", svec{"0/1", "0/2", "1/1"}},
                                            std::pair{"GQ", ivec{35, 17, 40}},
                                            std::pair{"DP", ivec{4, 2, 3}}},
                                 priv);

    return std::tuple{rec0, rec1, rec2, rec3, rec4};
}

auto example_records_novariant_vcf_style_genotypes()
{
    using namespace std::string_view_literals;

    bio::var_io::record_private_data priv{};
    constexpr int32_t mv = bio::var_io::missing_value<int32_t>;
    using ivec          = std::vector<int32_t>;
    using fvec          = std::vector<float>;
    using svec          = std::vector<std::string_view>;

    auto rec0 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 14370,
                                 "rs6054257",
                                 "G",
                                 svec{"A"},
                                 29.0,
                                 svec{"PASS"},
                                 std::tuple{std::pair{"NS", 3},
                                            std::pair{"DP", 14},
                                            std::pair{"AF", fvec{0.5f}},
                                            std::pair{"DB", true},
                                            std::pair{"H2", true}},
                                 std::pair{svec{"GT", "GQ", "DP", "HQ"},
                                           std::tuple{std::tuple{"0|0", 48,1,ivec{51,51}},
                                                      std::tuple{"1|0",48,8,ivec{51,51}},
                                                      std::tuple{"1/1",43,5,ivec{mv,mv}}}},
                                 priv);

    auto rec1 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 17330,
                                 ".",
                                 "T",
                                 svec{"A"},
                                 3.0,
                                 svec{"q10"},
                                 std::tuple{std::pair{"NS", 3},
                                            std::pair{"DP", 11},
                                            std::pair{"AF", fvec{0.017f}}},
                                 std::pair{svec{"GT", "GQ", "DP", "HQ"},
                                           std::tuple{std::tuple{"0|0",49,3,ivec{58,50}},
                                                      std::tuple{"0|1", 3,5,ivec{65, 3}},
                                                      std::tuple{"0/0",41,3}}},
                                 priv);

    auto rec2 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 1110696,
                                 "rs6040355",
                                 "A",
                                 svec{"G","T"},
                                 67,
                                 svec{"PASS"},
                                 std::tuple{std::pair{"NS",2},
                                            std::pair{"DP", 10},
                                            std::pair{"AF", fvec{0.333f,0.667f}},
                                            std::pair{"AA", "T"},
                                            std::pair{"DB", true}},
                                 std::pair{svec{"GT", "GQ", "DP", "HQ"},
                                           std::tuple{std::tuple{"1|2",21,6,ivec{23,27}},
                                                      std::tuple{"2|1", 2,0,ivec{18, 2}},
                                                      std::tuple{"2/2",35,4}}},
                                 priv);

    auto rec3 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 1230237,
                                 ".",
                                 "T",
                                 svec{},
                                 47,
                                 svec{"PASS"},
                                 std::tuple{std::pair{"NS",3},
                                            std::pair{"DP", 13},
                                            std::pair{"AA", "T"} },
                                 std::pair{svec{"GT", "GQ", "DP", "HQ"},
                                           std::tuple{std::tuple{"0|0",54,7,ivec{56,60}},
                                                      std::tuple{"0|0",48,4,ivec{51,51}},
                                                      std::tuple{"0/0",61,2}}},
                                 priv);
    auto rec4 = bio::make_record(bio::var_io::default_field_ids,
                                 "20",
                                 1234567,
                                 "microsat1",
                                 "GTC",
                                 svec{"G","GTCT"},
                                 50,
                                 svec{"PASS"},
                                 std::tuple{std::pair{"NS",3},
                                            std::pair{"DP", 9 },
                                            std::pair{"AA", "G"}},
                                 std::pair{svec{"GT", "GQ", "DP"},
                                           std::tuple{std::tuple{"0/1",35,4},
                                                      std::tuple{"0/2",17,2},
                                                      std::tuple{"1/1",40,3}}},
                                 priv);

    return std::tuple{rec0, rec1, rec2, rec3, rec4};
}
