#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <bio/var_io/writer.hpp>
#include <bio/var_io/reader.hpp>

#include "../../unit/format/bcf_data.hpp"

int main()
{
//================= SNIPPETS ======================

{
//![creation]
using namespace std::string_literals;
using namespace seqan3::literals;

// a plaintext header
std::string_view const text_header =
R"(##fileformat=VCFv4.3
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FILTER=<ID=q10,Description="Quality below 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
)";

bio::var_io::header header{text_header};

bio::var_io::writer writer{"example.vcf"};
writer.set_header(header);
//![creation]

//![simple_usage_file]
/* construction of writer as defined above */

bio::var_io::default_record rec;
rec.chrom()     = "20";
rec.pos()       = 14370;
rec.id()        = "my_var";
rec.ref()       = "ACT"_dna5;
rec.alt().push_back("A");                                           // vector of strings by default
rec.qual()      = 9.4;
rec.filter().push_back("q10");                                      // vector of strings by default

/* info is vector over bio::var_io::info_element */
rec.info().push_back({.id = "DP", .value = 30});                    // DP is single integer (see header)
rec.info().push_back({.id = "AF", .value = std::vector{0.5f}});     // AF is vector of float (see header)

/* genotypes is vector over bio::var_io::genotype_element */
rec.genotypes().push_back({ .id = "GT", .value = std::vector{"0|0"s, "1|0"s, "1/1"s}});
// value in genotype is always a vector of size == number of samples; see bio::var_io::dynamic_vector_type

writer.push_back(rec);

/* change record and repeat */
// rec.chrom() = ...
// rec.pos() = ...
// writer.push_back(rec);
//![simple_usage_file]

//![emplace_back]
/* construction of writer as defined above */

writer.emplace_back(
    "20",
    14370,
    "my_var",
    "ACT"_dna5,
    "A",
    9.4,
    "q10",
    std::vector{bio::var_io::info_element{.id = "DP", .value = 30},
                bio::var_io::info_element{.id = "AF", .value = std::vector{0.5f}}},
    std::vector{bio::var_io::genotype_element<bio::ownership::deep>{ .id = "GT", .value = std::vector{"0|0"s, "1|0"s, "1/1"s}}});
//![emplace_back]

//![emplace_back2]
/* construction of writer as defined above */
writer.emplace_back(bio::vtag<bio::field::chrom, bio::field::pos, bio::field::ref>,
                    "20", 14370, "ACT"_dna5);
//![emplace_back2]
}

{
//![options]
bio::var_io::writer writer{"example2.vcf",
                           bio::var_io::writer_options{ .windows_eol = true }};
//![options]

// add header so destructor works
std::string_view const text_header =
R"(##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
)";
writer.set_header(bio::var_io::header{text_header});
}

{
    {
        std::ofstream tmp{"example.bcf", std::ios::binary};
        tmp << example_from_spec_bcf;
    }

{
//![inout]
bio::var_io::reader{"example.bcf"} | bio::var_io::writer{"example.vcf"};
//![inout]
}

{
//![inout2]
auto pass = [] (auto & rec)
{
    return (rec.filter().empty() || (rec.filter().size() == 1 && rec.filter()[0] == "PASS"));
};

bio::var_io::reader r{"example.bcf"};
r | std::views::filter(pass) | bio::var_io::writer{"example.vcf"};
//![inout2]
//TODO collapse the above once https://gcc.gnu.org/bugzilla/show_bug.cgi?id=103904 is resolved
}

{
//![inout3]
bio::var_io::reader reader{"example.bcf"};
bio::var_io::writer writer{"example5.vcf"};

for (auto & rec : reader)
    if (rec.filter().empty() || (rec.filter().size() == 1 && rec.filter()[0] == "PASS"))
        writer.push_back(rec);
//![inout3]
}
}

//================= POST ==========================
    std::filesystem::remove("example.bcf");
    std::filesystem::remove("example.vcf");
    std::filesystem::remove("example2.vcf");
    std::filesystem::remove("example5.vcf");
}
