#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <bio/var_io/reader.hpp>

#include "../../unit/format/vcf_data.hpp"

int main()
{
//================= PRE ==========================
    {
        std::ofstream os{"example.vcf", std::ios::binary};
        os << example_from_spec_header;
        os << example_from_spec_records;
    }

    {
        std::ofstream os{"example.vcf.gz", std::ios::binary};
        os << example_from_spec_bgzipped;
    }

    {
        std::ofstream os{"example.vcf.gz.tbi", std::ios::binary};
        os << example_from_spec_bgzipped_tbi;
    }

    std::ifstream in{"example.vcf"};
    std::cin.rdbuf(in.rdbuf()); // rewire stdin

//================= SNIPPETS ======================

{
//![simple_usage_file]
bio::var_io::reader reader{"example.vcf"};

for (auto & rec : reader)
{
    seqan3::debug_stream << rec.chrom() << ':'
                         << rec.pos()   << ':'
                         << rec.ref()   << ':'
                         << rec.alt()   << '\n';
}
//![simple_usage_file]
}

std::cerr << "--\n";

{
//![simple_usage_stream]
bio::var_io::reader reader{std::cin, bio::vcf{}};

for (auto & rec : reader)
{
    seqan3::debug_stream << rec.chrom() << ':'
                         << rec.pos()   << ':'
                         << rec.ref()   << ':'
                         << rec.alt()   << '\n';
}
//![simple_usage_stream]
}

std::cerr << "--\n";

{
//![region]
bio::genomic_region reg{.chrom = "20", .beg = 17000, .end = 1230300};
bio::var_io::reader reader{"example.vcf.gz", bio::var_io::reader_options{.region = reg}};

// this will only print 3 records instead of 5
for (auto & rec : reader)
{
    seqan3::debug_stream << rec.chrom() << ':'
                         << rec.pos()   << ':'
                         << rec.ref()   << ':'
                         << rec.alt()   << '\n';
}
//![region]
}

std::cerr << "--\n";

{
//![views]
bio::var_io::reader reader{"example.vcf"};

auto min_qual = [](auto & rec) { return rec.qual() > 23; };

for (auto & rec : reader | std::views::filter(min_qual) | std::views::take(5))
{
    seqan3::debug_stream << rec.chrom() << ':'
                         << rec.pos()   << ':'
                         << rec.ref()   << ':'
                         << rec.alt()   << '\n';
}
//![views]
}

//================= POST ==========================
    std::filesystem::remove("example.vcf");
    std::filesystem::remove("example.vcf.gz");
    std::filesystem::remove("example.vcf.gz.tbi");
}
