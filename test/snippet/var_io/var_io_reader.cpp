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
}
