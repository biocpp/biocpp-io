#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <bio/ann_io/reader.hpp>

#include "../../unit/format/ann_data.hpp"

int main()
{
//================= PRE ==========================
    {
        std::ofstream os{"example.bed", std::ios::binary};
        os << full_example;
    }

    std::ifstream in{"example.bed"};
    std::cin.rdbuf(in.rdbuf()); // rewire stdin

//================= SNIPPETS ======================

{
//![simple_usage_file]
bio::ann_io::reader reader{"example.bed"};

for (auto & rec : reader)
{
    seqan3::debug_stream << rec.chrom()        << ':'
                         << rec.chromStart()   << ':'
                         << rec.chromEnd()     << '\n';
}

seqan3::debug_stream << reader.header().header_values << '\n';
//![simple_usage_file]
}

//================= POST ==========================
    std::filesystem::remove("example.bed");
}
