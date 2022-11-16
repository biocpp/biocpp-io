#include <filesystem>

#include <bio/io/seq/writer.hpp>
#include <bio/io/seq/reader.hpp>

#include "../../unit/seq/data.hpp"

int main()
{
//================= SNIPPETS ======================

{
//![creation]
bio::io::seq::writer writer{"example.fasta"};
//![creation]
}

{
//![simple_usage_file]
using namespace bio::alphabet::literals;

bio::io::seq::writer writer{"example.fastq"};

bio::io::seq::record_dna rec;
rec.id   = "my id";
rec.seq  = "ACT"_dna5;
rec.qual = "!!!"_phred42;
writer.push_back(rec);

/* change record and repeat */
// rec.id  = ...
// rec.seq = ...
// writer.push_back(rec);
//![simple_usage_file]
}

{
//![emplace_back]
using namespace bio::alphabet::literals;

bio::io::seq::writer writer{"example.fasta"};

writer.emplace_back("my id", "ACT"_dna5, ""); // last field can be empty for FastA format
//![emplace_back]
}

{
//![options]
bio::io::seq::writer writer{"example.fasta",
                            bio::io::seq::writer_options{ .max_seq_line_length = 0 }};
//![options]
}

{
    {
        std::ofstream tmp{"example.fastq", std::ios::binary};
        tmp << interleaved_fastq;
    }

{
//![inout]
bio::io::seq::reader{"example.fastq"} | bio::io::seq::writer{"example.fasta"};
//![inout]
}

{
//![inout2]
auto min_size = [] (auto & rec)
{
    return rec.seq.size() > 20;
};

bio::io::seq::reader r{"example.fastq"};
r | std::views::filter(min_size) | bio::io::seq::writer{"example.fasta"};

// In GCC>=12 and Clang>=16, you can write the entire expression in one line:
// bio::io::seq::reader{"example.bcf"} | std::views::filter(min_size) | bio::io::seq::writer{"example.vcf"};
//![inout2]
//TODO(gcc12): collapse the above
}

{
//![inout3]
bio::io::seq::reader reader{"example.fastq"};
bio::io::seq::writer writer{"example.fasta"};

for (auto & rec : reader)
    if (rec.seq.size() > 20)
        writer.push_back(rec);
//![inout3]
}
}

//================= POST ==========================
    std::filesystem::remove("example.fasta");
    std::filesystem::remove("example.fastq");
}
