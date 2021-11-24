#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <bio/seq_io/reader.hpp>

#include "../../unit/seq_io/data.hpp"

int main()
{
    //================= PRE ==========================
    {
        std::ofstream os{"example.fasta", std::ios::binary};
        os << input;
    }

    std::ifstream in{"example.fasta"};
    std::cin.rdbuf(in.rdbuf()); // rewire stdin

    //================= SNIPPETS ======================
{
//![simple_usage_file]
bio::seq_io::reader reader{"example.fasta"};

for (auto & rec : reader)
{
    seqan3::debug_stream << "ID:  " << rec.id() << '\n';
    seqan3::debug_stream << "Seq: " << rec.seq() << '\n';
}
//![simple_usage_file]
}

{
//![simple_usage_stream]
bio::seq_io::reader reader{std::cin, bio::fasta{}};

for (auto & rec : reader)
{
    seqan3::debug_stream << "ID:  " << rec.id() << '\n';
    seqan3::debug_stream << "Seq: " << rec.seq() << '\n';
}
//![simple_usage_stream]
}

{
//![decomposed]
bio::seq_io::reader reader{"example.fasta"};

for (auto & [ id, seq, qual ] : reader)
{
    seqan3::debug_stream << "ID:  " << id << '\n';
    seqan3::debug_stream << "Seq: " << seq << '\n';
}
//![decomposed]
}

{
//![views]
bio::seq_io::reader reader{"example.fasta"};

auto min_length = [](auto & rec) { return rec.seq().size() > 10; };

for (auto & rec : reader | std::views::filter(min_length) | std::views::take(5))
{
    seqan3::debug_stream << "ID:  " << rec.id() << '\n';
    seqan3::debug_stream << "Seq: " << rec.seq() << '\n';
}
//![views]
}

{
//![options]
bio::seq_io::reader reader{"example.fasta",
                           bio::seq_io::reader_options{ .field_types = bio::seq_io::field_types_protein,
                                                        .truncate_ids = true }};

for (auto & rec : reader)
{
    seqan3::debug_stream << "ID:  " << rec.id() << '\n';
    seqan3::debug_stream << "Seq: " << rec.seq() << '\n';
}
//![options]
}

    //================= POST ==========================
    std::filesystem::remove("example.fasta");
}
