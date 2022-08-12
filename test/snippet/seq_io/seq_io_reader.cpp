#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <bio/io/seq_io/reader.hpp>

#include "../../unit/seq_io/data.hpp"


template <typename arg1_t, typename arg2_t>
std::pair<arg1_t, arg2_t> get_arg_t(std::ranges::transform_view<arg1_t, arg2_t>);

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
bio::io::seq_io::reader reader{"example.fasta"};

for (auto & rec : reader)
{
    seqan3::debug_stream << "ID:  " << rec.id() << '\n';
    seqan3::debug_stream << "Seq: " << rec.seq() << '\n';
}
//![simple_usage_file]

using X1 = decltype(get_arg_t(reader.begin()->seq()).first);
using X2 = decltype(get_arg_t(reader.begin()->seq()).second);
using Y1 = decltype(get_arg_t(reader.begin()->qual()).first);
using Y2 = decltype(get_arg_t(reader.begin()->qual()).second);

static_assert(std::same_as<decltype(reader)::record_type,
//![simple_usage_file_type]
bio::io::record<bio::io::vtag_t<bio::io::field::id, bio::io::field::seq, bio::io::field::qual>,     // identifiers of the fields
            seqan3::type_list<std::string_view,                                 // type of the ID field
                              std::ranges::transform_view<X1, X2>,              // type of the SEQ field
                              std::ranges::transform_view<Y1, Y2>>>             // type of the QUAL field
//![simple_usage_file_type]
              >);

}

{
//![simple_usage_stream]
bio::io::seq_io::reader reader{std::cin, bio::io::fasta{}};

for (auto & rec : reader)
{
    seqan3::debug_stream << "ID:  " << rec.id() << '\n';
    seqan3::debug_stream << "Seq: " << rec.seq() << '\n';
}
//![simple_usage_stream]
}

{
//![decomposed]
bio::io::seq_io::reader reader{"example.fasta"};

for (auto & [ i, s, q ] : reader)
{
    seqan3::debug_stream << "ID:  " << i << '\n';
    seqan3::debug_stream << "Seq: " << s << '\n';
}
//![decomposed]
}

{
//![views]
bio::io::seq_io::reader reader{"example.fasta"};

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
bio::io::seq_io::reader reader{"example.fasta",
                           bio::io::seq_io::reader_options{.field_types = bio::io::seq_io::field_types_protein,
                                                       .truncate_ids = true }};

for (auto & rec : reader)
{
    seqan3::debug_stream << "ID:  " << rec.id() << '\n';
    seqan3::debug_stream << "Seq: " << rec.seq() << '\n';
}
//![options]
}

{
//![options2]
using namespace seqan3::literals;

bio::io::seq_io::reader reader{"example.fasta",
                           bio::io::seq_io::reader_options{.field_types = bio::io::seq_io::field_types<bio::io::ownership::deep>}};

for (auto & rec : reader)
{
    seqan3::debug_stream << "ID:   " << rec.id() << '\n';
    seqan3::debug_stream << "Seq:  " << rec.seq() << '\n';

    rec.seq().push_back('A'_dna5);                             // ← this is not possible with shallow records (default)
    seqan3::debug_stream << "SeqM: " << rec.seq() << '\n';
}
//![options2]

static_assert(std::same_as<decltype(reader)::record_type,
//![options2_type]
bio::io::record<bio::io::vtag_t<bio::io::field::id, bio::io::field::seq, bio::io::field::qual>,         // identifiers of the fields
            seqan3::type_list<std::string,                                          // type of the ID field
                              std::vector<seqan3::dna5>,                            // type of the SEQ field
                              std::vector<seqan3::phred63>>>                        // type of the QUAL field
//![options2_type]
              >);
}

    //================= POST ==========================
    std::filesystem::remove("example.fasta");
}
