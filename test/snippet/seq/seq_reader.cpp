#include <filesystem>

#include <bio/io/seq/reader.hpp>
#include <bio/alphabet/fmt.hpp>

#include "../../unit/seq/data.hpp"


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
bio::io::seq::reader reader{"example.fasta"};

for (auto & rec : reader)
{
    fmt::print("ID:  {}\n", rec.id);
    fmt::print("Seq: {}\n", rec.seq);
}
//![simple_usage_file]

using X1 = decltype(get_arg_t(reader.begin()->seq).first);
using X2 = decltype(get_arg_t(reader.begin()->seq).second);
using Y1 = decltype(get_arg_t(reader.begin()->qual).first);
using Y2 = decltype(get_arg_t(reader.begin()->qual).second);

static_assert(std::same_as<decltype(reader)::record_type,
//![simple_usage_file_type]
bio::io::seq::record<std::string_view,                    // type of the ID field
                     std::ranges::transform_view<X1, X2>, // type of the SEQ field
                     std::ranges::transform_view<Y1, Y2>> // type of the QUAL field
//![simple_usage_file_type]
              >);

}

{
//![simple_usage_stream]
bio::io::seq::reader reader{std::cin, bio::io::fasta{}};

for (auto & rec : reader)
{
    fmt::print("ID:  {}\n", rec.id);
    fmt::print("Seq: {}\n", rec.seq);
}
//![simple_usage_stream]
}

{
//![decomposed]
bio::io::seq::reader reader{"example.fasta"};

for (auto & [ i, s, q ] : reader)
{
    fmt::print("ID:  {}\n", i);
    fmt::print("Seq: {}\n", s);
}
//![decomposed]
}

{
//![views]
bio::io::seq::reader reader{"example.fasta"};

auto min_length = [](auto & rec) { return rec.seq.size() > 10; };

for (auto & rec : reader | std::views::filter(min_length) | std::views::take(5))
{
    fmt::print("ID:  {}\n", rec.id);
    fmt::print("Seq: {}\n", rec.seq);
}
//![views]
}

{
//![options]
bio::io::seq::reader reader{"example.fasta",
                            bio::io::seq::reader_options{
                                .record = bio::io::seq::record_protein_shallow{},
                                .truncate_ids = true }};

for (auto & rec : reader)
{
    fmt::print("ID:  {}\n", rec.id);
    fmt::print("Seq: {}\n", rec.seq);
}
//![options]
}

{
//![options2]
using namespace bio::alphabet::literals;

bio::io::seq::reader reader{"example.fasta",
                            bio::io::seq::reader_options{.record = bio::io::seq::record_dna_deep{}}};

for (auto & rec : reader)
{
    fmt::print("ID:   {}\n", rec.id);
    fmt::print("Seq:  {}\n", rec.seq);

    rec.seq.push_back('A'_dna5);           // ← this is not possible with shallow records (default)
    fmt::print("SeqM: {}\n", rec.seq);
}
//![options2]

static_assert(std::same_as<decltype(reader)::record_type,
//![options2_type]
bio::io::seq::record<std::string,                         // type of the ID field
                     std::vector<bio::alphabet::dna5>,    // type of the SEQ field
                     std::vector<bio::alphabet::phred42>> // type of the QUAL field
//![options2_type]
              >);


}

    //================= POST ==========================
    std::filesystem::remove("example.fasta");
}
