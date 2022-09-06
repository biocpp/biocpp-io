#include <filesystem>

#include <bio/alphabet/fmt.hpp>
#include <bio/io/var_io/reader.hpp>

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
//![field_types_deep]
// this results in the records becoming "copyable"
bio::io::var_io::reader_options options{ .record = bio::io::var_io::record_default{} };

// read the entire file, copy all records into a vector; immediately closes file again
std::vector records = bio::io::var_io::reader{"example.vcf", options} | bio::ranges::to<std::vector>();

/* do something else */

// process the records later-on
for (auto & rec : records)
    fmt::print("{}:{}:{}:{}\n", rec.chrom, rec.pos, rec.ref, rec.alt);
//![field_types_deep]
}

{
//![field_types_expert]
bio::io::var_io::record r{.id          = std::ignore,
                          .ref         = std::ignore,
                          .alt         = std::ignore,
                          .qual        = std::ignore,
                          .filter      = std::ignore,
                          .info        = std::ignore,
                          .genotypes   = std::ignore};

bio::io::var_io::reader_options options{.record = r};

bio::io::var_io::reader reader{"example.vcf", options};

for (auto & rec : reader)
{
    fmt::print("{}:{}\n", rec.chrom, rec.pos);
    // all other members are just placeholders / empty!
}
//![field_types_expert]
}

//================= POST ==========================
    std::filesystem::remove("example.vcf");
}
