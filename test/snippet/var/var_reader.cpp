#include <filesystem>

#include <bio/alphabet/fmt.hpp>
#include <bio/io/var/reader.hpp>
#include <bio/test/tmp_directory.hpp>

#include "../../unit/format/vcf_data.hpp"

int main()
{
//================= PRE ==========================

    bio::test::tmp_directory dir{};
    std::filesystem::current_path(dir.path());

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

    std::ifstream in{"example.vcf", std::ios::binary};
    std::cin.rdbuf(in.rdbuf()); // rewire stdin

//================= SNIPPETS ======================

{
//![simple_usage_file]
bio::io::var::reader reader{"example.vcf"};

for (auto & rec : reader)
    fmt::print("{}:{}:{}:{}\n", rec.chrom, rec.pos, rec.ref, rec.alt);

//![simple_usage_file]
}

fmt::print("--\n");

{
//![simple_usage_stream]
bio::io::var::reader reader{std::cin, bio::io::vcf{}};

for (auto & rec : reader)
    fmt::print("{}:{}:{}:{}\n", rec.chrom, rec.pos, rec.ref, rec.alt);
//![simple_usage_stream]
}

fmt::print("--\n");

{
//![region]
bio::io::genomic_region reg{.chrom = "20", .beg = 17000, .end = 1230300};
bio::io::var::reader reader{"example.vcf.gz", bio::io::var::reader_options{.region = reg}};

// this will only print 3 records instead of 5
for (auto & rec : reader)
    fmt::print("{}:{}:{}:{}\n", rec.chrom, rec.pos, rec.ref, rec.alt);
//![region]
}

fmt::print("--\n");

{
//![views]
bio::io::var::reader reader{"example.vcf"};

auto min_qual = [](auto & rec) { return rec.qual > 23; };

for (auto & rec : reader | std::views::filter(min_qual) | std::views::take(5))
    fmt::print("{}:{}:{}:{}\n", rec.chrom, rec.pos, rec.ref, rec.alt);
//![views]
}

//================= POST ==========================
    std::filesystem::remove("example.vcf");
    std::filesystem::remove("example.vcf.gz");
    std::filesystem::remove("example.vcf.gz.tbi");
}
