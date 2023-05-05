#include <bio/alphabet/nucleotide/dna4.hpp>
#include <bio/alphabet/quality/phred42.hpp>
#include <bio/io/seq/reader_options.hpp>

int main()
{

{
//![example_simple]
bio::io::seq::reader_options options
{
    .record         = bio::io::seq::record_protein_shallow{},
    .stream_options = { .threads = 1 },
    .truncate_ids   = true
};
//![example_simple]
}

{
//![example_advanced]
bio::io::seq::reader_options options
{
    .formats     = bio::meta::ttag<bio::io::fasta>,
    .record      = bio::io::seq::record{.seq = std::ignore, .qual = std::ignore }
};
//![example_advanced]
}
}
