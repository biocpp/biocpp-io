#include <bio/alphabet/nucleotide/dna4.hpp>
#include <bio/alphabet/quality/phred42.hpp>
#include <bio/io/seq_io/reader_options.hpp>

int main()
{

{
//![example_simple]
bio::io::seq_io::reader_options options
{
    .record         = bio::io::seq_io::record_protein_shallow{},
    .stream_options = bio::io::transparent_istream_options{ .threads = 1 },
    .truncate_ids   = true
};
//![example_simple]
}

// TODO(GCC11): omit bio::io::transparent_istream_options from the above example
// and bio::io::seq_io::record below

{
//![example_advanced]
bio::io::seq_io::reader_options options
{
    .formats     = bio::meta::ttag<bio::io::fasta>,
    .record      = bio::io::seq_io::record{.seq = std::ignore, .qual = std::ignore }
};
//![example_advanced]
}
}
