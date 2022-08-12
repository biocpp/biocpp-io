#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <bio/io/seq_io/reader_options.hpp>

int main()
{

{
//![example_simple]
bio::io::seq_io::reader_options options
{
    .field_types  = bio::io::seq_io::field_types_protein,
    .truncate_ids = true
};
//![example_simple]
}

{
//![example_advanced1]
bio::io::seq_io::reader_options options
{
    .field_types    = bio::io::seq_io::field_types<bio::io::ownership::shallow, seqan3::dna4, seqan3::phred42>,
    .stream_options = bio::io::transparent_istream_options{ .threads = 1 }
};
//![example_advanced1]
}

//TODO if https://gcc.gnu.org/bugzilla/show_bug.cgi?id=101803 gets backported to GCC10 we can omit
// bio::io::transparent_istream_options from the above example

{
//![example_advanced2]
bio::io::seq_io::reader_options options
{
    .field_ids   = bio::io::vtag<bio::io::field::seq>,
    .field_types = bio::io::ttag<std::string>,
    .formats     = bio::io::ttag<bio::io::fasta>
};
//![example_advanced2]
}
}
