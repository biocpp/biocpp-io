#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <bio/seq_io/reader_options.hpp>

int main()
{

{
//![example_simple]
bio::seq_io::reader_options options
{
    .field_types  = bio::seq_io::field_types_protein,
    .truncate_ids = true
};
//![example_simple]
}

{
//![example_advanced1]
bio::seq_io::reader_options options
{
    .field_types    = bio::seq_io::field_types<bio::ownership::shallow, seqan3::dna4, seqan3::phred42>,
    .stream_options = bio::transparent_istream_options{ .threads = 1 }
};
//![example_advanced1]
}

//TODO if https://gcc.gnu.org/bugzilla/show_bug.cgi?id=101803 gets backported to GCC10 we can omit
// bio::transparent_istream_options from the above example

{
//![example_advanced2]
bio::seq_io::reader_options options
{
    .field_ids   = bio::vtag<bio::field::seq>,
    .field_types = bio::ttag<std::string>,
    .formats     = bio::ttag<bio::fasta>
};
//![example_advanced2]
}
}
