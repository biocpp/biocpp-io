#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/quality/phred42.hpp>

#include <bio/io/record.hpp>

int main()
{
//![make_and_tie_record]
using namespace bio::alphabet::literals;

std::string                  id   = "seq1";
std::vector<bio::alphabet::dna5>    seq  = "ACGT"_dna5;
std::vector<bio::alphabet::phred42> qual = "!!!!"_phred42;

/* This creates a *deep* record; it contains a copy of the above strings/vectors */
auto rec1 = bio::io::make_record(bio::io::vtag<bio::io::field::id, bio::io::field::seq, bio::io::field::qual>,   // identifiers
                             id, seq, qual);                                                 // values


/* This creates a *shallow* record; it contains references to the above strings/vectors */
auto rec2 = bio::io::tie_record(bio::io::vtag<bio::io::field::id, bio::io::field::seq, bio::io::field::qual>,    // identifiers
                            id, seq, qual);                                                  // values
//![make_and_tie_record]

static_assert(std::same_as<decltype(rec1),
//![make_and_tie_record_type_rec1]
bio::io::record<bio::io::vtag_t<bio::io::field::id, bio::io::field::seq, bio::io::field::qual>,
            seqan3::type_list<std::string,
                              std::vector<bio::alphabet::dna5>,
                              std::vector<bio::alphabet::phred42>>>
//![make_and_tie_record_type_rec1]
             >);

static_assert(std::same_as<decltype(rec2),
//![make_and_tie_record_type_rec2]
bio::io::record<bio::io::vtag_t<bio::io::field::id, bio::io::field::seq, bio::io::field::qual>,
            seqan3::type_list<std::string &,
                              std::vector<bio::alphabet::dna5> &,
                              std::vector<bio::alphabet::phred42> &>>
//![make_and_tie_record_type_rec2]
             >);
}
