#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>

#include <bio/io/record.hpp>

int main()
{
//![make_and_tie_record]
using namespace seqan3::literals;

std::string                  id   = "seq1";
std::vector<seqan3::dna5>    seq  = "ACGT"_dna5;
std::vector<seqan3::phred42> qual = "!!!!"_phred42;

/* This creates a *deep* record; it contains a copy of the above strings/vectors */
auto rec1 = bio::make_record(bio::vtag<bio::field::id, bio::field::seq, bio::field::qual>,   // identifiers
                             id, seq, qual);                                                 // values


/* This creates a *shallow* record; it contains references to the above strings/vectors */
auto rec2 = bio::tie_record(bio::vtag<bio::field::id, bio::field::seq, bio::field::qual>,    // identifiers
                            id, seq, qual);                                                  // values
//![make_and_tie_record]

static_assert(std::same_as<decltype(rec1),
//![make_and_tie_record_type_rec1]
bio::record<bio::vtag_t<bio::field::id, bio::field::seq, bio::field::qual>,
            seqan3::type_list<std::string,
                              std::vector<seqan3::dna5>,
                              std::vector<seqan3::phred42>>>
//![make_and_tie_record_type_rec1]
             >);

static_assert(std::same_as<decltype(rec2),
//![make_and_tie_record_type_rec2]
bio::record<bio::vtag_t<bio::field::id, bio::field::seq, bio::field::qual>,
            seqan3::type_list<std::string &,
                              std::vector<seqan3::dna5> &,
                              std::vector<seqan3::phred42> &>>
//![make_and_tie_record_type_rec2]
             >);
}
