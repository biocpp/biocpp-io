#include <filesystem>

#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/quality/phred42.hpp>
#include <bio/io/seq/record.hpp>

int main()
{

{
//![make_record]
using namespace bio::alphabet::literals;

/* assume you already have existing variables */
std::string                         id   = "my id";
std::vector<bio::alphabet::dna5>    seq  = "ACT"_dna5;
std::vector<bio::alphabet::phred42> qual = "!!!"_phred42;

/* this creates a record with copies of your existing variables */
bio::io::seq::record r{id, seq, qual};

/* the template arguments are deduced from you variables */
static_assert(std::same_as<decltype(r),
                           bio::io::seq::record<std::string,
                                                std::vector<bio::alphabet::dna5>,
                                                std::vector<bio::alphabet::phred42>>>);
//![make_record]
}

{
//![tie_record]
using namespace bio::alphabet::literals;

/* assume you already have existing variables */
std::string                         id   = "my id";
std::vector<bio::alphabet::dna5>    seq  = "ACT"_dna5;
std::vector<bio::alphabet::phred42> qual = "!!!"_phred42;

/* this creates a record of references to your existing variables */
auto r = bio::io::seq::tie_record(id, seq, qual);

/* the template arguments are deduced as references to your variables */
static_assert(std::same_as<decltype(r),
                           bio::io::seq::record<std::string &,
                                                std::vector<bio::alphabet::dna5> &,
                                                std::vector<bio::alphabet::phred42> &>>);
//![tie_record]
}
}
