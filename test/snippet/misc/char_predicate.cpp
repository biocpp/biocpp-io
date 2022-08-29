#include <bio/alphabet/nucleotide/dna4.hpp>
#include <bio/io/misc/char_predicate.hpp>

int main()
{
    //! [is_eof]
    bio::io::is_eof(EOF); // returns true
    bio::io::is_eof('C'); // returns false
    //! [is_eof]

    //! [is_cntrl]
    bio::io::is_cntrl('\0'); // returns true.
    //! [is_cntrl]

    //! [is_print]
    bio::io::is_print(' '); // returns true.
    //! [is_print]

    //! [is_space]
    bio::io::is_space('\n'); // returns true.
    //! [is_space]

    //! [is_blank]
    bio::io::is_blank('\t'); // returns true.
    //! [is_blank]

    //! [is_graph]
    bio::io::is_graph('%'); // returns true.
    //! [is_graph]

    //! [is_punct]
    bio::io::is_punct(':'); // returns true.
    //! [is_punct]

    //! [is_alnum]
    bio::io::is_alnum('9'); // returns true.
    //! [is_alnum]

    //! [is_alpha]
    bio::io::is_alpha('z'); // returns true.
    //! [is_alpha]

    //! [is_upper]
    bio::io::is_upper('K'); // returns true.
    //! [is_upper]

    //! [is_lower]
    bio::io::is_lower('a'); // returns true.
    //! [is_lower]

    //! [is_digit]
    bio::io::is_digit('1'); // returns true.
    //! [is_digit]

    //! [is_xdigit]
    bio::io::is_xdigit('e'); // returns true.
    //! [is_xdigit]
}
