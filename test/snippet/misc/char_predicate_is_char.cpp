#include <bio/io/misc/char_predicate.hpp>

int main()
{
    bio::io::is_char<'C'>('C'); // returns true

    constexpr auto my_check = bio::io::is_char<'C'>;
    my_check('c'); // returns false, because case is different
}
