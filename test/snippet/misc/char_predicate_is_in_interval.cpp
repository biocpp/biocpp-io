#include <bio/io/misc/char_predicate.hpp>

int main()
{
    bio::io::is_in_interval<'A', 'G'>('C'); // returns true

    constexpr auto my_check = bio::io::is_in_interval<'A', 'G'>;
    my_check('H'); // returns false
}
