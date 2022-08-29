#include <iostream>
#include <bio/io/misc/char_predicate.hpp>

int main()
{
    char chr{'1'};
    constexpr auto my_cond = bio::io::is_char<'%'> || bio::io::is_digit;
    bool is_percent = my_cond(chr);
    std::cout << std::boolalpha << is_percent << '\n'; // true
}
