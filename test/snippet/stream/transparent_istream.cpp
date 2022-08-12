#include <fstream>

#include <bio/io/stream/transparent_istream.hpp>
#include <bio/io/stream/transparent_ostream.hpp>

int main()
{
    { // create my_file.txt and my_cin_file.txt so the code below works
        std::ofstream os{"my_file.txt"};
        os << "Foo";
        std::ofstream os2{"my_cin_file.txt"};
        os2 << "Foo";
    }

    { // create my_file.txt.gz so the code below works
        bio::transparent_ostream os{"my_file.txt.gz"};
        os << "Foo";
    }

    // reset std::cin buffer to be usable in this snippet
    std::ifstream in{"my_cin_file.txt"};
    std::cin.rdbuf(in.rdbuf());

    {
    //![construction]
    std::string buffer;

    bio::transparent_istream s1{"my_file.txt"};      // behaves like std::ifstream
    s1 >> buffer;

    bio::transparent_istream s2{"my_file.txt.gz"};   // file is transparently decompressed
    s2 >> buffer;

    bio::transparent_istream s3{std::cin};        // wrap standard input
    s3 >> buffer;
    //![construction]
    }

    {
    //![decompression]
    std::string buffer;

    bio::transparent_istream s1{"my_file.txt.gz", { .threads = 1} };
    s1 >> buffer;
    //![decompression]
    }
}
