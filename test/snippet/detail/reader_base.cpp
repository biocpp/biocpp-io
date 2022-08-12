#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <bio/io/seq_io/reader.hpp>

#include "../../unit/seq_io/data.hpp"

void process_read_pair(auto&&, auto&&) {}

int main()
{
    //================= PRE ==========================
    {
        std::ofstream os{"example.fastq", std::ios::binary};
        os << interleaved_fastq;
    }

    //================= SNIPPETS ======================
{
//![read_pair_processing]
// choose deep records so they can be copied/moved
bio::seq_io::reader_options options{ .field_types = bio::seq_io::field_types<bio::ownership::deep> };

// open an interleaved paired-end FastQ file
bio::seq_io::reader reader{"example.fastq", options};

// ask the reader for its record_type; create a variable to hold previous record
decltype(reader)::record_type last_record;

bool is_first_of_pair = true;
for (auto & current_record : reader)
{
    if (is_first_of_pair)
    {
        // backup the current record; only possible because it is deep
        std::swap(current_record, last_record);
    }
    else // is second of pair
    {
        // do something with the pair
        process_read_pair(last_record, current_record);
    }

    is_first_of_pair = !is_first_of_pair;
}
//![read_pair_processing]

}

    //================= POST ==========================
    std::filesystem::remove("example.fastq");
}
