// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <string_view>

inline constexpr std::string_view input =
  R"raw(>ID1
ACGTTTTTTTTTTTTTTT
>ID2
ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>ID3 lala
ACGTTTA
ACGTTTTTTTT

>ID4
ACGTTTA
>ID5 lala
ACGTTTA
ACGTTTTTTTT
)raw";

inline constexpr std::string_view input_bgzipped{
  "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x47\x00\xb3\xf3\x74\x31\xe4\x72\x74\x76\x0f\x41\x05"
  "\x5c\x76\x9e\x2e\x46\x58\xc4\x29\x05\x20\x73\x8d\x15\x72\x12\x73\x12\xa1\x86\x3b\x22\x5b\xc2\x05\x92\x36\x81\xcb\x00"
  "\x39\xa6\xb8\xd5\x02\x00\xcd\x3b\x57\x80\xba\x00\x00\x00\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02"
  "\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00",
  100};
