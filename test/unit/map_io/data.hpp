// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <string_view>

inline constexpr std::string_view input =
  R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	1H7M1D1M1S2H	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1P1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";

inline constexpr std::string_view input_bgzipped{
  "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
  "\xbc\x00\x55\x8e\xcb\x0a\xc2\x30\x14\x44\xd7\xd3\xbf\x28\x15\x1f"
  "\xb5\xd6\xde\x24\xa6\x90\x55\x63\x0b\xa9\x60\x83\x92\xe0\x5e\xb0"
  "\x82\xdb\xae\xf4\xef\x6d\xd3\x8d\x2e\x86\x0b\xc3\x3d\x87\xa9\xda"
  "\x06\x37\xab\x28\x97\x51\xe5\xae\x70\x56\x0d\xfd\x13\x67\xab\xb8"
  "\x88\x86\xfe\xfe\x20\x08\xc2\x54\x11\x24\x81\x1c\x75\xd4\x8c\x39"
  "\xcd\x5d\x01\x5e\x14\xd0\xb5\xf1\x88\x93\x64\x01\xed\xd4\x4b\x31"
  "\xd8\x6e\x3c\x65\xe0\x19\x04\x0b\xbf\x0c\x92\x81\xda\x72\xe6\x1d"
  "\x6b\xff\x0c\xc6\xd4\xde\x58\x6d\x82\x66\xb9\x5a\x6f\x52\xbc\x3f"
  "\xea\xa8\x5c\xc6\x33\x91\x1d\x82\x8a\x43\xf0\x00\x71\x48\x3e\x4f"
  "\xb9\x4c\x53\x42\x82\xf4\x57\x69\x8c\x36\x5e\x7b\x8d\x38\x4e\xb7"
  "\xd9\x2e\xdf\x47\x5f\x01\x82\x81\x27\xeb\x00\x00\x00\x1f\x8b\x08"
  "\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03"
  "\x00\x00\x00\x00\x00\x00\x00\x00\x00",
  217};