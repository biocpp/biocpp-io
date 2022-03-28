// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <ranges>
#include <string>

inline std::string const full_example =
  R"(browser position chr7:127471196-127495720
browser hide all
track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"
chr7	127471196	127472363	Pos1	0	+	127471196	127472363	255,0,0
chr7	127472363	127473530	Pos2	0	+	127472363	127473530	255,0,0
chr7	127473530	127474697	Pos3	0	+	127473530	127474697	255,0,0
chr7	127474697	127475864	Pos4	0	+	127474697	127475864	255,0,0
chr7	127475864	127477031	Neg1	0	-	127475864	127477031	0,0,255
chr7	127477031	127478198	Neg2	0	-	127477031	127478198	0,0,255
chr7	127478198	127479365	Neg3	0	-	127478198	127479365	0,0,255
chr7	127479365	127480532	Pos5	0	+	127479365	127480532	255,0,0
chr7	127480532	127481699	Neg4	0	-	127480532	127481699	0,0,255)";

inline std::string const minimal_example =
  R"(chr7	127471196	127472363
chr7	127472363	127473530
chr7	127473530	127474697
chr7	127474697	127475864
chr7	127475864	127477031
chr7	127477031	127478198
chr7	127478198	127479365
chr7	127479365	127480532
chr7	127480532	127481699)";

inline std::string const minimal_example_with_header =
  R"(browser position chr7:127471196-127495720
browser hide all
track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"
chr7	127471196	127472363
chr7	127472363	127473530
chr7	127473530	127474697
chr7	127474697	127475864
chr7	127475864	127477031
chr7	127477031	127478198
chr7	127478198	127479365
chr7	127479365	127480532
chr7	127480532	127481699)";

inline std::string const minimal_example_header_regenerated =
  R"(browser position chr7:127471196-127495720
browser hide all
track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On")";
