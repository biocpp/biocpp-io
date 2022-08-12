#include <bio/io/misc.hpp>

//![vtag]
void foo(bio::io::vtag_t<1>) { /* do one thing     */ }

void foo(bio::io::vtag_t<2>) { /* do another thing */ }

void bar()
{
    foo(bio::io::vtag<1>); // calls first overload
    foo(bio::io::vtag<2>); // calls second overload
}
//![vtag]

//![ttag]
void bax(seqan3::type_list<int>)    { /* do one thing     */ }

void bax(seqan3::type_list<float>)  { /* do another thing */ }

void bat()
{
    bax(bio::io::ttag<int>);   // calls first overload
    bax(bio::io::ttag<float>); // calls second overload
}
//![ttag]

int main()
{}
