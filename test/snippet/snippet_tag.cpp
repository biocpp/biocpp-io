#include <bio/io/misc.hpp>

//![vtag]
void foo(bio::vtag_t<1>) { /* do one thing     */ }

void foo(bio::vtag_t<2>) { /* do another thing */ }

void bar()
{
    foo(bio::vtag<1>); // calls first overload
    foo(bio::vtag<2>); // calls second overload
}
//![vtag]

//![ttag]
void bax(seqan3::type_list<int>)    { /* do one thing     */ }

void bax(seqan3::type_list<float>)  { /* do another thing */ }

void bat()
{
    bax(bio::ttag<int>);   // calls first overload
    bax(bio::ttag<float>); // calls second overload
}
//![ttag]

int main()
{}
