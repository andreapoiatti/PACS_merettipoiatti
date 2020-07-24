#include "dummy.h"

dummy foo(const dummy & d) {return d;}

int main()
{
std::vector<double> v(100000000,0.1);
dummy d(v),t;
box box(d);
std::cout<<"start test:"<<std::endl;
std::cout<<"pure copy"<<std::endl;

t=box.f1();
std::cout<<"start test print pure copy: ";
//t.printv();

std::cout<<"const ref"<<std::endl;
t=box.f2();
std::cout<<"start test print pure copy: ";
//t.printv();


std::cout<<"pointer to const"<<std::endl;
t=*box.f3();
std::cout<<"start test print pure copy: ";
//t.printv();
//const dummy* dd=box.f3();
//dd->setv();
//t=foo(d);
double pp=0.01;
const double * p=&pp;
//*p=3.3;
return 0;
};
