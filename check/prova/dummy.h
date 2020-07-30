#ifndef DUMMY_
#define DUMMY_
#include <iostream>
#include <vector>

class dummy
{
private:

std::vector<double> v;

public:
dummy(void){std::cout<<"default called"<<std::endl;};
dummy(const std::vector<double> &vv):v(vv){
        std::cout<<"Constructor in action"<<std::endl;
};
dummy(const dummy &dd)
{       this->v=dd.v;
        std::cout<<"Copy constructor in action"<<std::endl;
};
dummy& operator= (const dummy&d) {
        this->v=d.v;
        std::cout<<"Copy assignment called"<<std::endl;
        return *this;
}
void printv(void) {for (unsigned i=0; i<this->v.size();i++)
                     std::cout<<v[i]<<" ";
                    std::cout<<std::endl;}
void setv(std::vector<double> vv) {this->v=vv;};

};

class box
{
public:
dummy d_;
box(const dummy& dum): d_(dum) {};
dummy f1(void) {return this->d_;};
const dummy& f2(void) {return this->d_;};
const dummy * f3(void) {return &this->d_;};
};


#endif
