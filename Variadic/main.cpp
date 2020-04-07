#include <iostream>
#include "carrier.h"

int main(void)
{
        Ext1 ca;
        Ext2 cb;
        Carrier<> a;
        Carrier<Ext2> d(cb);
        Carrier<Ext1> c(ca);
        Carrier<Ext1, Ext2> e(ca,cb);
        User<Carrier<>> ua(a);
        User<Carrier<Ext1>> us(c);
        User<Carrier<Ext2>> ud(d);
        User<Carrier<Ext1, Ext2>> uk(e);

        us.print();
        std::cout<<std::endl;
        ud.print();
        std::cout<<std::endl;
        uk.print();
        std::cout<<std::endl;
        ua.print();

        return 0;
}
