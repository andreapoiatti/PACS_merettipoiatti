#ifndef __CARRIER__
#define __CARRIER__

#include <iostream>
#include <type_traits>

class Ext1;
class Ext2;

template<typename... Extensions>
class Carrier: public Extensions...
{
        private:
                //Basic Data
                int a=2;
        public:
                bool inherited_from_Ext1 = false;
                bool inherited_from_Ext2 = false;

                Carrier() = default;

                //! Constructor taking object Extensions and initializing with it the new class
                template<typename Dummy = typename std::enable_if<sizeof...(Extensions)!=0, void>, typename... Bricks>
                Carrier(Bricks && ... ext): Extensions(std::forward<Bricks>(ext)) ... {
                        if(std::is_base_of<Ext1, Carrier>::value)
                                inherited_from_Ext1 = true;
                        if(std::is_base_of<Ext2, Carrier>::value)
                                inherited_from_Ext2 = true;
                };

                void print_int(void){std::cout << a << std::endl;}
                void print_inh(void)
                {
                        std::cout << inherited_from_Ext1 << std::endl;
                        std::cout << inherited_from_Ext2 << std::endl;
                }
};

class Ext1
{
        private:
                //Basic Data
                char b='c';
        public:
                Ext1() = default;

                void print_char(void) const {std::cout << b << std::endl;}
};

class Ext2
{
        private:
                //Basic Data
                double c=3.0;
        public:
                Ext2() = default;

                void print_double(void) const {std::cout << c << std::endl;}
};

template <bool B1, bool B2>
struct multi_bool_type
{
    static const bool value1 = B1;
    static const bool value2 = B2;
};

typedef multi_bool_type<true, true> tt_type;
typedef multi_bool_type<false, true> ft_type;
typedef multi_bool_type<true, false> tf_type;
typedef multi_bool_type<false, false> ff_type;

class Utilities_User
{
        public:
                template<typename T>
                static void universal_print(const T & c, tt_type activation)
                {
                        c.print_double();
                        c.print_char();
                }

                template<typename T>
                static void universal_print(const T & c, tf_type activation)
                {
                        c.print_char();
                }

                template<typename T>
                static void universal_print(const T & c, ft_type activation)
                {
                        c.print_double();
                }

                template<typename T>
                static void universal_print(const T & c, ff_type activation)
                {
                        std::cout << "Print nothing" << std::endl;
                }

};


template<typename Container>
class User
{
        private:
                Container c;
        public:
                User(const Container & c_):c(c_) {};

        void print(void)
        {
                multi_bool_type<std::is_base_of<Ext1, Container>::value,std::is_base_of<Ext2, Container>::value> activation;
                Utilities_User::universal_print<Container>(c, activation);
        }

};





#endif
