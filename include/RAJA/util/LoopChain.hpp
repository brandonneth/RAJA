// Contains code relevant to the use of loop chaining to perform inter-loop optimizations
#ifndef RAJA_loopchain_generic_HPP
#define RAJA_loopchain_generic_HPP
#include "RAJA/pattern/forall.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
namespace RAJA
{

struct SymAccess;

struct SymIter {
    std::string name;
    int idx;
    std::vector<SymAccess> * accesses;

    SymIter(std::string str) : name(str), idx(0) {
        accesses = new std::vector<SymAccess>();
    }

    SymIter(int num) : name("placeholder"), idx(num) {
        accesses = new std::vector<SymAccess>();
    }

    SymIter(const SymIter & other) : name(other.name), idx(other.idx) {
        accesses = other.accesses;
    }

    SymIter operator + (const int & other) {
        std::stringstream b;

        b << name << "+" << other;

        SymIter newIter = SymIter(b.str());
        newIter.accesses = this->accesses;
        return newIter;
    }

    SymIter operator - (const int & other) {
        std::stringstream b;

        b << name << "-" << other;

        SymIter newIter = SymIter(b.str());
        newIter.accesses = this->accesses;
        return newIter;
    }
};

struct SymAccess {

    void * view;
    std::string type;

    std::vector<SymIter> iters;

    SymAccess(void * v, std::vector<SymIter>& is) {
        view = v;
        iters = is;
        type = "UNSET";

    }

    void set_read() {
        type = "READ";
    }

    void set_write() {
        type = "WRITE";
    }

    std::string access_string() {
        std::stringstream s;
        s << type << " " << view << " ";
        for(SymIter i : iters) {
            s << i.name << " ";
        }
        return s.str();
    }

    friend std::ostream&  operator<< (std::ostream& o, SymAccess a) {
        o << a.access_string();
        return o;
    }
};

struct SymAccessList {
    
    std::vector<SymAccess> accesses;
    
    SymAccessList() {
        accesses = std::vector<SymAccess>();
    }
    SymAccessList(SymAccess& a) {
        accesses = std::vector<SymAccess>();
        push_back(a);
    }
    
    void push_back(SymAccess& a) {
        accesses.push_back(a);
    }
    
    SymAccessList operator + (const SymAccessList& other) {
        SymAccessList newList = SymAccessList();
        
        for(SymAccess a : accesses) {
            newList.push_back(a);
        }
        
        for(SymAccess a : other.accesses) {
            newList.push_back(a);
        }
        
        return newList;
    }
    
    SymAccessList& arith_op () {
        return *this;
    }
    SymAccessList& operator + (const int &) {
        return arith_op();
    }
    
    SymAccessList operator + (const double &) {
        return arith_op();
    }
    
    SymAccessList  operator * (const int &) {
        return arith_op();
    }
    
    SymAccessList  operator * (const double &) {
        return arith_op();
    }
    
    SymAccessList operator = (const SymAccessList& other) {
        SymAccessList newList = SymAccessList();
        
        for(SymAccess& a : accesses) {
            a.set_write();
            newList.push_back(a);
        }
        
        for(SymAccess a : other.accesses) {
            a.set_read();
            newList.push_back(a);
        }
        
        for(SymAccess& a : newList.accesses) {
            for(SymIter& i : a.iters) {
                i.accesses->push_back(a);
            }
        }
        return newList;
    }
    
    operator int() const {
        for(SymAccess& a : accesses) {
            a.set_read();
        }
        
        for(SymAccess& a : accesses) {
            for(SymIter& i : a.iters) {
                i.accesses->push_back(a);
            }
        }
        return 0;
    }
    
    operator double() const {
        for(SymAccess& a : accesses) {
            a.set_read();
        }
        
        for(SymAccess& a : accesses) {
            for(SymIter& i : a.iters) {
                i.accesses->push_back(a);
            }
        }
        return 0;
    }
    SymAccessList& arith_assign() {
        for(SymAccess& a : accesses) {
            a.set_write();
        }
        
        for(SymAccess& a : accesses) {
            for(SymIter& i : a.iters) {
                i.accesses->push_back(a);
            }
        }
        return *this;
    }
    
    SymAccessList&  operator = (const int &) {
        return arith_assign();
    }
    
    SymAccessList&  operator = (const double &) {
        return arith_assign();
    }
};

template <typename ExecPolicy,
          typename Container,
          typename LoopBody>
struct ForAll {

    Container segment;
    LoopBody func;
    std::vector<SymAccess> symbolicAccesses;
    ForAll(Container && c, LoopBody && l) :
        segment(std::forward<Container>(c)),
        func(std::forward<LoopBody>(l))
    {
        execute_symbolically();
    }

    
    
    ForAll(const ForAll & f) = default;
    ForAll(ForAll && f) = default;
    
    
    void execute_symbolically() {
        SymIter i = SymIter("i");
        func(i);
        symbolicAccesses = *(i.accesses);
    }
    
    ForAll & operator=(const ForAll & f) = default;
    ForAll & operator=(ForAll && f) = default;
    
    void operator() () {
        forall<ExecPolicy>(
            std::forward<Container>(segment),
            std::forward<LoopBody>(func));
    }
    
};

template <typename ExecPolicy, typename Container, typename LoopBody>
ForAll<ExecPolicy,Container,LoopBody> makeForAll(Container&& c, LoopBody && l) {
    return ForAll<ExecPolicy,Container,LoopBody(
        std::forward<Container>(c),
        std::forward<LoopBody>(l));
}


template <typename ExecPolicy,
    typename Container,
    typename Func1,
    typename Func2>
auto fuse(RAJA::ForAll<ExecPolicy,Container,Func1> &forall1, RAJA::ForAll<ExecPolicy,Container,Func2> &forall2) {
    auto func = [=] (auto i) {forall1.func(i); forall2.func(i);};
    
    return RAJA::makeForAll<ExecPolicy>(std::forward<Container>(forall1.segment),std::forward<decltype(func)>(func));
}


#include "all-isl.h"

template <typename ExecPolicy, typename Container, typename Func1, typename Func2, typename...Args>
void chain(ForAll<ExecPolicy,Container,Func1> &forall1, RAJA::ForAll<ExecPolicy,Container,Func2> &forall2, Args...args) {
    if (can_fuse(forall1.symbolicAccesses, forall2.symbolicAccesses)) {
        auto newForAll = fuse(forall1, forall2);
        return chain(newForAll, args...);
    } else {
        forall1();
        return chain(forall2, args...);
    }
    
}

template <typename ExecPolicy, typename Container, typename Func1, typename Func2, typename...Args>
void chain(ForAll<ExecPolicy,Container,Func1> &forall1) {
    forall1();
}
} //namespace RAJA

#endif
