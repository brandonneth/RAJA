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


} //namespace RAJA

#endif
