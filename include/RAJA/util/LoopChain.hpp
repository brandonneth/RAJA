// Contains code relevant to the use of loop chaining to perform inter-loop optimizations

#ifndef RAJA_loopchain_generic_HPP
#define RAJA_loopchain_generic_HPP
#include "RAJA/config.hpp"



#include "RAJA/util/camp_aliases.hpp"

#include "RAJA/pattern/kernel.hpp"
#include "RAJA/pattern/forall.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

namespace RAJA
{
struct SymAccess;

int can_fuse(std::vector<SymAccess> const & a1, std::vector<SymAccess> const & a2);
struct SymIter {
    std::string name;
    int idx;
    std::shared_ptr<std::vector<SymAccess>> accesses;

    SymIter(std::string str) : name(str), idx(0) {
        accesses = std::make_shared<std::vector<SymAccess>>();
    }
    
    SymIter(int num) : name("placeholder"), idx(num) {
        accesses = std::make_shared<std::vector<SymAccess>>();
    }

    SymIter(const SymIter & other) : name(other.name), idx(other.idx) {
        //std::cout << "SymIter & constructor\n";
        accesses = other.accesses;
    }

    template <typename T>
    SymIter operator * (T const & other) {
        std::stringstream b;
        b << name << "*" << other;
        SymIter newIter = SymIter(b.str());
        newIter.accesses = this->accesses;
        return newIter;
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
    int isRead;
    int isWrite;
    SymAccess(void * v, std::vector<SymIter>& is) {
        view = v;
        iters = is;
        type = "UNSET";

    }

    void set_read() {
        type = "READ";
        isRead = 1;
    }

    void set_write() {
        type = "WRITE";
        isWrite = 1;
    }

    std::string access_string() {
        std::stringstream s;
        for(auto i : iters) {
            s << i.name << ",";
        }
        std::string res = s.str();
        res.pop_back();
        return res;
    }

    friend std::ostream&  operator<< (std::ostream& s, SymAccess a) {
        s << a.type << " " << a.view << " ";
        for(SymIter i : a.iters) {
            s << i.name << " ";
        }
        return s;
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
    
    SymAccessList& operator + (const SymAccessList& other) {
        for(SymAccess a : other.accesses) {
            accesses.push_back(a);
        }
        return *this;
    }
    SymAccessList& operator * (SymAccessList const & other) {
        return operator+(other);
    }
    
/*    template <typename T>
    SymAccessList& operator * (T const &) {
        return arith_op();
    }

    
    template <typename T>
    int operator <= (T const &) {
        for(SymAccess& a : accesses) {
            a.set_read();
        }
        
        for(SymAccess& a : accesses) {
            for(SymIter& i : a.iters) {
                i.accesses->push_back(a);
            }
        }
        return 1;


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
    
    
    */
    SymAccessList operator = (const SymAccessList& other) {
        std::cout << "SymAccessList operator =\n";

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

    void arith_cast() {
        for(SymAccess& a : accesses) {
            a.set_read();
        }
        
        for(SymAccess& a : accesses) {
            for(SymIter& i : a.iters) {
                i.accesses->push_back(a);
            }
        }

    }
/*
    operator long int() {
        arith_cast();
        return 1;
    }
    operator int() {
        arith_cast();
        return 1;
    }
    */
    operator double() {
        arith_cast();
        return 1.0;
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

template <typename PolicyType, typename SegmentTuple, typename... Bodies> 
struct KernelW {
  
  SegmentTuple segments;
  std::tuple<Bodies...> bodies;

  KernelW(SegmentTuple &&s, Bodies &&... b) : segments(s), bodies(b...) {
  
  }
  
  
  void execute_symbolically() {
    SymIter i = SymIter("i0");
    SymIter j = SymIter("i1");

    auto func = std::get<0>(bodies);


    func(i,j);

    auto accesses1 = (i.accesses);

    std::cout << std::size(*accesses1) << "\n";
    std::cout << "access1: " << accesses1->at(0) << "\n";
    std::cout << "access2: " << accesses1->at(1) << "\n";
    
    auto accesses2 = (j.accesses);
    
    std::cout << std::size(*accesses2) << "\n";
    std::cout << "access1: " << accesses2->at(0) << "\n";
    std::cout << "access2: " << accesses2->at(1) << "\n";



  }

  




 
  template <std::size_t... Is>
  void execute(std::index_sequence<Is...>) {
    kernel_param<PolicyType>(segments,RAJA::make_tuple(), std::get<Is>(bodies)...);
    
  }
  void operator () () {
 	auto seq = std::index_sequence_for<Bodies...>{};
    execute(seq);
  }

};

template <typename PolicyType, typename SegmentTuple, typename... Bodies> 
KernelW<PolicyType,SegmentTuple,Bodies...> makeKernel(SegmentTuple&& segment, Bodies &&... bodies) {

   return KernelW<PolicyType,SegmentTuple,Bodies...>(std::forward<SegmentTuple>(segment), std::forward<Bodies>(bodies)...);
}


template <typename PolicyType, typename SegmentTuple, typename... Bodies1, typename... Bodies2, std::size_t... Is1, std::size_t... Is2>
KernelW<PolicyType,SegmentTuple,Bodies1...,Bodies2...> expandKernels(SegmentTuple segments, std::tuple<Bodies1...> bodies1, std::tuple<Bodies2...> bodies2, std::index_sequence<Is1>, std::index_sequence<Is2>) {
  
  return makeKernel<PolicyType>(segments, std::get<Is1>(bodies1)..., std::get<Is2>(bodies2)...);
}

template <typename PolicyType, typename SegmentTuple, typename... Bodies1, typename... Bodies2>
KernelW<PolicyType,SegmentTuple,Bodies1...,Bodies2...>  fuseKernels(
   KernelW<PolicyType,SegmentTuple,Bodies1...> knl1,
   KernelW<PolicyType,SegmentTuple,Bodies2...> knl2)
{
  std::cout << "fusing kernels\n";
  
  static bool canFuse = 1;

  if(canFuse) {
    auto bodies1 = knl1.bodies;
    auto bodies2 = knl2.bodies;
    
    auto seq1 = std::index_sequence_for<Bodies1...>{};
    auto seq2 = std::index_sequence_for<Bodies2...>{};

    auto newKernel =  expandKernels<PolicyType>(knl1.segments, bodies1, bodies2, seq1, seq2 );
    return newKernel;
  }
  
}

/*
 *
 *
 * ForAll Wrapper Class and Functions
 *
 *
*/

template <typename ExecPolicy,
          typename Container,
          typename LoopBody>
struct ForAll {

    const Container &segment;
    LoopBody func;
    constexpr ForAll(Container const & c, LoopBody const & l) :
        segment((c)),
        func(l)
    {
    }

    
    
    ForAll(const ForAll & f) = default;
    ForAll(ForAll && f) = default;
    
    
    std::vector<SymAccess> execute_symbolically() const {
        
        SymIter i = SymIter("i0");
        func(i);
        return  *(i.accesses);
    }
    
    ForAll & operator=(const ForAll & f) = default;
    ForAll & operator=(ForAll && f) = default;
    
    void operator() () const {
        forall<ExecPolicy>(
            segment,
            func);
    }
    
};


template <typename ExecPolicy, typename Container, typename LoopBody>
ForAll<ExecPolicy,Container,LoopBody> makeForAll(Container const & c, LoopBody const & l) {
    return ForAll<ExecPolicy,Container,LoopBody>(
        (c),
        (l));
}


template <typename ExecPolicy,
    typename Container,
    typename Func1,
    typename Func2>
 auto fuse(RAJA::ForAll<ExecPolicy,Container,Func1> forall1, RAJA::ForAll<ExecPolicy,Container,Func2> forall2) {
    auto func = [=] (auto i) {forall1.func(i); forall2.func(i);};
    
    return RAJA::makeForAll<ExecPolicy>(forall2.segment,(func));
}


#include "all-isl.h"

template <typename ExecPolicy, typename Container, typename Func1, typename Func2, typename...Args>
void chain(ForAll<ExecPolicy,Container,Func1> forall1, RAJA::ForAll<ExecPolicy,Container,Func2>  forall2, Args&&...args) {

    static bool fusable = can_fuse(forall1.execute_symbolically(), forall2.execute_symbolically());

    if (fusable) {
        auto newForAll = fuse(forall1, forall2);
        chain(newForAll, std::forward<Args>(args)...);
    } else {
        forall1();
        chain(forall2, std::forward<Args>(args)...);
    }
    
}

template <typename ExecPolicy, typename Container, typename Func1, typename...Args>
void chain(ForAll<ExecPolicy,Container,Func1> forall1) {
    forall1();
}
} //namespace RAJA

#endif
