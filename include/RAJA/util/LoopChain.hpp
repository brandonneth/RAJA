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

int can_fuse_1d(std::vector<SymAccess> const & a1, std::vector<SymAccess> const & a2);
int can_fuse_2d(std::vector<SymAccess> const & a1, std::vector<SymAccess> const & a2);
int amount_to_shift(std::vector<SymAccess> const & a1, std::vector<SymAccess> const & a2);




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

    friend std::ostream&  operator<< (std::ostream& s, SymIter i) {
        s << i.name;
        return s;
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
        isRead = 0;
        isWrite = 0;
    }

    void set_read() {
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
        //std::cout << "SymAccessList operator =\n";

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

    SymAccessList operator += (double) {
        //std::cout << "SymAccessList operator =\n";

        SymAccessList newList = SymAccessList();
         
        for(SymAccess& a : accesses) {
            a.set_write();
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
  
  const SegmentTuple &segments;
  std::tuple<Bodies...> bodies;

  KernelW(SegmentTuple const & s, Bodies const &... b) : segments(s), bodies(b...) {
  
  }
 
  KernelW(const KernelW & k) = default; 
  KernelW(KernelW && k) = default;
  KernelW & operator=(const KernelW & k) = default;
  KernelW & operator=(KernelW && k) = default;


  std::vector<SymAccess> execute_symbolically() {
    SymIter i = SymIter("i0");
    SymIter j = SymIter("i1");

    auto func = std::get<0>(bodies);


    func(i,j);

    auto accesses1 = (i.accesses);
    auto accesses2 = (j.accesses);

    std::vector<SymAccess> allAccesses = std::vector<SymAccess>();

    for(long unsigned int i = 0; i < accesses1->size(); i++) {allAccesses.push_back(accesses1->at(i));}
    for(long unsigned int i = 0; i < accesses2->size(); i++) {allAccesses.push_back(accesses2->at(i));}
    return allAccesses;
  }

  




 
  template <std::size_t... Is>
  void execute(std::index_sequence<Is...>) {
    RAJA::kernel<PolicyType>(segments, std::get<Is>(bodies)...);
    
  }
  void operator () () {
 	auto seq = std::index_sequence_for<Bodies...>{};
    execute(seq);
  }

};

template <typename PolicyType, typename SegmentTuple, typename... Bodies> 
KernelW<PolicyType,SegmentTuple,Bodies...> makeKernel(SegmentTuple const & segment, Bodies const &... bodies) {

   return KernelW<PolicyType,SegmentTuple,Bodies...>((segment), (bodies)...);
}

template <typename PolicyType, typename SegmentTuple, typename... Bodies1, typename... Bodies2>
auto fuse_kernels(
   KernelW<PolicyType,SegmentTuple,Bodies1...> knl1,
   KernelW<PolicyType,SegmentTuple,Bodies2...> knl2)
{
  //std::cout << "fusing kernels\n";
  
  //static bool canFuse = 1;
    
    
  static auto newlambda = [=] (auto i, auto j) {
     std::get<0>(knl1.bodies)(i,j);
     std::get<0>(knl2.bodies)(i,j);
  };

  static auto newKernel = makeKernel<PolicyType>(knl1.segments,newlambda);
  
  //std::cout << "\nexecuting fused kernel from within fuseKernels\n";

  return newKernel;
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

    Container segment;
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

template <typename ExecPolicy, typename Container, typename LoopBody1>
auto shift_forall(ForAll<ExecPolicy,Container,LoopBody1> knl, const int & amount) {
   
   auto shifted_func = [=](auto i) {knl.func(i-amount);};
   
   auto segment = knl.segment;
   
   long int lb = *segment.begin();
   long int ub = *segment.end();
   
   auto shifted_lb = lb + amount;
   auto shifted_ub = ub + amount;
   auto shifted_segment = RAJA::RangeSegment(shifted_lb,shifted_ub); 
  
   return RAJA::makeForAll<ExecPolicy>(RAJA::RangeSegment(shifted_lb,shifted_ub),shifted_func);

}


template <typename ExecPolicy, typename Container, typename Func1, typename Func2>
auto fuse_unequal_bounds(RAJA::ForAll<ExecPolicy,Container,Func1> forall1, RAJA::ForAll<ExecPolicy,Container,Func2> forall2) {
   auto l1 = *forall1.segment.begin();
   auto u1 = *forall1.segment.end();
 
   auto l2 = *forall2.segment.begin();
   auto u2 = *forall2.segment.end();
  
   


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
        static auto newForAll = fuse(forall1, forall2);
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


std::string read_string(std::vector<SymAccess> a1, std::vector<SymAccess> a2, int dimensions);
std::string write_string(std::vector<SymAccess> a1, std::vector<SymAccess> a2, int dimensions);

isl_stat print_point_val(isl_point * pnt, void * user);
isl_stat collect_point_val(isl_point * pnt, void * user);

template <typename ExecPol, typename Container, typename Func1, typename Func2>
std::vector<int> amount_to_shift_no_nesting(ForAll<ExecPol,Container,Func1> forall1, ForAll<ExecPol,Container,Func2> forall2) {
   
  auto a1 = forall1.execute_symbolically();
  auto a2 = forall2.execute_symbolically();

  auto readString = read_string(a1,a2,1);//.c_str();
  auto writeString = write_string(a1,a2,1);//.c_str();
 
  std::cout << "Calculating shift amount\n";
  std::cout << "read string: " << readString << "\n";
  std::cout << "write string: "<< writeString << "\n";
  
  isl_ctx* ctx = isl_ctx_alloc();
  isl_printer * p = isl_printer_to_file(ctx, stdout);

  const char * domainString = "[l0,u0] -> { L1[i] : i < u0 and i >= l0; L2[i] : i < u0 and i >= l0;}";
  isl_union_set * domain = isl_union_set_read_from_str(ctx,domainString);

  
  isl_union_map * reads = isl_union_map_read_from_str(ctx, readString.c_str());
  isl_union_map * writes = isl_union_map_read_from_str(ctx, writeString.c_str());
  

  //reads = isl_union_map_intersect_domain(reads, isl_union_set_copy(domain));
  //writes = isl_union_map_intersect_domain(writes, isl_union_set_copy(domain));

  isl_union_map * reads_inverse = isl_union_map_reverse(isl_union_map_copy(reads));
  isl_union_map * writes_inverse = isl_union_map_reverse(isl_union_map_copy(writes));
  std::cout << "\nreads map\n";
  p = isl_printer_print_union_map(p, reads);

  std::cout << "\nwrites map\n";
  p = isl_printer_print_union_map(p, writes);

  std::cout << "\nreads map inverse \n";
  p = isl_printer_print_union_map(p, reads_inverse);

  std::cout << "\nwrites map inverse \n";
  p = isl_printer_print_union_map(p, writes_inverse);

  isl_union_map * raw = isl_union_map_apply_range(isl_union_map_copy(writes), isl_union_map_copy(reads_inverse));
  isl_union_map * waw = isl_union_map_apply_range(isl_union_map_copy(writes), isl_union_map_copy(writes_inverse));
  isl_union_map * war = isl_union_map_apply_range(isl_union_map_copy(reads), isl_union_map_copy(writes_inverse));

  std::cout << "\nraw\n";
  p = isl_printer_print_union_map(p, raw);

  std::cout << "\nwar\n";
  p = isl_printer_print_union_map(p, war);

  std::cout << "\nwaw\n";
  p = isl_printer_print_union_map(p, waw);

  isl_union_map * deps = isl_union_map_union(raw, waw);
  deps = isl_union_map_union(deps, war);

  std::cout << "\ndeps\n";
  p = isl_printer_print_union_map(p,deps);
  
  isl_union_set * l1_0 = isl_union_set_read_from_str(ctx, "{L1[0]}");

  std::cout << "\nL1 iteration 0 set\n";
  p = isl_printer_print_union_set(p,l1_0);
 
  isl_union_set * depVectors = isl_union_set_apply(l1_0, deps);

  std::cout << "\ndep vectors\n";
  p = isl_printer_print_union_set(p,depVectors);  

  std::vector<int> * l1_to_l2_deps = new std::vector<int>();
  
  isl_union_set * l2Points = isl_union_set_read_from_str(ctx, "{L2[n]}");

  isl_union_set * l1_to_l2_set = isl_union_set_intersect(isl_union_set_copy(depVectors), l2Points);
  std::cout << "\nl1 to l2 set\n";
  p = isl_printer_print_union_set(p,l1_to_l2_set);
 
  isl_set * set =  isl_set_from_union_set(l1_to_l2_set);

  std::cout <<"\nl1 to l2 as a set set\n";
  p = isl_printer_print_set(p, set);
 
  isl_set_foreach_point(set, collect_point_val, (void*) l1_to_l2_deps);

  
  std::vector<int> depVals = *l1_to_l2_deps;  

  delete l1_to_l2_deps;

  int minAmount = 0;

  for(int val : depVals) {
    if(val < minAmount) {
      minAmount = val;
    }
  } 
 
   std::cout << "\n";
   return std::vector<int>{0,minAmount};
} //amount_to_shift_no_nesting
} //namespace RAJA

#endif
