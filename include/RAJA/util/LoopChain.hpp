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


std::string read_string(std::vector<std::vector<SymAccess>> accessLists, int dimensions);
std::string write_string(std::vector<std::vector<SymAccess>> accessLists, int dimensions);

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


struct Default{};
template <auto> struct Fuse {};
template <auto> struct Shift {};



template <typename ExecPol,typename Container,typename Func, typename ...ForAlls>
auto apply_shift(std::vector<int> shiftAmounts, ForAll<ExecPol,Container,Func> knl, ForAlls&&... knls) {

  std::vector<int> remainingShifts = std::vector<int>();
  for(unsigned long i = 1; i < shiftAmounts.size(); i++ ){
    remainingShifts.push_back(shiftAmounts[i]);
  }

  int shiftAmount = shiftAmounts.at(0);

  auto remainingForAllTuple = apply_shift(remainingShifts, std::forward<ForAlls>(knls)...);
  
  auto newKnl = shift_forall(knl, shiftAmount);

  auto currForAllTuple = std::make_tuple(newKnl);

  auto fullTuple = std::tuple_cat(currForAllTuple,remainingForAllTuple);

  return fullTuple;
}

template <typename ExecPol,typename Container,typename Func>
auto apply_shift(std::vector<int> shiftAmounts, ForAll<ExecPol,Container,Func> knl) {

  int shiftAmount = shiftAmounts.at(0);

  auto newKnl = shift_forall(knl, shiftAmount);
  return std::make_tuple(newKnl);
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



template <typename ExecPolicy, typename Container, typename Func1, typename...Args>
std::vector<int> fusion_overlap_bounds_forall(ForAll<ExecPolicy,Container,Func1> forall1, Args&&... args) {
  std::vector<int> bounds = std::vector<int>();

  bounds.push_back(*forall1.segment.begin());
  bounds.push_back(*forall1.segment.end());

  return fusion_overlap_bounds_forall(bounds, std::forward<Args>(args)...);
}

template <typename ExecPolicy, typename Container, typename Func1, typename...Args>
std::vector<int> fusion_overlap_bounds_forall(std::vector<int> bounds, ForAll<ExecPolicy,Container,Func1> forall1, Args&&... args) {
  
  auto lower = *forall1.segment.begin();
  auto upper = *forall1.segment.end();
  std::cout << "bounds: " << lower << ", " << upper << "\n";
  if (lower > bounds[0]) {
    bounds[0] = lower;
  }
  if(upper < bounds[1]) {
    bounds[1] = upper;
  }
  return fusion_overlap_bounds_forall(bounds, std::forward<Args>(args)...);

}
template <typename ExecPolicy, typename Container, typename Func1, typename...Args>
std::vector<int> fusion_overlap_bounds_forall(std::vector<int> bounds, ForAll<ExecPolicy,Container,Func1> lastKnl) {
  auto lower = *lastKnl.segment.begin();
  auto upper = *lastKnl.segment.end();

  std::cout << "bounds: " << lower << ", " << upper << "\n";
  if (lower > bounds[0]) {
    bounds[0] = lower;
  }
  if(upper < bounds[1]) {
    bounds[1] = upper;
  }

  return bounds;

}
template <typename ExecPol, typename Container, typename Func, typename ...Rest>
std::vector<int> amount_to_shift_foralls(ForAll<ExecPol,Container,Func> forall1, Rest&&... rest){

  std::cout << "In user-called amount to shift foralls\n";
  
  std::vector<std::vector<SymAccess>> accessLists = std::vector<std::vector<SymAccess>>();


  return amount_to_shift_foralls(accessLists, forall1, std::forward<Rest>(rest)...);

}

template<typename ExecPol, typename Container, typename Func, typename ...Rest>
std::vector<int> amount_to_shift_foralls(const std::vector<std::vector<SymAccess>> &accessLists, ForAll<ExecPol,Container,Func> forall1, Rest&&... rest) {
 
  std::cout << "In intermediate amount to shift foralls, size = " << accessLists.size() << "\n"; 
  std::vector<SymAccess> accesses = forall1.execute_symbolically();
  
  std::vector<std::vector<SymAccess>> newAccessLists = std::vector<std::vector<SymAccess>>();

  for(auto lst : accessLists) {
    newAccessLists.push_back(lst);
  }
  newAccessLists.push_back(accesses);

  return amount_to_shift_foralls(newAccessLists, std::forward<Rest>(rest)...);

}

template <typename...Rest>
std::vector<int> amount_to_shift_foralls(const std::vector<std::vector<SymAccess>> &accessLists) {

  std::cout << "Finished collecting access information, performing shift calculation\n";

  std::cout << "Gathered access information for " << accessLists.size() << " kernels\n";

  auto numKernels = accessLists.size();

  auto readString = read_string(accessLists,1);
  auto writeString = write_string(accessLists,1);
  
  std::cout << "reads: " << readString;
  std::cout << "\nwrites: " << writeString;

  isl_ctx* ctx = isl_ctx_alloc();
  isl_printer * p = isl_printer_to_file(ctx, stdout);

  isl_union_map * reads = isl_union_map_read_from_str(ctx, readString.c_str());
  isl_union_map * writes = isl_union_map_read_from_str(ctx, writeString.c_str());

  isl_union_map * reads_inverse = isl_union_map_reverse(isl_union_map_copy(reads));
  isl_union_map * writes_inverse = isl_union_map_reverse(isl_union_map_copy(writes));

  isl_union_map * raw = isl_union_map_apply_range(isl_union_map_copy(writes), isl_union_map_copy(reads_inverse));
  isl_union_map * waw = isl_union_map_apply_range(isl_union_map_copy(writes), isl_union_map_copy(writes_inverse));
  isl_union_map * war = isl_union_map_apply_range(isl_union_map_copy(reads), isl_union_map_copy(writes_inverse));

  isl_union_map * deps = isl_union_map_union(raw, waw);
  deps = isl_union_map_union(deps, war);

  std::vector<std::vector<int>> depOffsetLists = std::vector<std::vector<int>>();
  for(int i = 0; i < numKernels; i++) {
    for(int j = i+1; j < numKernels; j++) {
      std::string srcLoop = "{L" + std::to_string(i) + "[0]}";
      std::string dstLoop = "{L" + std::to_string(j) + "[n]}";
 
      isl_union_set * srcSet = isl_union_set_read_from_str(ctx, srcLoop.c_str());
      isl_union_set * dstSet = isl_union_set_read_from_str(ctx, dstLoop.c_str());

      
      isl_union_set * fromSrc = isl_union_set_apply(srcSet, isl_union_map_copy(deps));
       
      isl_union_set * src2dstUnion = isl_union_set_intersect(isl_union_set_copy(fromSrc), dstSet);
      
      std::cout << "\ndeps from " << i << " to " << j << "\n";
      p = isl_printer_print_union_set(p,src2dstUnion);

     
            std::vector<int> * depOffsetPtr = new std::vector<int>();
      
      if(!isl_union_set_is_empty(src2dstUnion)) {
        isl_set * src2dst =  isl_set_from_union_set(src2dstUnion);
 
        isl_set_foreach_point(src2dst, collect_point_val, (void*) depOffsetPtr);
      }
      std::vector<int> depOffsets = *depOffsetPtr;
      
      delete depOffsetPtr;
      
      depOffsetLists.push_back(depOffsets);
      
    }//j
  }//i

  std::string objFncDmn = "";
  std::string objFncRng = "";
  for(int i = 0; i < numKernels; i++) {
    objFncDmn += "S" + std::to_string(i);
    objFncRng += "S" + std::to_string(i);// + "*S" + std::to_string(i);
    if(i != numKernels-1) {
      objFncDmn += ",";
      objFncRng += "+";
    }
  }
  std::string objFncStr = "{ [" + objFncDmn + "] -> [" + objFncRng + "] }";
  
  std::cout << "\nobjective function string: " << objFncStr << "\n";

  isl_aff * objFnc = isl_aff_read_from_str(ctx, objFncStr.c_str());
  
  std::string startingSetStr = "{ [" + objFncDmn + "] }";
  
  isl_set * shiftSet = isl_set_read_from_str(ctx, startingSetStr.c_str());

  
  for(int i = numKernels-1; i >= 0; i--) {
    for(int j = numKernels-1; j > i; j--) {
      auto depOffsets = depOffsetLists.back();
      depOffsetLists.pop_back();
      if(depOffsets.size() == 0) {continue;}
      
      std::string srcName = "S" + std::to_string(i);
      std::string dstName = "S" + std::to_string(j);
      
      int minOffset = depOffsets[0];
     
      for(auto offset : depOffsets) {
        if(offset < minOffset) {minOffset = offset;}
      }
      
      std::string constraintString = "{ [" + objFncDmn + "]" + ": " + dstName + "-" + srcName + "= " + std::to_string(-1 * minOffset) + "}";
      
      std::cout << "constraint: " << constraintString << "\n";

      isl_set * constraintSet = isl_set_read_from_str(ctx, constraintString.c_str());
      shiftSet = isl_set_intersect(isl_set_copy(shiftSet), constraintSet);
    }
  }

  if(isl_set_is_empty(shiftSet)) {
    return std::vector<int>();
  }

  
  isl_point * point = isl_set_sample_point(shiftSet);


  
  std::cout << "\nsample point\n";
  p = isl_printer_print_point(p,point);

  std::vector<int> shifts = std::vector<int>();

  for(int i = 0; i < numKernels; i++) {
    isl_val * shift = isl_point_get_coordinate_val(point, isl_dim_set, i);

    std::cout << "\ndimension " << i << "\n";
    p = isl_printer_print_val(p, shift);

    auto shiftAmount = isl_val_get_num_si(shift);

    shifts.push_back(shiftAmount);
  }

  //isl_val * minVal = isl_set_min_val(shiftSet, objFnc);

  //std::cout << "minVal\n";
  //p = isl_printer_print_val(p,minVal);
  return shifts;

}


template <typename ExecPol, typename Container, typename Func, typename ...Rest>
void shift_fuse_exec_foralls(ForAll<ExecPol, Container,Func> forall1, Rest&&...rest) {

  std::vector<int> shiftAmounts = amount_to_shift_foralls(forall1, std::forward<Rest>(rest)...);
  if(shiftAmounts.size() == 0) {
    std::cout << "Cannot shift. Executing First kernel and trying again\n";
    forall1();
    //shift_fuse_exec_foralls(std::forward<Rest>(rest)...);
  }

  for(auto amount : shiftAmounts) {
    std::cout << "shift: " << amount << "\n";
  }

  apply_shifts_fuse_exec_foralls(shiftAmounts, forall1, std::forward<Rest>(rest)..., 0);


}


template <typename ExecPol, typename Container, typename Func,  typename...Rest>
void apply_shifts_fuse_exec_foralls(const std::vector<int> & shiftAmounts, ForAll<ExecPol, Container,Func> forall1, Rest&&... rest) {

  std::cout << "applying shifts, " << shiftAmounts.size() << " remaining\n";
  
  std::vector<int> remainingShifts = std::vector<int>();
  for(unsigned long i = 1; i < shiftAmounts.size(); i++ ){
    remainingShifts.push_back(shiftAmounts[i]);
  }

  int shiftAmount = shiftAmounts.at(0);

  auto shiftedKernel = shift_forall(forall1, shiftAmount);

  apply_shifts_fuse_exec_foralls(remainingShifts, std::forward<Rest>(rest)..., shiftedKernel);

}

template < typename...Rest>
void apply_shifts_fuse_exec_foralls(const std::vector<int> & shiftAmounts, int seperator, Rest&&... rest) {
  
  std::cout << "finished applying shifts, determining fusion overlap\n";

  std::vector<int> fusionBounds = fusion_overlap_bounds_forall(std::forward<Rest>(rest)...);

  std::cout << "fusion bounds: " << fusionBounds[0] << ", " << fusionBounds[1] << "\n";


  std::cout << "\nEXECUTING\n\n";

  exec_pre_bounds_foralls(fusionBounds, std::forward<Rest>(rest)...);
  apply_fuse_exec_foralls(fusionBounds, std::forward<Rest>(rest)...);  
  exec_post_bounds_foralls(fusionBounds, std::forward<Rest>(rest)...);

}


template <typename ExecPol, typename Container, typename Func,  typename...Rest>
void exec_pre_bounds_foralls(std::vector<int> bounds, ForAll<ExecPol,Container,Func> knl, Rest&&... rest) {
  std::cout << "executing prebounds\n";

  auto knlLower = *knl.segment.begin();

  if(knlLower < bounds[0]) {
    RAJA::forall<ExecPol>(RAJA::RangeSegment(knlLower,bounds[0]), knl.func);
  }
  
  exec_pre_bounds_foralls(bounds, std::forward<Rest>(rest)...);

}

template <typename ExecPol, typename Container, typename Func,  typename...Rest>
void exec_pre_bounds_foralls(std::vector<int> bounds, ForAll<ExecPol,Container,Func> knl) {
  std::cout << "executing last prebounds\n";

  auto knlLower = *knl.segment.begin();

  if(knlLower < bounds[0]) {
    RAJA::forall<ExecPol>(RAJA::RangeSegment(knlLower,bounds[0]), knl.func);
  }
  

}

template <typename ExecPol, typename Container, typename Func,  typename...Rest>
void exec_post_bounds_foralls(std::vector<int> bounds, ForAll<ExecPol,Container,Func> knl, Rest&&... rest) {
  std::cout << "executing postbounds\n";


  auto knlUpper = *knl.segment.end();

  if(knlUpper > bounds[1]) {
    RAJA::forall<ExecPol>(RAJA::RangeSegment(bounds[1], knlUpper), knl.func);
  }
  exec_post_bounds_foralls(bounds, std::forward<Rest>(rest)...);
}

template <typename ExecPol, typename Container, typename Func,  typename...Rest>
void exec_post_bounds_foralls(std::vector<int> bounds, ForAll<ExecPol,Container,Func> knl) {
  std::cout << "executing last postbounds\n";

  auto knlUpper = *knl.segment.end();

  if(knlUpper > bounds[1]) {
    RAJA::forall<ExecPol>(RAJA::RangeSegment(bounds[1], knlUpper), knl.func);
  }
  

}


template <typename ExecPol, typename Container, typename Func, typename Func2, typename...Rest>
void apply_fuse_exec_foralls(std::vector<int> fusionBounds, ForAll<ExecPol,Container,Func> knl1, ForAll<ExecPol,Container,Func2> knl2,  Rest&&... rest) {
 
  std::cout << "fusing foralls\n"; 
  auto fusedFunc = [=](auto i) {knl1.func(i); knl2.func(i);};

  auto fusedKnl = makeForAll<ExecPol>(knl1.segment, fusedFunc);
  apply_fuse_exec_foralls(fusionBounds, fusedKnl, std::forward<Rest>(rest)...);  
}

template <typename ExecPol, typename Container, typename Func, typename Func2>
void apply_fuse_exec_foralls(std::vector<int> fusionBounds, ForAll<ExecPol,Container,Func> knl1, ForAll<ExecPol,Container,Func2> knl2 ){
  std::cout << "fusing foralls\n"; 
  auto fusedFunc = [=](auto i) {knl1.func(i); knl2.func(i);};

  auto fusedKnl = makeForAll<ExecPol>(knl1.segment, fusedFunc);
  apply_fuse_exec_foralls(fusionBounds, fusedKnl);  
}





template <typename ExecPol, typename Container, typename Func>
void apply_fuse_exec_foralls(std::vector<int> fusionBounds, ForAll<ExecPol,Container,Func> knl1) {

  std::cout << "done fusing, executing fused kernel on bounds " << fusionBounds[0] << "," << fusionBounds[1] << "\n";
  auto newBounds = RAJA::RangeSegment(fusionBounds[0], fusionBounds[1]);

  RAJA::forall<ExecPol>(newBounds, knl1.func); 
}


template <typename Tuple, std::size_t... Is>
std::vector<int> fusion_bounds_foralls(Tuple knls, std::index_sequence<Is...>) {

  std::cout << "Determining the fusion bounds using a tuple of foralls\n";
  int lower = *std::get<0>(knls).segment.begin();
  int upper = *std::get<0>(knls).segment.end();
  
  std::cout << "starting bounds " << lower << "," << upper << "\n";
  ((lower = (*std::get<Is>(knls).segment.begin() > lower ? *std::get<Is>(knls).segment.begin() : lower)), ...);
  ((upper = (*std::get<Is>(knls).segment.end() < upper ? *std::get<Is>(knls).segment.end() : upper)), ...);
  return std::vector<int>{lower,upper};

}

template <typename ExecPol, typename Container, typename Func>
int exec_forall_new_bounds(ForAll<ExecPol, Container, Func> knl, auto bnd1, auto bnd2) {
  
  RAJA::forall<ExecPol>(RAJA::RangeSegment(bnd1,bnd2), knl.func);
  return 1;
}


template <typename Tuple, std::size_t... Is>
auto pre_fuse_foralls(Tuple knls, std::index_sequence<Is...>, std::vector<int> bounds) {
  
  ((*std::get<Is>(knls).segment.begin() < bounds[0] ? exec_forall_new_bounds(std::get<Is>(knls), *std::get<Is>(knls).segment.begin(), bounds[0]) : 0), ...);


}

template <typename Tuple, std::size_t... Is>
auto post_fuse_foralls(Tuple knls, std::index_sequence<Is...>, std::vector<int> bounds) {
  
  ((*std::get<Is>(knls).segment.end() > bounds[1] ? exec_forall_new_bounds(std::get<Is>(knls), bounds[1], *std::get<Is>(knls).segment.end()) : 0), ...);


}



template <typename ExecPol, typename Container, typename Func1, typename Func2, typename... ForAlls>
auto fuse_and_exec_unpacked(std::vector<int> bounds, ForAll<ExecPol,Container, Func1> knl1, ForAll<ExecPol,Container,Func2> knl2, ForAlls&&... knls) {
  std::cout << "fusing and executing unpacked tuple of foralls\n";

  auto newKnl = makeForAll<ExecPol>(knl1.segment, [=](auto i){ knl1.func(i); knl2.func(i);});

  return fuse_and_exec_unpacked(bounds, newKnl, std::forward<ForAlls>(knls)...);
  
}

template <typename ExecPol, typename Container, typename Func>
auto fuse_and_exec_unpacked(std::vector<int> bounds, ForAll<ExecPol,Container, Func> knl) {
  return makeForAll<ExecPol>(RAJA::RangeSegment(bounds[0], bounds[1]), knl.func);
}
template <typename Tuple, std::size_t... Is>
auto fuse_and_exec(Tuple knls, std::index_sequence<Is...>, std::vector<int>bounds) {
  
  return fuse_and_exec_unpacked(bounds, (std::get<Is>(knls))...);



}

//arguments are the kernels, transformations are template arguments
template <typename ...ForAlls>
void chain_foralls(ForAlls&&... knls) {
  
  
  static std::vector<int> shiftAmounts = amount_to_shift_foralls(std::forward<ForAlls>(knls)...);

  static auto shiftedForAlls = apply_shift(shiftAmounts, std::forward<ForAlls>(knls)...);
 
  static auto seq = std::index_sequence_for<ForAlls...>{};



  static auto bounds = fusion_bounds_foralls(shiftedForAlls, seq);
  
  pre_fuse_foralls(shiftedForAlls, seq, bounds);


  static auto fusedKnl = fuse_and_exec(shiftedForAlls, seq, bounds);
  fusedKnl();
  
  post_fuse_foralls(shiftedForAlls, seq, bounds);
}



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
