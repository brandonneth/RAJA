#ifndef _RAJA_MODELDECISIONVARIABLE_HPP
#define _RAJA_MODELDECISIONVARIABLE_HPP

#include <sstream>
#include "RAJA/loopchain/FormatDecisions/Model/ConstraintGrammar.hpp"


#include <chrono>
namespace RAJA {

  
template <typename Args, camp::idx_t...Is, typename ViewType>
void f_helper(Args&& args_tuple, camp::idx_seq<Is...>, ViewType & v, ViewType &  v2) {
  v(camp::get<Is>(args_tuple)...) = v2(camp::get<Is>(args_tuple)...);

}

template <typename ViewType>
auto f_lambda(ViewType & v, ViewType & v2) {
  return [&](auto...is) {
    f_helper(camp::forward_as_tuple(is...), idx_seq_from_to<0,sizeof...(is)>(), v, v2);
  };
}


 

  template <camp::idx_t I>
  class FormatDecision : public DecisionVariable<I> {
    public:
    std::string viewName;
    camp::idx_t idx;
    std::array<camp::idx_t, I> permutation;
    
    FormatDecision(std::string viewName_, camp::idx_t idx_, std::vector<camp::idx_t> perm) 
      : viewName(viewName_), idx(idx_)
    {
      for(camp::idx_t i = 0; i < I; i++) {
        permutation[i] = perm.at(i);
      }
    }
    FormatDecision(std::string viewName_, camp::idx_t idx_, std::array<camp::idx_t,I> perm) 
      : viewName(viewName_), idx(idx_), permutation(perm)
    {}


    std::string rep() const {
      std::stringstream s;
      s << viewName << "_f" << idx;
      for(size_t i = 0; i < permutation.size(); i++) {
        s << "_" << permutation[i];
      }
      return s.str();

    }
virtual camp::idx_t estimate_cost() {
  return 0; 
} 

  };
    
  template <camp::idx_t I>
  class ConversionDecision : public DecisionVariable<I> {
   public:

    ConversionDecision(std::string viewName_, camp::idx_t idx_, 
                       std::vector<camp::idx_t> in_perm_, std::vector<camp::idx_t> out_perm_) : 
      viewName(viewName_), 
      idx(idx_) 
    {
      
 
    for(camp::idx_t i = 0; i < I; i++) {
        in_perm[i] = in_perm_.at(i);
        out_perm[i] = out_perm_.at(i);
      }

    }

    std::string viewName;
    camp::idx_t idx;
    std::array<camp::idx_t,I> in_perm;
    std::array<camp::idx_t,I> out_perm;


    virtual std::string rep() const {

      std::stringstream s;
      s << viewName << "_conv" << idx << "_in";
      for(size_t i = 0; i < in_perm.size(); i++) {
        s << "_" << in_perm[i];
      }
      s << "_out";
      for(size_t i = 0; i < out_perm.size(); i++) {
        s << "_" << out_perm[i];
      }
      return s.str();
    }
    std::string key() const {
      std::stringstream s;
      for(size_t i = 0; i < in_perm.size(); i++) {
        s << "_" << in_perm[i];
      }
      s << "_out";
      for(size_t i = 0; i < out_perm.size(); i++) {
        s << "_" << out_perm[i];
      }
      return s.str();

    }

virtual camp::idx_t estimate_cost() {
  static std::map<std::string,idx_t> cache;
  
  if(cache[key()] != 0) {
    return cache[key()];
  }
  if (in_perm == out_perm) {
    cache[key()] = 1;
    return cache[key()];
  }

  using V =  View<double, Layout<I>>;

  camp::idx_t n = 128;

  camp::idx_t s = 1;
  std::array<camp::idx_t, I> sizes;
  for(camp::idx_t i = 0; i < I; i++) {
    sizes[i] = n;
    s *= n;
  }
  Layout<I> in_layout = make_permuted_layout(sizes, in_perm);
  Layout<I> out_layout = make_permuted_layout(sizes, out_perm);
  
  using kpol = typename nest<I>::pol;

  auto start =std::chrono::high_resolution_clock::now(); 
  V in(new double[s], in_layout);

  V out(new double[s], out_layout);

  auto l = f_lambda(in, out);
  auto seg = RangeSegment(0,n);

  auto segs = tuple_repeat<I>(seg);

  auto conversion = make_kernel<kpol>(segs, l);

  conversion(); 
  using T = double;
T * _v = in.get_data();   
  T * _v2 = out.get_data();
  memcpy(_v, _v2, sizeof(double) *in.size());
    free(_v2);
    free(_v);
  auto stop =std::chrono::high_resolution_clock::now();

  auto t =std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count();
  cache[key()] = t;
  //std::cerr << "convcost " << key() << ": " << t << "\n";
  return cache[key()];

} 

  };


template <camp::idx_t NumDims>
TermList<NumDims> * decision_variables(std::string viewName, camp::idx_t idx) {
  auto vars = new TermList<NumDims>();
  auto all_perms = enumerate_perms<NumDims>();
  for(auto perm : all_perms) {
    std::vector<camp::idx_t> p;
    for(int i = 0; i < NumDims; i++) {
      p.push_back(perm[i]);
    }
    vars->add(new FormatDecision<NumDims>(viewName, idx, p));
  }

  //format only for idx=0
  if(idx != 0) {
    for(auto perm1 : all_perms) {
      for(auto perm2 : all_perms) {
        std::vector<camp::idx_t> p1;
        std::vector<camp::idx_t> p2;
        for(int i = 0; i < NumDims; i++) {
          p1.push_back(perm1[i]);
          p2.push_back(perm2[i]);
        }   
        vars->add(new ConversionDecision<NumDims>(viewName, idx, p1, p2));
      }
    }
  }

 


  return vars;
      


}

template <camp::idx_t NumDims>
TermList<NumDims> * all_decision_variables(std::string viewName, camp::idx_t num_knls) {
  auto vars = new TermList<NumDims>();
  auto all_perms = enumerate_perms<NumDims>();

  //FormatDecisions 
  for(int idx = 0; idx <= num_knls+1; idx++) {
    for(auto perm : all_perms) {
      std::vector<camp::idx_t> p;
      for(int i = 0; i < NumDims; i++) {
        p.push_back(perm[i]);
      }
      vars->add(new FormatDecision<NumDims>(viewName, idx, p));
    }
  }

  //ConversionDecisions
  for(int idx = 1; idx <= num_knls+1; idx++) {
    for(auto perm1 : all_perms) {
      for(auto perm2 : all_perms) {
        std::vector<camp::idx_t> p1;
        std::vector<camp::idx_t> p2;
        for(int i = 0; i < NumDims; i++) {
          p1.push_back(perm1[i]);
          p2.push_back(perm2[i]);
        }   
        vars->add(new ConversionDecision<NumDims>(viewName, idx, p1, p2));
      }
    }
  }

  return vars;  
}





} //namespace RAJA

#endif

