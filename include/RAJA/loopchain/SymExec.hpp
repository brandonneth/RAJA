//Contains code for the symbolic execution of raja kernels

#ifndef RAJA_SymExec_HPP
#define RAJA_SymExec_HPP

#include "RAJA/config.hpp"
#include "RAJA/util/Operators.hpp"
#include "RAJA/util/Layout.hpp"
#include "RAJA/loopchain/PolicyToNestingOrder.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <set>

namespace RAJA
{
struct SymIterator;

struct SymAccess;


struct SymAccessList;


struct SymIterator {
  
  std::string name;
  long int idx;
  std::shared_ptr<std::vector<SymAccess>> accesses; // accesses that use this iterator

  SymIterator() = default;
  SymIterator(std::string str) : name(str), idx(0) {
    accesses = std::make_shared<std::vector<SymAccess>>();
  }

   /*
  SymIterator(long int num) : name("placeholder"), idx(num) {
    accesses = std::make_shared<std::vector<SymAccess>>();
  }
  */
  SymIterator(const SymIterator & other) : name(other.name), idx(other.idx) {
    accesses = other.accesses;
  }

  camp::idx_t index() {
    std::size_t pos = name.find("i");
    std::string indexStart = name.substr(pos+1);
    return std::stoi(indexStart);
  }

  int clear() {
    accesses->clear();
    return 0;
  }
  bool operator == (const int) {
    return 0;
  }
  template <typename T>
  SymIterator operator * (const T & other) {
    std::stringstream b;

    b << name << "*" << other;
    SymIterator newIterator = SymIterator(b.str());
    newIterator.accesses = this->accesses;

    return newIterator;
  }

  template <typename T>
  SymIterator operator + (const T & other) {
    std::stringstream b;

    b << name << "+" << other;
    SymIterator newIterator = SymIterator(b.str());
    newIterator.accesses = this->accesses;

    return newIterator;
  }

  template <typename T>
  SymIterator operator - (const T & other) {
    std::stringstream b;

    b << name << "-" << other;
    SymIterator newIterator = SymIterator(b.str());
    newIterator.accesses = this->accesses;

    return newIterator;
  }  

  template <typename T>
  SymIterator operator / (const T & other) {
    std::stringstream b;

    b << name << "/" << other;
    SymIterator newIterator = SymIterator(b.str());
    newIterator.accesses = this->accesses;

    return newIterator;
  }

  template <typename T>
  SymIterator operator % (const T & other) {
    std::stringstream b;

    b << name << "%" << other;
    SymIterator newIterator = SymIterator(b.str());
    newIterator.accesses = this->accesses;

    return newIterator;
  }

  friend std::ostream& operator<< (std::ostream& s, SymIterator i) {
    s << i.name;
    return s;
  }

}; //SymIterator

template <typename T>
SymIterator operator + (const T & other, const SymIterator & iterator);
template <typename T>
SymIterator operator - (const T & other, const SymIterator & iterator);
template <typename T>
SymIterator operator / (const T & other, const SymIterator & iterator);
template <typename T>
SymIterator operator % (const T & other, const SymIterator & iterator);

template <typename T>
SymIterator operator * (const T & other, const SymIterator & iterator) {
  std::stringstream b;

  b << other << "*" << iterator.name;
  SymIterator newIterator = SymIterator(b.str());
  newIterator.accesses = iterator.accesses;

  return newIterator;
}

template <typename LayoutType> 
std::vector<idx_t> layout_to_perm(const LayoutType & layout) {

  const size_t n_dims = layout.n_dims;
  
  std::vector<idx_t> perm = std::vector<idx_t>(n_dims);

  auto strides = layout.strides;
  std::set<int> strided = std::set<int>();
  for(size_t i = 0; i < n_dims; i++) {
    camp::idx_t biggestIndex = -1;
    camp::idx_t biggestStride = -1;
    for(size_t j = 0; j < n_dims; j++) {
      if(strided.find(strides[j]) != strided.end()) { continue; }
      if(strides[j] > biggestStride) {
        biggestIndex = j;
        biggestStride = strides[j];
      }
    }
    strided.insert(biggestStride);
    perm[i] = biggestIndex;
  }

  return perm;

}//layout_to_perm

struct SymAccess {
  
  void * view;
  std::vector<SymIterator> iterators;
  std::vector<idx_t> layout_permutation;

  bool isRead;
  bool isWrite;
  bool isSparse;
  size_t numIterations;

  template <typename LayoutType>
  SymAccess(void * _view, const LayoutType & layout, std::vector<SymIterator>& _iterators) {
    view = _view;
    layout_permutation = layout_to_perm(layout);
    iterators = _iterators;
    isRead = false;
    isWrite = false;
    isSparse = false;
  }
  
  SymAccess(void * _view, const std::vector<camp::idx_t> & perm, std::vector<SymIterator>& _iterators) {
    view = _view;
    layout_permutation = perm;
    iterators = _iterators;
    isRead = false;
    isWrite = false;
    isSparse =  false;
  }

  void set_read(); 
  
  void set_write();
  void link_to_iterators();

  std::string access_string();
  std::string permuted_access_string();
  operator SymAccessList();

  friend std::ostream& operator<< (std::ostream& s, SymAccess a);
  
  template <typename T>
  bool is_view(T v) {
    return view == v.get_data();
  }

  
  int operator < (const SymAccess & other) {
    std::stringstream a, b;
    a << *this;
    b << other;
   
    return a.str().compare(b.str());
  }

  int operator == (const SymAccess & other) {
    std::stringstream a, b;
    a << *this;
    b << other;
    return 0 == a.str().compare(b.str());
  }

  auto argument_order() {
    std::vector<camp::idx_t> argument_indices;
    for(auto iterator : iterators) {
      argument_indices.push_back(iterator.index());
    }
    return argument_indices;
  }

  template <typename KernelPolicy, camp::idx_t ViewDims, camp::idx_t SearchIdx>
  auto normalized_layout_perm(std::array<camp::idx_t, ViewDims> layout_perm) {
    auto policy_order = policy_to_nesting_order<KernelPolicy, SearchIdx>();
    //std::cerr << "policy to nesting order: ";
    //for(auto i : policy_order) { std::cerr << i << " ";}
    //std::cerr << "\n" << std::flush;
    std::vector<camp::idx_t> policy_applied;
    for(auto index : argument_order()) {
      auto indexof = std::find(policy_order.begin(), policy_order.end(), index) - policy_order.begin();
      policy_applied.push_back(indexof);
    }

    std::vector<camp::idx_t> normalized_access;
    for(auto index : layout_perm) {
      normalized_access.push_back(policy_applied.at(index));
    }

    //Now need to make stuff like (0,3,2) be (0,2,1) and such. 
    
    std::array<camp::idx_t, ViewDims> final_perm;
    for(int i = 0; i < ViewDims; i++)  {
      //count how many elements its bigger than.
      int biggerThan = 0;
      for(int j = 0; j < ViewDims; j++) {
        if(normalized_access.at(i) > normalized_access.at(j)) {
          biggerThan++;
        }
      }
      final_perm[i] = biggerThan;
    }
    return final_perm;
  } 
}; //SymAccess



bool operator < (const SymAccess a, const SymAccess b);
void print_access_list(std::ostream&s, std::vector<SymAccess> accesses, int indent); 
void print_permuted_access_list(std::ostream&s, std::vector<SymAccess> accesses, int indent); 




struct SymAccessList {

  std::vector<SymAccess> accesses;

  SymAccessList() {
    accesses = std::vector<SymAccess>();
  }

  SymAccessList(const SymAccess & a) {
    accesses = std::vector<SymAccess>();
    accesses.push_back(a);
  }

  SymAccessList(const SymAccessList& other) {
    accesses = std::vector<SymAccess>();
    for(SymAccess a : other.accesses) {
      accesses.push_back(a);
    }
  } 
  void push_back(const SymAccess & a) {
    accesses.push_back(a);
  }

  //rhs operations 



   //for "a(i) + 2" like statements
  void num_cast() {
    for(SymAccess& a : accesses) {
      a.set_read();
      a.link_to_iterators();
    }
  }

  /*//for "a(i) + i" like statements
  SymAccessList & operator + (const SymIterator &) {
    for(SymAccess& a : accesses) {
      a.set_read();
    }
    return *this;
  } */
  operator double() {num_cast(); return 1.0;}
  
  
  //assignment operations

  SymAccessList operator = (const SymAccessList & other) {
    SymAccessList newList = SymAccessList();

    for(SymAccess & a : accesses) {
      a.set_write();
      newList.push_back(a);
    }

    for(SymAccess a : other.accesses) {
      a.set_read();
      newList.push_back(a);
    }

    for(SymAccess& a : newList.accesses) {
      a.link_to_iterators();
    }

    return newList;
  }

  SymAccessList symbolic_update_equals (const SymAccessList & other) {
    SymAccessList newList = SymAccessList();

    for(SymAccess& a : accesses) {
      a.set_write();
      a.set_read();
      newList.push_back(a);
    }

    for(SymAccess a : other.accesses) {
      a.set_read();
      newList.push_back(a);
    }

    for(SymAccess& a : newList.accesses) {
      a.link_to_iterators();
    }

    return newList;
  }


  SymAccessList num_assign() {
    SymAccessList newList = SymAccessList();

    for(SymAccess& a : accesses) {
      a.set_write();
      newList.push_back(a);
    }

    for(SymAccess& a : newList.accesses) {
      a.link_to_iterators();
    }

    return newList;
  }

  SymAccessList operator = (int) { return num_assign(); }
  SymAccessList operator = (long int) { return num_assign(); }
  SymAccessList operator = (float) { return num_assign(); }
  SymAccessList operator = (double) { return num_assign(); }
  SymAccessList operator = (SymIterator) { return num_assign(); }
  
  SymAccessList num_update_equals() {
    SymAccessList newList = SymAccessList();

    for(SymAccess& a : accesses) {
      a.set_write();
      a.set_read();
      newList.push_back(a);
    }

    for(SymAccess& a : newList.accesses) {
      a.link_to_iterators();
    }

    return newList;
  }

  SymAccessList operator += (int) { return num_update_equals(); }
  SymAccessList operator += (long int) { return num_update_equals(); }
  SymAccessList operator += (float) { return num_update_equals(); }
  SymAccessList operator += (double) { return num_update_equals(); }
  
  SymAccessList operator -= (int) { return num_update_equals(); }
  SymAccessList operator -= (long int) { return num_update_equals(); }
  SymAccessList operator -= (float) { return num_update_equals(); }
  SymAccessList operator -= (double) { return num_update_equals(); }
 
  SymAccessList operator *= (int) { return num_update_equals(); }
  SymAccessList operator *= (long int) { return num_update_equals(); }
  SymAccessList operator *= (float) { return num_update_equals(); }
  SymAccessList operator *= (double) { return num_update_equals(); }
 
  SymAccessList operator /= (int) { return num_update_equals(); }
  SymAccessList operator /= (long int) { return num_update_equals(); }
  SymAccessList operator /= (float) { return num_update_equals(); }
  SymAccessList operator /= (double) { return num_update_equals(); }
 
  SymAccessList operator %= (int) { return num_update_equals(); }
  SymAccessList operator %= (long int) { return num_update_equals(); }
  SymAccessList operator %= (float) { return num_update_equals(); }
  SymAccessList operator %= (double) { return num_update_equals(); }


}; //SymAccessList


SymAccessList list_op_int(const SymAccessList & l0, int scalar);
SymAccessList list_op_double(const SymAccessList & l0, double scalar);
SymAccessList list_op_iterator(const SymAccessList & l0, const SymIterator & i1);
SymAccessList int_op_list(int scalar, const SymAccessList & l1);
SymAccessList double_op_list(double scalar, const SymAccessList & l1);
SymAccessList iterator_op_list(const SymIterator & i0, const SymAccessList & l1);

// Code for expressions involving two symbolic access lists. 
//This should include operations that involve a symoblic access bc it has a conversion function to symaccesslist.
SymAccessList operator + (const SymAccessList & l0, const SymAccessList & l1);
SymAccessList operator - (const SymAccessList & l0, const SymAccessList & l1);
SymAccessList operator * (const SymAccessList & l0, const SymAccessList & l1);
SymAccessList operator / (const SymAccessList & l0, const SymAccessList & l1);
SymAccessList operator % (const SymAccessList & l0, const SymAccessList & l1);


// functions for list op int
SymAccessList operator + (const SymAccessList & l0, int scalar);
SymAccessList operator - (const SymAccessList & l0, int scalar);
SymAccessList operator * (const SymAccessList & l0, int scalar);
SymAccessList operator / (const SymAccessList & l0, int scalar);
SymAccessList operator % (const SymAccessList & l0, int scalar);

//functions for list op double
SymAccessList operator + (const SymAccessList & l0, double scalar);
SymAccessList operator - (const SymAccessList & l0, double scalar);
SymAccessList operator * (const SymAccessList & l0, double scalar);
SymAccessList operator / (const SymAccessList & l0, double scalar);
SymAccessList operator % (const SymAccessList & l0, double scalar);

//functions for list op const SymIterator
SymAccessList operator + (const SymAccessList & l0, const SymIterator & i1);
SymAccessList operator - (const SymAccessList & l0, const SymIterator & i1);
SymAccessList operator * (const SymAccessList & l0, const SymIterator & i1);
SymAccessList operator / (const SymAccessList & l0, const SymIterator & i1);
SymAccessList operator % (const SymAccessList & l0, const SymIterator & i1);


//functions for int op list
SymAccessList operator + (int scalar, const SymAccessList & l1);
SymAccessList operator - (int scalar, const SymAccessList & l1);
SymAccessList operator * (int scalar, const SymAccessList & l1);
SymAccessList operator * (long int scalar, const SymAccessList & l1);
SymAccessList operator / (int scalar, const SymAccessList & l1);
SymAccessList operator % (int scalar, const SymAccessList & l1);

//functions for double op list
SymAccessList operator + (double scalar, const SymAccessList & l1);
SymAccessList operator - (double scalar, const SymAccessList & l1);
SymAccessList operator * (double scalar, const SymAccessList & l1);
SymAccessList operator / (double scalar, const SymAccessList & l1);


//functions for SymIterator op list
SymAccessList operator + (const SymIterator & i0, const SymAccessList & l1);
SymAccessList operator - (const SymIterator & i0, const SymAccessList & l1);
SymAccessList operator * (const SymIterator & i0, const SymAccessList & l1);
SymAccessList operator / (const SymIterator & i0, const SymAccessList & l1);
SymAccessList operator % (const SymIterator & i0, const SymAccessList & l1);


// Assignments

SymAccessList list_assign_list(const SymAccessList &l0, const SymAccessList & l1);
SymAccessList list_assign_int(const SymAccessList & l0, int scalar);
SymAccessList list_assign_double(const SymAccessList & l0, double scalar);
SymAccessList list_assign_iterator(const SymAccessList & l0, const SymIterator & i1);

SymAccessList list_update_list(const SymAccessList &l0, const SymAccessList & l1);
SymAccessList list_update_int(const SymAccessList & l0, int scalar);
SymAccessList list_update_double(const SymAccessList & l0, double scalar);
SymAccessList list_update_iterator(const SymAccessList & l0, const SymIterator & i1);


} // namespace RAJA






#endif
