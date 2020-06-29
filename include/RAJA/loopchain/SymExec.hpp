//Contains code for the symbolic execution of raja kernels

#ifndef RAJA_SymExec_HPP
#define RAJA_SymExec_HPP

#include "RAJA/config.hpp"

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

namespace RAJA
{
struct SymIterator;
struct SymAccess;
struct SymAccessList;

struct SymIterator {
  
  std::string name;
  long int idx;
  std::shared_ptr<std::vector<SymAccess>> accesses; // accesses that use this iterator

  SymIterator(std::string str) : name(str), idx(0) {
    accesses = std::make_shared<std::vector<SymAccess>>();
  }

  SymIterator(long int num) : name("placeholder"), idx(num) {
    accesses = std::make_shared<std::vector<SymAccess>>();
  }

  SymIterator(const SymIterator & other) : name(other.name), idx(other.idx) {
    accesses = other.accesses;
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


struct SymAccess {
  
  void * view;
  std::vector<SymIterator> iterators;
  bool isRead;
  bool isWrite;

  SymAccess(void * _view, std::vector<SymIterator>& _iterators) {
    view = _view;
    iterators = _iterators;
    isRead = false;
    isWrite = false;
  }

  void set_read() {
    isRead = true;
  }

  void set_write() {
    isWrite = true;
  }

  void link_to_iterators() {
    for(SymIterator i : iterators) {
      i.accesses->push_back(*this);
    }
  }

  std::string access_string() {
    std::stringstream s;
    for(auto i : iterators) {
      s << i.name << ",";
    }
    
    std::string res = s.str();
    res.pop_back();
    return res;
  }

  friend std::ostream& operator<< (std::ostream& s, SymAccess a) {
    s << " " << a.view << " ";
    s << a.access_string();
    return s;
  }
}; //SymAccess
  

struct SymAccessList {

  std::vector<SymAccess> accesses;

  SymAccessList() {
    accesses = std::vector<SymAccess>();
  }

  SymAccessList(const SymAccess & a) {
    accesses = std::vector<SymAccess>();
    accesses.push_back(a);
  }

  
  void push_back(const SymAccess & a) {
    accesses.push_back(a);
  }

  //rhs operations 

  SymAccessList & arith_operator(const SymAccessList & other) {

    for(SymAccess a : other.accesses) {
      accesses.push_back(a);
    }
    return *this;
  }

  SymAccessList & operator + (const SymAccessList& other) {
    return arith_operator(other);
  }  
 
  SymAccessList & operator - (const SymAccessList & other) {
    return arith_operator(other);
  }

  SymAccessList & operator * (const SymAccessList & other) {
    return arith_operator(other);
  }

  SymAccessList & operator / (const SymAccessList & other) {
    return arith_operator(other);
  }

  SymAccessList & operator % (const SymAccessList & other) {
    return arith_operator(other);
  }

  //for "a(i) + 2" like statements
  void num_cast() {
    for(SymAccess& a : accesses) {
      a.set_read();
      a.link_to_iterators();
    }
  }

  //for "a(i) + i" like statements
  SymAccessList & operator + (const SymIterator & i) {
    for(SymAccess& a : accesses) {
      a.set_read();
    }
    return *this;
  } 
  operator int() {num_cast(); return 1;}
  operator long int() {num_cast(); return 1;}
  operator float() {num_cast(); return 1.0;}
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


} // namespace RAJA






#endif
