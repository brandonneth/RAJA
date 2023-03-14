// Contains functions that extract isl polyhedral sets from kernel objects.

#ifndef RAJA_sparse_view_HPP
#define RAJA_sparse_view_HPP
#include <cassert>


namespace RAJA {


template <typename ElmType,size_t NumDims>
struct SparseImpl {

using IdxType = camp::idx_t;
virtual ElmType & operator() (std::vector<IdxType> v) = 0;
  
  virtual std::vector<IdxType> get_permutation() = 0;  
  virtual void* get_data() = 0;  
};

template <typename ElmType, size_t NumDims>
struct COO : public SparseImpl<ElmType, NumDims> {

  using IdxType = camp::idx_t;
  using VectorType = std::vector<IdxType>;
  std::vector<VectorType> indices;
 

  std::vector<ElmType> val;
  std::vector<IdxType> permutation;
  ElmType fill = 0;
  const size_t len;

  COO(std::vector<VectorType> t, std::vector<ElmType> v, std::vector<IdxType> p) 
    : indices(t), val(v), permutation(p), len(val.size()) {
    verify();
  }

  void verify() {
    //check lengths
    //check that its sorted the right way
    //check the permutation is right
    assert(val.size() == row.size());
    assert(val.size() == col.size());
    sort();
  }

  template <IdxType...Is>
  auto make_entry(IdxType n, camp::idx_seq<Is...>) {
    std::vector<IdxType> i{indices[Is][n]...};
    return camp::make_tuple(i, val[n]);
  }

  void sort() {
    auto valTuple = camp::make_tuple(val[0]);
    auto indexTuple = tuple_repeat<NumDims>(indices[0][0]);
    auto entryTuple = tuple_cat(indexTuple, valTuple);
 
    using EntryType = camp::tuple<std::vector<IdxType>, ElmType>;

    std::vector<EntryType> entries;

    //std::cout << "Starting the sort method.len = " << len << "\n";
   
    //doPrint("Row: ", indices[0], "\n");
    //doPrint("Col: ", indices[1], "\n");
    //doPrint("Val: ", val, "\n");

    for(size_t i = 0 ; i < len; i++) {
      EntryType next = make_entry(i, idx_seq_from_to<0,NumDims>());
      //doPrint("entry ", i, ": ", camp::get<0>(next), "\n");
      entries.push_back(next);
    }
 
    auto comp = [=](auto a, auto b) {
      //std::cout  << "comparing\n";
 
      auto aa = camp::get<0>(a);
      auto bb = camp::get<0>(b);
      for(IdxType i = 0; i < NumDims; i++) {
        //std::cout << "Dim " << i << ".\n";
        //doPrint(aa);
        //std::cout << "\n";
        //doPrint(bb);
        //std::cout << "\n";
        if (aa[permutation[i]] != bb[permutation[i]]) {
          return aa[permutation[i]] < bb[permutation[i]];
        }
      }
      return false;
    } ;
    //std::cerr << "sorting\n";
    std::sort(entries.begin(), entries.end(), comp);
    
    for(size_t i = 0; i < len; i++) {
      auto newIndices = camp::get<0>(entries[i]);
      auto newValue = camp::get<1>(entries[i]);
      //doPrint("new entry: ", newIndices, "\n");
      for(size_t j = 0; j < NumDims; j++) {
        indices[j][i] = newIndices[j];
      }
      val[i] = newValue;
    }
  
  }

  ElmType & operator() (std::vector<IdxType> v) { 
    //doPrint("Access:", v, "\n");   
    return find<0>(0,len,v);
  }
  
  
 
  template<size_t CurrDepth>
  ElmType & find(size_t begin, size_t end, std::vector<IdxType> v) { 
    auto currCoord = indices[permutation[CurrDepth]];
    
    //std::cout << "find at depth" << CurrDepth << ". permutation = " << permutation[CurrDepth] << "\n"; 
    auto currBegin = currCoord.begin() + begin;
    auto currEnd = currCoord.begin() + end;
    
    
    auto currFind = std::find(currBegin, currEnd, v[permutation[CurrDepth]]);
    if (currFind == currEnd) {return fill;}
   

    auto currFindEnd = std::find_end(currFind, currEnd, currFind, currFind + 1) + 1;
    
    auto nextBeginDist = begin + std::distance(currBegin, currFind);
    auto nextEndDist = begin + std::distance(currBegin, currFindEnd);

    if constexpr (CurrDepth == NumDims - 1) {
      camp::sink(nextEndDist); // silence unused by set warning
      return val[nextBeginDist];
    } else {
      return find<CurrDepth+1>(nextBeginDist, nextEndDist, v);
    }
  }

  std::vector<IdxType> get_permutation(){ return permutation;};  
  void* get_data(){ return val.data();};  
};


template <typename ElmType, size_t NumDims>
struct  SparseView {
 
  
  SparseImpl<ElmType, NumDims> * impl;

  SparseView(SparseImpl<ElmType, NumDims>* i) : impl(i) {}

  void* get_data() { return impl->get_data();}

  template <typename... Args>
  SymAccessList operator () (SymIterator symIterator, Args... args) const
  {
    std::vector<SymIterator> allIterators{{symIterator, args...}};

    SymAccess thisAccess = SymAccess(impl->get_data(), impl->get_permutation(), allIterators);
    thisAccess.isSparse = true;
    SymAccessList l = SymAccessList(thisAccess);
    return l ;
  }
  template <typename...Indices>
  ElmType & operator() (camp::idx_t i, Indices... is) {

    auto v = std::vector<camp::idx_t>{{i, is...}};
    return impl->operator()(v);
  } 
}; // SparseView<NumDims>


template <typename Tuple, camp::idx_t...Is>
auto tuple_to_vector(Tuple t, camp::idx_seq<Is...>) {
  auto vec = std::vector<std::vector<camp::idx_t>>{{camp::get<Is>(t)...}};
  return vec;
}

template <typename ElmType, typename...Vectors>
auto make_sparse_view(Vectors...vecs) {
  auto tuple = camp::make_tuple(vecs...);
  size_t constexpr numDims = sizeof...(Vectors) - 1;

  auto val = camp::get<numDims>(tuple);
  auto indices = tuple_slice<0,numDims>(tuple);
  
  auto indicesVector = tuple_to_vector(indices, idx_seq_for(indices));

  auto permutation = std::vector<idx_t>();
  for(size_t i = 0; i < numDims; i++) {
    permutation.push_back(i);
  }
  auto impl = new COO<ElmType, numDims>(indicesVector, val, permutation);  

  return SparseView<ElmType, numDims>(impl);
}

template <typename ElmType, typename...Vectors>
auto make_sparse_view_permuted(Vectors...vecs) {
  auto tuple = camp::make_tuple(vecs...);
  size_t constexpr numDims = sizeof...(Vectors) - 2;

  auto val = camp::get<numDims>(tuple);
  auto permutation = camp::get<numDims+1>(tuple);
  auto indices = tuple_slice<0,numDims>(tuple);
  
  auto indicesVector = tuple_to_vector(indices, idx_seq_for(indices));

  auto impl = new COO<ElmType, numDims>(indicesVector, val, permutation);  

  return SparseView<ElmType, numDims>(impl);

}



template <typename T>
RAJA::SparseView<T,2> make_random_sparse_view2(camp::idx_t Ni, camp::idx_t Nj, size_t perRowL, size_t perRowH ) {
  using namespace RAJA;
  std::vector<camp::idx_t> row;
  std::vector<camp::idx_t> col;

  //create a vector that we'll shuffle to use fill each row's entries
 
  auto fullCol = std::vector<camp::idx_t>(Nj);
  for(idx_t i = 0; i < Nj; i++) {fullCol[i] = i;}

  for(idx_t i = 0; i < Ni; i++) {
    //for each row, how many elements?
    auto numInRow = rand() % (perRowH-perRowL) + perRowL;
    auto thisRow = std::vector<camp::idx_t>(numInRow, i);
    
    std::random_shuffle(fullCol.begin(), fullCol.end());
    auto thisCol = std::vector<camp::idx_t>(fullCol.begin(), fullCol.begin() + numInRow);

    std::sort(thisCol.begin(), thisCol.end());

    row.insert(row.end(), thisRow.begin(), thisRow.end());
    col.insert(col.end(), thisCol.begin(), thisCol.end());

    //std::cout << "After inserting the " << i << "th row, we have:\n";
    for(size_t j = 0; j < row.size(); j++) {
    //  std::cout << "(" << row[j] << "," << col[j] << ")  ";
    }
    //std::cout << "\n";
  }
  
  std::vector<T> val(row.size());
  for(size_t i = 0; i < row.size(); i++) {
    T x = (T) rand() / (T) RAND_MAX;
    val[i] = x;
  }
  //std::cout << "creating sparse view. Here's whats passed in:\n";
  for(size_t i = 0; i < col.size(); i++) {
  //    std::cout << "(" << row[i] << "," << col[i] << ", " << val[i] << ")  ";
  }
  //std::cout << "\n";
  return make_sparse_view<T>(row,col,val);
}


} //namespace RAJA
#endif
