// Contains functions that extract isl polyhedral sets from kernel objects.

#ifndef RAJA_sparse_view_HPP
#define RAJA_sparse_view_HPP

namespace RAJA {




template <typename ElmType>
class SparseView {
  using IdxType = camp::idx_t;

  std::vector<IdxType> row;
  std::vector<IdxType> col;
  std::vector<ElmType> val;
  ElmType fill = 0;
  IdxType len;

  ElmType & operator() (IdxType i0, IdxType i1) {
    auto rowFind = std::find(begin(row), end(row), i0);
    
    if (rowFind == std::end(row)) {return fill;}
    auto rowFindStart = std::find(row.begin(), row.end(), i0);
    if(rowFindStart == row.end()) {return fill;}

    auto rowFindEnd = std::find_end(rowFindStart, row.end(), rowFindStart, rowFindStart+1);

    auto startDist = std::distance(row.begin(), rowFindStart);
    auto endDist = std::distance(row.begin(), rowFindEnd);

    auto colFind = std::find(col.begin() + startDist, col.begin() + endDist, i1);

    if (colFind == col.begin() + endDist) {return fill;}

    auto idx = std::distance(col.begin(), colFind);

    return val(idx);
  }
  
} // class SparseView


} //namespace RAJA
#endif
