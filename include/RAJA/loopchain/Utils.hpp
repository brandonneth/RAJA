//Contains camp utilities mostly
#ifndef RAJA_LoopChainUtils_HPP
#define RAJA_LoopChainUtils_HPP



namespace RAJA
{

// camp::idx_seq utility functions
//   - idx_seq_cat
//   - idx_seq_shift
//   - idx_seq_for
//   - idx_seq_from_to

// idx_seq_cat: Concatenates an arbitrary number of sequences

template <camp::idx_t...Is, camp::idx_t...Js>
auto idx_seq_cat(camp::idx_seq<Is...>, camp::idx_seq<Js...>) {
  return camp::idx_seq<Is...,Js...>{};
}

template <camp::idx_t...Is, typename...Seqs>
auto idx_seq_cat(camp::idx_seq<Is...> s1, Seqs&&...seqs) {
  return idx_seq_cat(s1, idx_seq_cat(std::forward<Seqs>(seqs)...));
}

// idx_seq_shift: Shifts each value in the sequence by the amount

template <camp::idx_t amount, camp::idx_t...Is>
auto idx_seq_shift(camp::idx_seq<Is...>) {
  return camp::idx_seq<(Is+amount)...>{};
}

// idx_seq_for: Returns the idx seq for a tuple

template <typename... ElemTypes>
auto idx_seq_for(camp::tuple<ElemTypes...> t) {
  return camp::make_idx_seq_t<sizeof...(ElemTypes)>{};
}

// idx_seq_from_to: Returns an idx seq from the start value up to the end value

template <camp::idx_t start, camp::idx_t end>
auto idx_seq_from_to() {
  auto zeroSeq = camp::make_idx_seq_t<end-start>{};
  
  return idx_seq_shift<start>(zeroSeq);
}





// camp::tuple utility functions
//  - tuple_cat
//  - tuple_slice
//  - tuple_len 


// tuple_cat: Concatenates an arbitrary number of tuples

template <typename Tuple1, typename Tuple2>
auto tuple_cat(Tuple1 t1, Tuple2 t2) {
  return camp::tuple_cat_pair(t1,t2);
}

template <typename Tuple1, typename...Tuples>
auto tuple_cat(Tuple1 t1, Tuples&&...tuples) {
  return tuple_cat(t1, tuple_cat(std::forward<Tuples>(tuples)...));
}

// tuple_slice: Returns a slice of a tuple, with [) bounds as template arguments or with a idx sequence

template <typename Tuple, camp::idx_t...Is>
auto tuple_slice(Tuple t, camp::idx_seq<Is...>) {
  return make_tuple(camp::get<Is>(t)...);
}

template <camp::idx_t start, camp::idx_t end, typename Tuple>
auto tuple_slice(Tuple t) {
  if constexpr (start >= end) {
    return make_tuple();
  } else {
    return tuple_slice(t, idx_seq_from_to<start,end>());
  }
}

// tuple_len: Returns the length of a tuple
template <typename...Ts>
RAJA_INLINE
constexpr auto tuple_len(camp::tuple<Ts...>) {
  return sizeof...(Ts);
}



// Vararg utilities, implemented for both tuples and parameter packs
//   - max
//   - min

template <typename T>
auto max(T val) {
  return val;
}

template <typename T, typename...Ts>
auto max(T val, Ts...rest) {
  return std::max(val, max(rest...));
}

template <typename T>
auto min(T val) {
  return val;
}

template <typename T, typename...Ts>
auto min(T val, Ts...rest) {
  return std::min(val, min(rest...));
}


} //namespace RAJA

#endif
