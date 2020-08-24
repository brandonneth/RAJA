//Contains camp utilitisie mostly
#ifndef RAJA_LoopChainUtils_HPP
#define RAJA_LoopChainUtils_HPP

namespace RAJA
{

RAJA_INLINE
auto tuple_cat(auto tuple1, auto tuple2) {
  return camp::tuple_cat_pair(tuple1, tuple2);
}
template <typename...T>
RAJA_INLINE
auto tuple_cat(auto tuple1, T&&...tuples) {
  return camp::tuple_cat_pair(tuple1, tuple_cat(std::forward<T>(tuples)...));
}

template <camp::idx_t amount, camp::idx_t...Seq>
RAJA_INLINE
auto add_to_idx_seq(camp::idx_seq<Seq...>) {
  return camp::idx_seq<(Seq+amount)...>{};

}
template <camp::idx_t start, camp::idx_t end>
RAJA_INLINE
auto idx_seq_from_to() {
  auto zeroSeq = camp::make_idx_seq_t<end-start>{};

  return add_to_idx_seq<start>(zeroSeq);
}

template <camp::idx_t...UpTo, camp::idx_t...After>
RAJA_INLINE
auto replace_index_with_helper(auto originalTuple, auto replacementTuple, camp::idx_seq<UpTo...>, camp::idx_seq<After...>) {


  auto upToTuple = camp::make_tuple(camp::get<UpTo>(originalTuple)...);

  auto afterTuple = camp::make_tuple(camp::get<After>(originalTuple)...);

  auto tempTuple = camp::tuple_cat_pair(upToTuple, replacementTuple);

  return camp::tuple_cat_pair(tempTuple, afterTuple);
}


template <camp::idx_t ToReplace>
RAJA_INLINE
auto replace_index_with(auto originalTuple, auto replacementTuple) {

  constexpr auto size = camp::tuple_size<decltype(originalTuple)>().value;

  auto upToSeq = idx_seq_from_to<0,ToReplace>();
  auto afterSeq = idx_seq_from_to<ToReplace+1, size>();


  return replace_index_with_helper(originalTuple, replacementTuple, upToSeq, afterSeq);;

}

template <camp::idx_t... UpTo, camp::idx_t...After>
RAJA_INLINE
auto remove_index_helper(auto originalTuple, camp::idx_seq<UpTo...>, camp::idx_seq<After...>) {

  auto upToTuple = camp::make_tuple(camp::get<UpTo>(originalTuple)...);

  auto afterTuple = camp::make_tuple(camp::get<After>(originalTuple)...);

  return camp::tuple_cat_pair(upToTuple, afterTuple);

}


template <camp::idx_t ToRemove>
RAJA_INLINE
auto remove_index(auto originalTuple) {
  constexpr auto size = camp::tuple_size<decltype(originalTuple)>().value;

  auto upToSeq = idx_seq_from_to<0,ToRemove>();
  auto afterSeq = idx_seq_from_to<ToRemove+1, size>();

  return remove_index_helper(originalTuple, upToSeq, afterSeq);
}



template <camp::idx_t... Is>
RAJA_INLINE
auto slice_tuple(auto tuple, camp::idx_seq<Is...>) {
  return make_tuple(camp::get<Is>(tuple)...);
}


template <camp::idx_t start, camp::idx_t end> 
RAJA_INLINE
auto slice_tuple(auto tuple) {

  if constexpr (start == end) {
    return make_tuple();
  } else {
    auto seq = idx_seq_from_to<start,end>();

    return slice_tuple(tuple, seq);
  }
}

template <typename...Args>
RAJA_INLINE
auto idx_seq_for(camp::tuple<Args...>) {

  return camp::make_idx_seq_t<sizeof...(Args)>{};
}

template <camp::idx_t...Is, camp::idx_t...Js>
RAJA_INLINE
auto idx_seq_cat(camp::idx_seq<Is...>, camp::idx_seq<Js...>) {
  return camp::idx_seq<Is...,Js...>{};
}//idx_seq_cat


template <camp::idx_t I>
RAJA_INLINE
auto tuple_repeat(auto repeated) {
  if constexpr (I <= 0) {

    return make_tuple();
  } else {
    return tuple_cat(make_tuple(repeated), tuple_repeat<I-1>(repeated));
  }
}


template <typename T>
RAJA_INLINE
auto tuple_repeat(auto repeated) {
  return make_tuple(repeated);
}

template <typename T, typename... Ts> 
RAJA_INLINE
auto tuple_repeat(auto repeated) {
  auto subTuple = tuple_repeat<Ts...>(repeated);

  return tuple_cat(subTuple, make_tuple(repeated));
}


template <typename T>
RAJA_INLINE
auto flatten_tuple(camp::tuple<T> tupleOfTuples) {
  return camp::get<0>(tupleOfTuples);
}

template <typename T, typename...Ts>
RAJA_INLINE
auto flatten_tuple(camp::tuple<T,Ts...> tupleOfTuples) {
   
   return tuple_cat(camp::get<0>(tupleOfTuples), flatten_tuple(remove_index<0>(tupleOfTuples)));
  

}

template <typename...Ts>
RAJA_INLINE
constexpr auto tuple_len(camp::tuple<Ts...>) {
  return sizeof...(Ts);
}
} //namespace RAJA

#endif