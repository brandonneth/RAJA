//Contains camp utilitisie mostly
#ifndef RAJA_LoopChainUtils_HPP
#define RAJA_LoopChainUtils_HPP

namespace RAJA
{


auto tuple_cat(auto tuple1, auto tuple2) {
  return camp::tuple_cat_pair(tuple1, tuple2);
}
template <typename...T>
auto tuple_cat(auto tuple1, T&&...tuples) {
  return camp::tuple_cat_pair(tuple1, tuple_cat(std::forward<T>(tuples)...));
}

template <camp::idx_t amount, camp::idx_t...Seq>
auto add_to_idx_seq(camp::idx_seq<Seq...>) {
  return camp::idx_seq<(Seq+amount)...>{};

}
template <camp::idx_t start, camp::idx_t end>
auto idx_seq_from_to() {
  auto zeroSeq = camp::make_idx_seq_t<end-start>{};

  return add_to_idx_seq<start>(zeroSeq);
}

template <camp::idx_t...UpTo, camp::idx_t...After>
auto replace_index_with_helper(auto originalTuple, auto replacementTuple, camp::idx_seq<UpTo...>, camp::idx_seq<After...>) {


  auto upToTuple = camp::make_tuple(camp::get<UpTo>(originalTuple)...);

  auto afterTuple = camp::make_tuple(camp::get<After>(originalTuple)...);

  auto tempTuple = camp::tuple_cat_pair(upToTuple, replacementTuple);

  return camp::tuple_cat_pair(tempTuple, afterTuple);
}


template <camp::idx_t ToReplace>
auto replace_index_with(auto originalTuple, auto replacementTuple) {

  constexpr auto size = camp::tuple_size<decltype(originalTuple)>().value;

  auto upToSeq = idx_seq_from_to<0,ToReplace>();
  auto afterSeq = idx_seq_from_to<ToReplace+1, size>();


  return replace_index_with_helper(originalTuple, replacementTuple, upToSeq, afterSeq);;

}

template <camp::idx_t... UpTo, camp::idx_t...After>
auto remove_index_helper(auto originalTuple, camp::idx_seq<UpTo...>, camp::idx_seq<After...>) {

  auto upToTuple = camp::make_tuple(camp::get<UpTo>(originalTuple)...);

  auto afterTuple = camp::make_tuple(camp::get<After>(originalTuple)...);

  return camp::tuple_cat_pair(upToTuple, afterTuple);

}


template <camp::idx_t ToRemove>
auto remove_index(auto originalTuple) {
  constexpr auto size = camp::tuple_size<decltype(originalTuple)>().value;

  auto upToSeq = idx_seq_from_to<0,ToRemove>();
  auto afterSeq = idx_seq_from_to<ToRemove+1, size>();

  return remove_index_helper(originalTuple, upToSeq, afterSeq);
}



template <camp::idx_t... Is>
auto slice_tuple(auto tuple, camp::idx_seq<Is...>) {
  return make_tuple(camp::get<Is>(tuple)...);
}


template <camp::idx_t start, camp::idx_t end> 
auto slice_tuple(auto tuple) {

  if constexpr (start == end) {
    return make_tuple();
  } else {
    auto seq = idx_seq_from_to<start,end>();

    return slice_tuple(tuple, seq);
  }
}

} //namespace RAJA


template <typename...Args>
auto idx_seq_for(camp::tuple<Args...>) {

  return camp::make_idx_seq_t<sizeof...(Args)>{};
}

template <camp::idx_t...Is, camp::idx_t...Js>
auto idx_seq_cat(camp::idx_seq<Is...>, camp::idx_seq<Js...>) {
  return camp::idx_seq<Is...,Js...>{};
}//idx_seq_cat
#endif
