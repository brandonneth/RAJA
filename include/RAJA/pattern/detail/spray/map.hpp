#ifndef RAJA_DETAIL_REDUCE_SPRAY_MAP_HPP
#define RAJA_DETAIL_REDUCE_SPRAY_MAP_HPP

#include <map>
namespace RAJA {



template <typename ReductionPolicy, typename ContentType>
struct MapReducer {};


template <typename ReductionPolicy, typename ContentType> 
struct MapReduction 
{
using MapType = std::map<unsigned, ContentType>;
public:
  ContentType * denseContent = nullptr;
  mutable MapType mapContent = MapType();
  MapReduction<ReductionPolicy, ContentType> const * parent = nullptr; 

  MapReduction() {}

  MapReduction(ContentType *orig) 
    : denseContent(orig)
  {
    
  }

  MapReduction(const MapReduction &other) 
    : parent(other.parent ? other.parent : &other)
  {
  }

  ~MapReduction() 
  {
    if(parent) {
      MapReducer<ReductionPolicy,ContentType>{}(parent, this);
    }
  }

  ContentType &operator()(int idx) const {
    if (denseContent) { return (denseContent)[idx]; }

    auto it = this->mapContent.find(idx);
    if(it == this->mapContent.end()) {
      it = this->mapContent.insert(it, {idx, ContentType()});
    }
    return it->second;
  }


}; // MapReduction

} // namespace RAJA





#endif
