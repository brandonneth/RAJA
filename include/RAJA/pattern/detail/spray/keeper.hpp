#ifndef RAJA_DETAIL_REDUCE_SPRAY_KEEPER_HPP
#define RAJA_DETAIL_REDUCE_SPRAY_KEEPER_HPP
namespace RAJA {

#if defined(RAJA_ENABLE_OPENMP)

template <typename ContentType, unsigned blocksize> 
class UpdateChunk
{
public:
  ContentType content[blocksize];
  int indices[blocksize];
  int top = 0;
  UpdateChunk<ContentType, blocksize> * next = nullptr;

  UpdateChunk() {}
  
  ContentType * nextRef(int idx) {
    if(this->top < blocksize) {
      ContentType * ret = &this->content[this->top];
      this->content[this->top] = ContentType();
      this->indices[this->top] = idx;
      this->top++;
      return ret;
    } else {
      return nullptr;
    }
  } //nextRef

}; //UpdateChunk


template <typename ReductionPolicy, typename ContentType, unsigned blocksize = 256> 
class KeeperReduction : public virtual SprayReduction<ReductionPolicy,ContentType> 
{
public:

  mutable ContentType * content;
  bool initialized = false;
  int size;
  bool isOrig;
  UpdateChunk<ContentType,blocksize>** myOutgoingUpdates;
  UpdateChunk<ContentType,blocksize>*** allIncomingUpdates;
  int numThreads;
  int threadID;

  KeeperReduction() : threadID(omp_get_thread_num()) {}

  KeeperReduction(int size, ContentType * orig, int numThreads = -1) {
    if(numThreads == -1) {
      this->numThreads = omp_get_max_threads();
    }

    this->content = orig;
    this->initialized = true;
    this->size = size;
    this->isOrig = true;
    this->allIncomingUpdates = new UpdateChunk<ContentType, blocksize>**[this->numThreads];
    auto upd = new UpdateChunk<ContentType, blocksize>*[this->numThreads * this->numThreads];
    for(int i = 0; i < this->numThreads; i++) {
      this->allIncomingUpdates[i] = &upd[i*this->numThreads];
      #pragma omp simd
      for(int j = 0; j < this->numThreads; j++) {
        this->allIncomingUpdates[i][j] = nullptr;
      }
    }
  }

  KeeperReduction(const KeeperReduction<ReductionPolicy, ContentType, blocksize> & other) : 
    content(other.content),
    initialized(true),
    size(other.size),
    isOrig(false),
    numThreads(other.numThreads),
    allIncomingUpdates(other.allIncomingUpdates)
  {
    this->threadID = omp_get_thread_num();
    this->myOutgoingUpdates = new UpdateChunk<ContentType, blocksize>*[other.numThreads]; 
    for(int i = 0; i < other.numThreads; i++) {
      this->myOutgoingUpdates[i] = nullptr;
    }
  }

  ~KeeperReduction() {
    if(this->isOrig) {
      delete[] this->allIncomingUpdates[0];
      delete[] this->allIncomingUpdates;
    } else {
      #pragma omp barrier
      for(int i = 0; i < this->numThreads; i++) {
        auto myIncomingUpdates = this->allIncomingUpdates[this->threadID][i];
        while(myIncomingUpdates) {
          for(int i = 0; i < myIncomingUpdates->top; i++) {
            this->content[myIncomingUpdates->indices[i]] += myIncomingUpdates->content[i];
          }
          auto oldIncomingUpdates = myIncomingUpdates;
          myIncomingUpdates = myIncomingUpdates->next;
          delete oldIncomingUpdates;
        }
      }
      delete[] myOutgoingUpdates;
    }
  }

  ContentType& operator()(int idx) const {
    int owner = idx * this->numThreads / this->size;
    if(owner == this->threadID) {
      return this->content[idx];
    } else {
      if(!this->myOutgoingUpdates[owner]) {
        auto newChunk = new UpdateChunk<ContentType, blocksize>;
        this->myOutgoingUpdates[owner] = newChunk;
        this->allIncomingUpdates[owner][this->threadID] = newChunk;
      }
      ContentType* ret = this->myOutgoingUpdates[owner]->nextRef(idx);
      if(!ret) {
        auto newChunk = new UpdateChunk<ContentType, blocksize>;
        this->myOutgoingUpdates[owner]->next = newChunk;
        this->myOutgoingUpdates[owner] = newChunk;
        ret = this->myOutgoingUpdates[owner]->nextRef(idx);
      }
      return (*ret);
    } 
  }



}; // KeeperReduction


#endif // ompenabled
} // namespace RAJA





#endif
