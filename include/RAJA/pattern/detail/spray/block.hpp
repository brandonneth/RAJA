#ifndef RAJA_DETAIL_REDUCE_SPRAY_BLOCK_HPP
#define RAJA_DETAIL_REDUCE_SPRAY_BLOCK_HPP
#define BSIZE 64
namespace RAJA {



template <typename ReductionPolicy, typename ContentType>
class BlockReducer {};


template <typename ReductionPolicy, typename ContentType> 
class BlockReduction 
{
public:
  mutable ContentType ** blocks;
  #if defined(RAJA_ENABLE_OPENMP)
  omp_lock_t *writelocks;
  #endif
  int nblocks;
  int totalsize;
  bool initialized;
  bool useLocks;
  BlockReduction<ReductionPolicy, ContentType> const * orig = nullptr;
  
  BlockReduction(int totalsize, ContentType * orig, bool useLocks = false) {
    int nblocks = (totalsize - 1) / BSIZE + 1;
    this->blocks = (ContentType **) malloc(sizeof(ContentType *));
    if (this->blocks == NULL) {
      printf("Failed to alloc blocks.\n");
      exit(1);
    }
    if (useLocks) {
      this->writelocks = (omp_lock_t *)malloc(nblocks * sizeof(omp_lock_t));
      if (this->writelocks == NULL) {
        printf("Failed to alloc %d writelocks.\n", nblocks);
        exit(1);
      }
      for (int i = 0; i < nblocks; i++) {
        omp_init_lock(&(this->writelocks[i]));
      }
    }
    this->nblocks = 1;
    this->blocks[0] = orig;
    this->orig = NULL;
    this->totalsize = totalsize;
    this->initialized = true;
    this->useLocks = useLocks;
  }

  BlockReduction() { this->initialized = false; }

  ~BlockReduction() {
    if(orig) {
      BlockReducer<ReductionPolicy,ContentType>{}(orig, this);
    } else {
      free(this->blocks);
      #if defined(RAJA_ENABLE_OPENMP)
      if(orig == nullptr && this->useLocks) {
        for(int i = 0; i < nblocks; i++) {
          omp_destroy_lock(&(this->writelocks[i]));
        }
        free(this->writelocks);
      }
      #endif
    }
  }

  BlockReduction(const BlockReduction &other) 
    : orig(other.orig ? other.orig : &other)
  {
    int totalsize = other.totalsize;
    this->nblocks = (totalsize - 1) / BSIZE + 1;
    this->blocks =
        (ContentType **)malloc(this->nblocks * sizeof(ContentType *));
    if (this->blocks == NULL) {
      printf("Failed to alloc blocks.\n");
      exit(1);
    }
    for (int i = 0; i < this->nblocks; i++) {
      this->blocks[i] = NULL;
    }
    this->totalsize = totalsize;
    this->initialized = true;
  }
  

  static void addBlock(ContentType *curOut, ContentType *curIn, int block, int totalsize, bool outAligned) {
    int startIdx = block * BSIZE;
    int j;
    if(startIdx + BSIZE > totalsize) {
      if (outAligned) {
#pragma omp simd aligned(curIn, curOut : 64)
        for (j = 0; j < totalsize - startIdx; j++) {
          curOut[j] += curIn[j];
        }
      } else {
#pragma omp simd aligned(curIn : 64)
        for (j = 0; j < totalsize - startIdx; j++) {
          curOut[j] += curIn[j];
        }
      }
    } else {
     if (outAligned) {
#pragma omp simd aligned(curIn, curOut : 64)
        for (j = 0; j < BSIZE; j++) {
          curOut[j] += curIn[j];
        }
      } else {
#pragma omp simd aligned(curIn : 64)
        for (j = 0; j < BSIZE; j++) {
          curOut[j] += curIn[j];
        }
      }
    }
  }//addBlock

  void createBlockIfNeeded(int block) const {
    if(this->blocks[block] == NULL) {
      if(this->orig->useLocks && omp_test_lock(&(this->orig->writelocks[block]))) {
        this->blocks[block] = &(this->orig->blocks[0][block*BSIZE]);
      } else {
        this->blocks[block] = 
          (ContentType *)aligned_alloc(64, BSIZE * sizeof(ContentType));
        ContentType *curOut = this->blocks[block];
        for(int j = 0; j < BSIZE; j++) {
          curOut[j] = 0;
        }
      }
    }
  } //createBlockIfNeeded

  ContentType &operator()(int idx) const {
    int idxInBlock = idx % BSIZE;
    int block = idx / BSIZE;
    this->createBlockIfNeeded(block);
    return this->blocks[block][idxInBlock];
  }







}; // BlockReduction

} // namespace RAJA





#endif
