#pragma once

#include "Config.hpp"

#if 1 // store pointer to startadress in memoryblock
inline void* aligned_malloc(size_t _size, small_t alignment)
{
  void* p1;
  void** p2;
  size_t offset = alignment - 1 + sizeof(void*);
  if ((p1 = (void*)malloc(_size + offset)) == NULL)
    return NULL;

  p2 = (void**)(((size_t)(p1) + offset) & ~(alignment-1));
  p2[-1] = p1;
  return p2;
}

inline void aligned_free(void* p)
{
  void* p1 = ((void**)p)[-1];
  free(p1);
}
#define AMDIS_ALIGNED_SIZE(size, alignment) ((size) + (alignment) - 1 + sizeof(void*))

#else // store offset to startadress in memoryblock
inline void* aligned_malloc(size_t _size, small_t alignment)
{
  size_t offset = alignment - 1 + sizeof(small_t);
  void* mem = malloc(_size + offset);
  if (mem == NULL)
    return NULL;

  small_t* ptr = reinterpret_cast<small_t*>(((size_t)(mem) + offset) & ~(alignment-1));
  small_t shift = (size_t)(ptr) - (size_t)(mem);
  ptr[-1] = shift;
  return ptr;
}

inline void aligned_free(void* ptr)
{
  small_t* tmp = reinterpret_cast<small_t*>(ptr);
  void* addr = reinterpret_cast<void*>((size_t)(ptr) - tmp[-1]);
  free(addr);
}
#define AMDIS_ALIGNED_SIZE(size, alignment) ((size) + (alignment) - 1 + sizeof(small_t))
#endif
