/** \file Forward.h */

#pragma once

#include "Config.h"		// small_t

namespace AMDiS 
{
  // some forwards declaration
  
  // memory-policies
  template <class T, small_t N, small_t M>   struct MemoryBaseStatic;
  template <class T, bool aligned>          struct MemoryBaseDynamic;
  template <class T, small_t N, small_t M>   struct MemoryBaseHybrid;
    
  // size-policies
  struct DefaultSizePolicy;
  template <size_t S> struct StaticSizePolicy;
  
  // matrix-vector types
  template <class Model, class MemoryPolicy>   struct MatrixVectorBase;
  template <class MemoryPolicy, class SizePolicy>    struct VectorBase;
  template <class MemoryPolicy, class SizePolicy>    struct MatrixBase;
  
} // end namespace AMDiS
