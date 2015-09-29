/** \file resize.hpp */

#pragma once

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>
#include <traits/concepts_base.hpp>

namespace AMDiS
{
  namespace concepts
  {
    HAS_MEMBER_GENERATE(change_dim)
    HAS_MEMBER_GENERATE(resize)
    
    template <class T>
    using MtlVecResize 
      = bool_<HAS_MEMBER(change_dim)<Decay_t<T>, void(size_t)>::value>;
      
    template <class T>
    using MtlMatResize 
      = bool_<HAS_MEMBER(change_dim)<Decay_t<T>, void(size_t, size_t)>::value>;
    
    template <class T>
    using StlVecResize 
      = bool_<HAS_MEMBER(resize)<Decay_t<T>, void(size_t)>::value>;
      
  } // end namespace concepts
  
  namespace traits
  {
    template <class T>
    struct resize<T, Requires_t<concepts::Arithmetic<T>> >
    {
      void operator()(T&, size_t, size_t) const
      {
        ERROR_EXIT("Scalars can not be resized!\n");
      }
    };

    // == mtl-vectors ===
    template <class T>
    struct resize<T, Requires_t<concepts::MtlVecResize<T>> >
    {
      void operator()(T& v, size_t r) const
      {
        v.change_dim(r);
      }
    };


    /// change_dim implementation for STL vectors
    template <class T>
    struct resize<T, Requires_t<concepts::StlVecResize<T>> >
    {
      void operator()(T& v, size_t r)
      {
        v.resize(r);
      }
    };


    // === matrices ===
    template <class T>
    struct resize<T, Requires_t<concepts::MtlMatResize<T>> >
    {
      void operator()(T& m, size_t r, size_t c)
      {
        m.change_dim(r,c);
      }
    };

  } // end namespace traits


  /// resize function for vectors
  template <class Collection>
    Requires_t<concepts::Vector<Collection>>
  resize(Collection& c, size_t rows)
  {
    traits::resize<Collection>()(c, rows);
  }

  /// resize function for matrices
  template <class Collection>
    Requires_t<concepts::Matrix<Collection>>
  resize(Collection& c, size_t rows, size_t cols)
  {
    traits::resize<Collection>()(c, rows, cols);
  }

} // end namespace AMDiS

#endif // AMDIS_TYPE_TRAITS_RESIZE_HPP
