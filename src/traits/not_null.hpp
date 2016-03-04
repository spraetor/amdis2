/** \file not_null.hpp */

#pragma once

// std c++ headers
#include <type_traits>
#include <memory>
#include <cassert>
#include <cstddef>

// AMDiS headers
#include "traits/basic.hpp"

// TODO: replace assert by TEST_EXIT_DBG

namespace AMDiS
{  
  /// Helper class that can be used in function argument lists to
  /// indicate that the assigned pointer should be not NULL. This can
  /// be better checked in this class, rather than in the calling function.
  ///
  /// Usage: 
  ///   void foo(not_null<T*> arg); 
  ///   foo(NULL) // --> error
  ///   foo(nullptr) // --> error
  ///   T value; foo(&value) // --> OK, if &value != NULL
  ///
  template <class T, 
    class = Requires_t<std::is_pointer<T>> >
  struct not_null
  {
    using Self    = not_null;
    using PtrType = typename std::pointer_traits<T>::pointer;
    using Element = typename std::pointer_traits<T>::element_type;
    
    // store a pointer and check whether it is NULL
    not_null(PtrType p)
      : pointer(p)
    {
      assert( pointer != NULL );
    }
    
    // copy a pointer that is marked not_null with potential derived 
    // pointer type
    template <class U,
      class = Requires_t<std::is_convertible<T, U>> >
    explicit not_null(not_null<U> const& other)
    {
      *this = other; // redirect to assignment operator
    }
    
    not_null(Self const& other) = default;
    not_null(Self&& other)      = default;
    
    // delete some other possible constructors
    not_null()               = delete;
    not_null(std::nullptr_t) = delete;
    not_null(long)           = delete;
    
    // --------------------------------------------------------------
    
    // assignment of raw pointer with !NULL checking
    Self& operator=(PtrType const& p)
    {
      assert( p != NULL );
      pointer = p;
      return *this;
    }
    
    // copy assignment without !NULL check
    Self& operator=(Self const& other) = default;
    
    template <class U,
      class = Requires_t<std::is_convertible<T, U>> >
    Self& operator=(not_null<U> const& other)
    {
      pointer = other.get();
      return *this;
    }
    
    // move assignment without !NULL check
    Self& operator=(Self&& other)
    {
      pointer = other.get();
      return *this;
    }
    
    template <class U,
      class = Requires_t<std::is_convertible<T, U>> >
    Self& operator=(not_null<U>&& other)
    {
      pointer = other.get();
      return *this;
    }
    
    // don't allow to assign a NULL
    Self& operator=(std::nullptr_t) = delete;
    Self& operator=(int)            = delete;
    
    // --------------------------------------------------------------
    
    Element& operator*()  const { return *pointer; }
    Element* operator->() const { return &*pointer; }
    
    PtrType get() const 
    { 
#ifdef _MSC_VER
      __assume(pointer != NULL);
#endif
      return pointer;
    }
    
    // --------------------------------------------------------------
    
    operator PtrType() { return get(); }
    
    // don't allow comparison for NULL, since this is checked on 
    // construction
    operator bool() const = delete;
    
  private:
    PtrType pointer;
  };

} // end namespace AMDiS
