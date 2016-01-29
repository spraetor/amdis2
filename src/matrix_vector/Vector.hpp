/** \file Vector.hpp */

#pragma once

// std c++ headers
#include <algorithm>		// std::copy
#include <utility>		// std::swap
#include <initializer_list>	// std::initializer_list

// AMDiS headers
#include <Log.h>
#include <traits/traits_fwd.hpp>
#include <traits/basic.hpp>
#include <matrix_vector/expr/base_expr.hpp> // VectorExpr
#include "MatrixVectorBase.hpp"	            // MatrixVectorBase

namespace AMDiS
{
  /// Base class for all vectors.
  /** Provide a MemoryPolicy \p MemoryPolicy and a \p SizePolicy for
   *  automatic size calculation.
   **/
  template <class MemoryPolicy, class SizePolicy = DefaultSizePolicy>
  struct VectorBase
    : public MatrixVectorBase<VectorBase<MemoryPolicy, SizePolicy>, MemoryPolicy>,
      public VectorExpr<VectorBase<MemoryPolicy, SizePolicy>>
  {
    using Self  = VectorBase;
    using Super = MatrixVectorBase<Self, MemoryPolicy>;

    using value_type     = Value_t<Super>;
    using size_type      = Size_t<Super>;

    using pointer        = value_type*;
    using const_pointer  = value_type const*;
    using iterator       = pointer;
    using const_iterator = const_pointer;

  protected:
    using Super::_elements;
    using Super::_size;
    using S = SizePolicy;

  public:
    using Super::set;

    // ----- constructors / assignment -----------------------------------------
  public:
    /// \brief Default constructor.
    /// allocates memory for a vector of size \p s
    explicit VectorBase(size_type s = 0)
      : Super(S::eval(s))
    {}

    /// \brief Constructor with initializer.
    /// allocates memory for a vector of size \p s and sets all
    /// entries to \p value0
    explicit VectorBase(size_type s, value_type value0)
      : Super(S::eval(s))
    {
      set(value0);
    }

    /// Copy constructor
    VectorBase(Self const& other)
      : Super(other._size)
    {
      std::copy(other._elements, other._elements + _size, _elements);
    }

    /// \brief Constructor based on an expression.
    /// Use the assignment operator for expressions to copy values elementwise
    template <class Expr>
    VectorBase(VectorExpr<Expr> const& expr)
      : Super(size(expr))
    {
      this->operator=(expr);
    }

    /// constructor using initializer list
    VectorBase(std::initializer_list<value_type> l)
      : Super(l.size())
    {
      TEST_EXIT_DBG( _size == l.size() )
	("size(Initializer-list) should be size(vector)!");
      std::copy(l.begin(), l.end(), _elements);
    }

    /// destructor
    ~VectorBase() {}

#ifndef _MSC_VER // bug in MSVC <= 2013
    /// copy assignment operator
    Self& operator=(Self const& other)
    {
      this->resize( other._size );
      std::copy(other.begin(), other.end(), _elements);
      return *this;
    }
#endif

    using Super::operator= ;
    using Super::operator+= ;
    using Super::operator-= ;
    using Super::operator*= ;
    using Super::operator/= ;

    // need non-templated arguments in order to eliminate a friend declaration
    // warning in gcc
    friend void swap(VectorBase& first, VectorBase& second)
    {
      using std::swap; // enable ADL
      swap(first._size, second._size);
      swap(first._elements, second._elements);
    }

    // ----- element access functions  -----------------------------------------
  public:
    // import operator() from Super-class
    using Super::operator() ;

    /// Access to the i-th vector element.
    value_type& operator[](size_type i) // [[expects: i < _size]]
    {
      return _elements[i];
    }

    /// Access to the i-th vector element. (const variant)
    value_type const& operator[](size_type i) const // [[expects: i < _size]]
    {
      return _elements[i];
    }

    /// Access to the i-th vector element with index checking.
    value_type& at(size_type i)
    {
      TEST_EXIT_DBG(i < _size)
	("Index " << i << " out of range [0, " << _size << ")!\n");
      return _elements[i];
    }

    /// Access to the i-th vector element with index checking. (const variant)
    value_type const& at(size_type i) const
    {
      TEST_EXIT_DBG(i < _size)
	("Index " << i << " out of range [0, " << _size << ")!\n");
      return _elements[i];
    }

    template <class T1, class T2>
    void setMidpoint(T1 const& t1, T2 const& t2)
    {
      this->operator=(0.5*(t1 + t2));
    }
  };


  /// Size of VectorBase
  template <class M, class S>
  size_t size(VectorBase<M,S> const& vec)
  {
    return vec.getSize();
  }

  /// number of rows of VectorBase
  template <class M, class S>
  size_t num_rows(VectorBase<M,S> const& vec)
  {
    return vec.getSize();
  }

  /// number of columns of VectorBase
  template <class M, class S>
  size_t num_cols(VectorBase<M,S> const& vec)
  {
    return 1;
  }

} // end namespace AMDiS
