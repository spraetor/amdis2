/** \file Matrix.hpp */

#pragma once

// std c++ headers
#include <algorithm>		// std::copy
#include <utility>		// std::swap
#include <initializer_list>	// std::initializer_list

// AMDiS headers
#include <Log.h>
#include <traits/traits_fwd.hpp>
#include <matrix_vector/expr/base_expr.hpp> // MatrixExpr
#include "MatrixVectorBase.hpp"	            // MatrixVectorBase


namespace AMDiS 
{ 
  /// Base class for all matrices.
  /** Provide a MemoryPolicy \p MemoryPolicy and a \p SizePolicy for
   *  automatic size calculation. 
   **/
  template <class MemoryPolicy, class SizePolicy = DefaultSizePolicy>
  struct MatrixBase 
      : public MatrixVectorBase< MatrixBase<MemoryPolicy, SizePolicy>, MemoryPolicy >,
	public MatrixExpr< MatrixBase<MemoryPolicy, SizePolicy> >
  {
    typedef MatrixBase                           self;
    typedef MatrixVectorBase<self, MemoryPolicy>  super;
    
    typedef typename super::value_type     value_type;
    typedef typename super::size_type       size_type;
    
    typedef value_type*                       pointer;
    typedef value_type const*           const_pointer;
    typedef pointer                          iterator;
    typedef const_pointer              const_iterator;
    
  protected:
    using super::_elements;
    using super::_size;
    using S = SizePolicy;
    
  public:
    using super::set;
  
  // ----- constructors / assignment -------------------------------------------
  public:
    /// \brief Default constructor. 
    /// allocates memory for a matrix of size \p r x \p c
    explicit MatrixBase(size_type r = 0, size_type c = 0)
      : super(S::eval(r) * S::eval(c == 0 ? r : c)),
        _rows(S::eval(r)),
        _cols(S::eval(c == 0 ? r : c))
    { }
    
    /// \brief Constructor with initializer.
    /// allocates memory for a matrix of size \p r x \p c and sets all 
    /// entries to \p value0
    explicit MatrixBase(size_type r, size_type c, value_type value0)
      : super(S::eval(r) * S::eval(c)),
        _rows(S::eval(r)),
        _cols(S::eval(c))
    {
      set(value0);
    }
    
    /// Copy constructor
    MatrixBase(self const& other)
      : super(other._size),
        _rows(other._rows),
        _cols(other._cols)
    {
      std::copy(other.begin(), other.end(), _elements);
    }
    
    /// \brief cConstructor based on an expression.
    /// Use the assignment operator for expressions to copy values elementwise
    template <class Expr>
    MatrixBase(MatrixExpr<Expr> const& expr)
      : super(size(expr)),
        _rows(num_rows(expr)),
        _cols(num_cols(expr))
    {
      this->operator=(expr);
    }

    /// destructor
    ~MatrixBase() { }
    
#ifndef _MSC_VER // bug in MSVC <= 2013
    /// copy assignment operator
    self& operator=(self const& other)
    {
      std::copy(other.begin(), other.end(), _elements);
      return *this;
    }
#endif
    
    // import (compound)-assignment operators from super-class
    using super::operator= ;
    using super::operator+= ;
    using super::operator-= ;
    using super::operator*= ;
    using super::operator/= ;
    
    // need non-templated arguments in order to eliminate a friend declaration 
    // warning in gcc
    friend void swap(MatrixBase& first, MatrixBase& second)
    {
      using std::swap; // enable ADL
      swap(first._size, second._size);
      swap(first._rows, second._rows);
      swap(first._cols, second._cols);
      swap(first._elements, second._elements);
    }
  
    /// resize matrix
    void resize(size_type r, size_type c)
    {
      if (r != _rows || c != _cols) {
        super::resize(r*c);
        _rows = r;
        _cols = c;
      }
    }
    
  // ----- element access functions  -------------------------------------------
  public:   
    /// Access to i-th matrix row.
    pointer operator[](size_type i) 
    {
      return _elements + _cols * i;
    }

    /// Access to i-th matrix row for constant matrices.
    const_pointer operator[](size_type i) const 
    {
      return _elements + _cols * i;
    }
    
    /// Access to the i-th vector element.
    value_type& operator()(size_type i, size_type j) 
    {
      return _elements[i * _cols + j];
    }
    
    /// Access to the i-th vector element. (const variant)
    const value_type& operator()(size_type i, size_type j) const 
    {
      return _elements[i * _cols + j];
    }
    
    // contiguous memory access (used by expressions)
    using super::operator() ;
    
    /// Access to the i-th vector element with index checking.
    value_type& at(size_type i, size_type j) 
    {
      TEST_EXIT_DBG(i < _rows && j < _cols)
	("Index (" << i << ", " << j << ") out of range [0, " << _rows << ")x[0," << _cols << ")!\n");
      return this->operator()(i,j); 
    }
    
    /// Access to the i-th vector element with index checking. (const variant)
    const value_type& at(size_type i, size_type j) const
    { 
      TEST_EXIT_DBG(i < _rows && j < _cols)
	("Index (" << i << ", " << j << ") out of range [0, " << _rows << ")x[0," << _cols << ")!\n");
      return this->operator()(i,j); 
    }
    
    /// Return the number of \ref _rows
    size_type getNumRows() const { return _rows; }
    
    /// Return the number of \ref _cols
    size_type getNumCols() const { return _cols; }
    
#if 1
    /// initialize the matrix as a diagonal matrix
    void setDiag(value_type v)
    {
      TEST_EXIT_DBG( _rows == _cols )("Only for square matrices\n");
      set((value_type)(0));
      for (size_t i = 0; i < _rows; ++i)
	this->operator()(i,i) = v;
    }
#endif
  
  // ---------------------------------------------------------------------------
  private:
    size_type _rows;
    size_type _cols;
  };
  
  /// Size of MatrixBase
  template <class M, class S>
  size_t size(MatrixBase<M,S> const& mat)
  {
    return mat.getSize();
  }
  
  /// number of rows of MatrixBase
  template <class M, class S>
  size_t num_rows(MatrixBase<M,S> const& mat)
  {
    return mat.getNumRows();
  }
  
  /// number of columns of MatrixBase
  template <class M, class S>
  size_t num_cols(MatrixBase<M,S> const& mat)
  {
    return mat.getNumCols();
  }
  
} // end namespace AMDiS
