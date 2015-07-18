/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file MatrixVector.h */

#ifndef AMDIS_MATRIXVECTOR_H
#define AMDIS_MATRIXVECTOR_H

#include "Global.h"
#include "Serializable.h"
#include "traits/basic.hpp" // not dependenciess on other AMDiS types

namespace AMDiS {

  /// Class for efficient vector operations of fixed size numerical vectors.
  template<typename T>
  class Vector : public Serializable
  {
  public:
    
    typedef T       value_type;
    typedef Vector  self;
    
    /// Constructor.
    Vector(int s = 0) 
      : size(s),
	valArray(size ? new T[size] : NULL)
    { }

    /// Copy constructor.
    Vector(self const& other) 
      : Serializable(),
	size(other.size),
	valArray(size ? new T[size] : NULL)
    {
      setValues(other.valArray);
    }

    /// Destructor.
    ~Vector() 
    {
      if (valArray) {
	delete [] valArray; 
	valArray = NULL;
      }
    }

    inline bool used() const 
    {
      return (valArray != NULL);
    }

    /// Change the size of the vector to newSize.
    inline void resize(int s) 
    {
      if (size != s && s > 0) {
	if (valArray) 
	  delete [] valArray;
	valArray = new T[s];
	size = s;
      }
    }

    /// Assignement operator
    self& operator=(self const& other)
    {
      resize( other.getSize() );
      setValues(other.getValArray());
      return *this;
    }

    /// Assignement operator for scalars
    template <typename S>
    typename enable_if< boost::is_convertible<S, T>, Vector<T> >::type &
    operator=(S value) 
    {
      set(value);
      return *this;
    }

    /// Assignement operator
    inline self& operator=(const T* vec) 
    {
      setValues(vec);
      return *this;
    }

    /// Sets all entries to scal.
    template <typename S>
    typename enable_if< boost::is_convertible<S, T> >::type
    set(S value) 
    {
      std::fill(begin(), end(), value);
    }

    /// Sets all entries.
    template <typename S>
    inline void setValues(const S* values) 
    {
      std::copy(values, values + size, begin());
    }

    /// Access to the i-th vector element.
    inline T& operator[](int i) 
    {
      TEST_EXIT_DBG((unsigned)i < (unsigned)size)("Invalid index %d!\n", i);
      return valArray[i];
    }

    /// Access to the i-th vector element for const vectors.
    inline const T& operator[] (int i) const 
    {
      TEST_EXIT_DBG((unsigned)i < (unsigned)size)("Invalid index %d!\n", i);
      return valArray[i];
    }

    /// Returns pointer to the first vector element.
    inline T* begin() { return valArray; }
    inline T const* begin() const { return valArray; }

    /// Returns pointer after the last vector element.
    inline T* end() { return valArray + size; }
    inline T const* end() const { return valArray + size; }

    /// Returns \ref size.
    inline int getSize() const 
    { 
      return size; 
    }

    /// Returns \ref valArray as T-array
    inline T* getValArray() { return valArray; }    
    inline T const* getValArray() const { return valArray; }

    void print() const 
    {
      std::cout << this->size << " vector: " << std::endl;
      for (int i = 0; i < size; i++) 
	std::cout << this->valArray[i] << " ";
      std::cout << std::endl;
    }

    void serialize(std::ostream &out) 
    {
      out.write(reinterpret_cast<const char*>(&size), sizeof(int));
      out.write(reinterpret_cast<const char*>(valArray), size * sizeof(T));
    }

    void deserialize(std::istream &in) 
    {
      in.read(reinterpret_cast<char*>(&size), sizeof(int));
      in.read(reinterpret_cast<char*>(valArray), size * sizeof(T));
    }

    std::string getTypeName() const 
    { 
      return "Vector"; 
    }

  protected:
    /// Size of the vector.
    int size;

    /// Internal storage pointer.
    T *valArray;
  };


  /// Class for efficient matrix operations of fixed size numerical matrices.
  template<typename T>
  class Matrix : public Vector<T>
  {
  public:
    typedef Matrix     self;
    typedef Vector<T>  super;
  
    /// Default constructor.
    Matrix()
      : super(0), rows(0), cols(0)
    {}
    
    /// Constructor.
    Matrix(int r, int c)
      : super(r * c), 
	rows(r),
	cols(c)
    {}
    
    /// copy constructor
    Matrix(self const& other)
      : super(other),
	rows(other.getNumRows()),
	cols(other.getNumCols())
    { }

    /// Changes the size of the matrix to newRows x newCols.
    inline void resize(int newRows, int newCols) 
    {
      if ((newRows != rows) || (newCols != cols)) {
	super::resize(newRows * newCols);
	rows = newRows;
	cols = newCols;
      }
    }
        
    self& operator=(self const& other)
    {
      resize( other.getNumRows(), other.getNumCols() );
      this->setValues(other.getValArray());
      return *this;
    }

    /// Assignement operator for scalars.
    template <typename S>
    typename enable_if< boost::is_convertible<S, T>, self>::type &
    operator=(S value) 
    {
      this->set(value);
      return *this;
    }

    /// Access to i-th matrix row.
    inline T *operator[](int i) 
    {
      return this->valArray + cols * i;
    }

    /// Access to i-th matrix row for constant matrices.
    inline const T *operator[](int i) const 
    {
      return this->valArray + cols * i;
    }

    /// Access to i-th matrix row.
    inline T& operator()(int i, int j) 
    {
      return this->operator[](i)[j];
    }

    /// Access to i-th matrix row for constant matrices.
    inline const T& operator()(int i, int j) const 
    {
      return this->operator[](i)[j];
    }

    /// Returns \ref rows.
    inline int getNumRows() const 
    { 
      return rows; 
    }

    /// Return \ref cols.
    inline int getNumCols() const 
    { 
      return cols; 
    }

    void print() const 
    {
      std::cout << this->rows << " x " << this->cols << " matrix: " << std::endl;
      for (int i = 0; i < rows; i++) {
	for (int j = 0; j < cols; j++)
	  std::cout << this->valArray[i * cols + j] << " ";
	std::cout << std::endl;
      }
    }

  protected:
    /// Number of matrix rows.
    int rows;

    /// Number of matrix columns.
    int cols;
  };
  
  template<typename T>
  inline size_t num_rows(const Vector<T>& vec)
  {
    return static_cast<size_t>(vec.getSize());
  }
  
  template<typename T>
  inline void resize(Vector<T>& vec, size_t newSize)
  {
    vec.resize(newSize);
  }
  
}

#endif // AMDIS_MATRIXVECTOR_H