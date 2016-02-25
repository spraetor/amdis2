/** \file FixVec.h */

#pragma once

#include <iostream>

#include "AMDiS_fwd.hpp"
#include "AMDiS_base.hpp"
#include "Global.hpp"
#include "Math.hpp"
#include "MatrixVector.hpp"

namespace AMDiS
{

  /// determines how to initialize a FixVec when constructed.
  enum InitType
  {
    NO_INIT = 0,       /**< no initialisation */
    VALUE_LIST = 1,    /**< a complete value list is given */
    DEFAULT_VALUE = 2, /**< all values ar set to a given default value */
    UNIT_VECTOR = 3,   /**< the i-th value is 1, all other values are 0 */
    UNIT_MATRIX = 4    /**< values at the diagonal of a matrix are set to one */
  };

  /** \ingroup Common
   * \brief
   * Contains an vector of FixVecs of same type. To simply allocate an array of
   * FixVecs
   * isn't possible, because the FixVec constructor normally needs at least
   * the corresponding dimension. So you must create an array of FixVec pointers
   * and call the constructor of each FixVec manually. When you use
   * VectorOfFixVecs, this work is done by VectorOfFixVecs's constructor.
   */
  template <class FixVecType>
  class VectorOfFixVecs
  {
  public:
    /// constructs a VectorOfFixVecs without initialisation. dim is passed to
    /// FixVec's constructors. size_ is the number of contained FixVecs. initType
    /// must be NO_INIT.
    VectorOfFixVecs(int d, int s, InitType DBG_VAR( initType ))
      : size(s),
        dim(d)
    {
      TEST_EXIT_DBG(initType == NO_INIT)("wrong initType or wrong initializer\n");

      vec.resize(size);
      for (int i = 0; i < size; i++)
        vec[i] = new FixVecType(dim);
    }

    /// constructs a VectorOfFixVecs via an value list.  dim is passed to
    /// FixVec's constructors. size_ is the number of contained FixVecs. initType
    /// must be VALUE_LIST. ini contains the initialisation values.
    VectorOfFixVecs(int d, int s, InitType DBG_VAR( initType ), FixVecType const* ini)
      : size(s),
        dim(d)
    {
      TEST_EXIT_DBG(initType == VALUE_LIST)("wrong initType or wrong initializer\n");

      vec.resize(size);
      for (int i = 0; i < size; i++)
        vec[i] = new FixVecType(ini[i]);
    }

    /// constructs a VectorOfFixVecs with an default value.  dim is passed to
    /// FixVec's constructors. size_ is the number of contained FixVecs. initType
    /// must be DEFAULT_VALUE. All entries are set to ini.
    VectorOfFixVecs(int d, int s, InitType DBG_VAR( initType ), const FixVecType& ini)
      : size(s),
        dim(d)
    {
      TEST_EXIT_DBG(initType == DEFAULT_VALUE)("wrong initType or wrong initializer\n");

      vec.resize(size);
      for (int i = 0; i < size; i++)
        vec[i] = new FixVecType(ini);
    }

    /// Copy constructor
    VectorOfFixVecs(const VectorOfFixVecs<FixVecType>& rhs)
    {
      size = rhs.getSize();
      dim = rhs.getDim();

      vec.resize(size);
      for (int i = 0; i < size; i++)
        vec[i] = new FixVecType(*(rhs.vec[i]));
    }

    /// Destructor
    ~VectorOfFixVecs()
    {
      for (int i = 0; i < size; i++)
        delete vec[i];

      vec.clear();
    }

    /// Allows the access like in a normal array via []. Used for const objects.
    inline const FixVecType& operator[](int index) const
    {
      TEST_EXIT_DBG(index >= 0 && index < size)("invalid index\n");
      return *(vec[index]);
    }

    /// Allows the access like in a normal array via [].
    inline FixVecType& operator[](int index)
    {
      TEST_EXIT_DBG(index >= 0 && index < size)("invalid index\n");
      return *(vec[index]);
    }

    /// Assignment operator
    VectorOfFixVecs<FixVecType>&
    operator=(const VectorOfFixVecs<FixVecType>& rhs)
    {
      TEST_EXIT_DBG(size == rhs.size)("vectors of different size\n");
      if (this != &rhs)
      {
        for (int i = 0; i < size; i++)
          *(vec[i]) = *(rhs.vec[i]);
      }
      return *this;
    }

    /// Resize the vector
    inline void resize(int newsize)
    {
      vec.resize(newsize);
      for (int i = size; i < newsize; i++)
        vec[i] = new FixVecType(dim);
      size = newsize;
    }

    /// Returns the \ref size of this VectorOfFixVecs
    inline int getSize() const
    {
      return size;
    }

    /// Returns \ref dim
    inline int getDim() const
    {
      return dim;
    }

    /// Returns the size of the contained FixVecs
    inline int getSizeOfFixVec() const
    {
      return vec[0]->getSize();
    }

  protected:
    /// number of contained FixVecs
    int size;

    /// dimension of world
    int dim;

    /// array of pointers to FixVecs
    std::vector<FixVecType*> vec;
  };


  /** \ingroup Common
   * \brief
   * Like the class VectorOfFixVecs contains a vector of FixVecs, this class
   * contains a matrix of FixVecs of same type.
   */
  template<typename FixVecType>
  class MatrixOfFixVecs
  {
  public:
    /// Constructs the matrix without initialisation. r is the number of rows,
    /// c is the number of columns. The other parameters are analog to the
    /// VectorOfFixVecs constructors.
    MatrixOfFixVecs(int dim, int r, int c, InitType DBG_VAR( initType ))
      : rows(r),
        columns(c)
    {
      TEST_EXIT_DBG(initType == NO_INIT)("wrong initType or wrong initializer\n");
      vec = new VectorOfFixVecs<FixVecType>* [rows];
      for (VectorOfFixVecs<FixVecType>** i = &vec[0]; i < &vec[rows]; i++)
        *i = new VectorOfFixVecs<FixVecType>(dim, columns, NO_INIT);
    }

    /// destructor
    ~MatrixOfFixVecs()
    {
      for (VectorOfFixVecs<FixVecType>** i = &vec[0]; i < &vec[rows]; i++)
        delete *i;

      delete [] vec;
    }

    /// assignment operator
    inline VectorOfFixVecs<FixVecType>& operator[](int index)
    {
      TEST_EXIT_DBG(index >= 0 && index < rows)("invalid index\n");
      return *(vec[index]);
    }

    /// assignment operator
    inline const VectorOfFixVecs<FixVecType>& operator[](int index) const
    {
      TEST_EXIT_DBG(index >= 0 && index < rows)("invalid index\n");
      return *(vec[index]);
    }

    /// Returns \ref rows
    inline int getNumberOfRows() const
    {
      return rows;
    }

    /// Returns \ref columns
    inline int getNumberOfColumns() const
    {
      return columns;
    }

  private:
    /// number of rows of the matrix
    int rows;

    /// number of columns of the matrix
    int columns;

    /// array of pointers to VectorOfFixVecs
    VectorOfFixVecs<FixVecType>** vec;
  };

  /// creates and inits a VectorOfFixVecs<DimVec<double> >
  VectorOfFixVecs<DimVec<double>>* createAndInit(int dim, int size, ...);

  /// creates and inits and double array
  double* createAndInitArray(int size, ...);

  template <class T>
  struct GradientType
  {
    using type = WorldVector<T>;

    static DenseVector<double> getValues(type const& t)
    {
      DenseVector<double> result(t.getSize());
      for (size_t i = 0; i < t.getSize(); i++)
        result[i] = static_cast<T>(t[i]);
      return result;
    }
  };

  template <class T>
  struct GradientType<WorldVector<T>>
  {
    using type = WorldVector<WorldVector<T>>;

    static DenseVector<double> getValues(type const& t)
    {
      DenseVector<double> result(sqr(t.getSize()));
      for (size_t i = 0; i < t.getSize(); i++)
        for (size_t j = 0; j < t.getSize(); j++)
          result[i*(t.getSize()) + j] = static_cast<T>(t[i][j]);
      return result;
    }
  };

  template <class T>
  using Gradient_t = typename GradientType<T>::type;

  template <class T>
  struct D2Type
  {
    using type = WorldMatrix<T>;

    static DenseVector<double> getValues(type const& t)
    {
      DenseVector<double> result(t.numRows() * t.numCols());
      for (size_t i = 0; i < t.numRows(); i++)
        for (size_t j = 0; j < t.numCols(); j++)
          result[i*(t.numCols()) + j] = static_cast<T>(t[i][j]);
      return result;
    }
  };

  template <class T>
  using D2_t = typename D2Type<T>::type;

  template <class T,GeoIndex d>
  inline void resize(FixVec<T,d>& /*vec*/, size_t /*newSize*/) {}

} // end namespace AMDiS


#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>

namespace mtl
{
  template <typename T>
  struct Collection<AMDiS::WorldVector<AMDiS::WorldVector<T>>>
  {
    typedef T      value_type;
    typedef int    size_type;
  };

  namespace ashape
  {
    template <typename T>
    struct ashape_aux<AMDiS::WorldVector<AMDiS::WorldVector<T>>>
    {
      typedef mat<typename ashape<T>::type> type;  // nonscal
    };
  } // end namespace ashape

  namespace traits
  {
    template <typename T>
    struct category<AMDiS::WorldVector<AMDiS::WorldVector<T>>>
    {
      typedef tag::dense2D type;
    };

    template <typename T>
    struct is_row_major<AMDiS::WorldVector<AMDiS::WorldVector<T>>>
      : boost::mpl::true_ {};

  } // end namespace traits

} // end namespace mtl
