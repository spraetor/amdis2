/** \file FixVec.h */

#pragma once

#include <iostream>

#include "AMDiS_fwd.h"
#include "AMDiS_base.h"
#include "Global.h"
#include "Math.h"
#include "MatrixVector.h"

namespace AMDiS 
{

  /// determines how to initialize a FixVec when constructed.
  enum InitType {
    NO_INIT = 0,       /**< no initialisation */
    VALUE_LIST = 1,    /**< a complete value list is given */
    DEFAULT_VALUE = 2, /**< all values ar set to a given default value */
    UNIT_VECTOR = 3,   /**< the i-th value is 1, all other values are 0 */
    UNIT_MATRIX = 4    /**< values at the diagonal of a matrix are set to one */
  };

  /** \ingroup Common
   * \brief
   * A FixVec is a template vector of a fixed size. 
   *
   * The size is determined at construction time and depends on the dimension
   * and the template parameter d. So a FixVec<int, VERTEX> is a integer vector
   * with 3 entries in 2d and with 4 entries in 3d. The dimension and the way
   * the vector should be initialized are specified by the constructor call.
   */
  template <class T, GeoIndex d>
  class FixVec : public Vector<T>
  {
    typedef FixVec     self;
    typedef Vector<T>  super;
  public:
    
    /// Constructor without initialisation. initType must be NO_INIT. If dim is
    /// not spezified, a FixVec for DIM_OF_WORLD is created.
    FixVec(int dim = -1, InitType initType = NO_INIT)
      : Vector<T>(calcSize(dim))
    { 
      TEST_EXIT_DBG(initType == NO_INIT)("wrong initType or missing initializer\n");
    }

    /// constructor with value list initialisation. initType must be VALUE_LIST.
    /// ini is an array which contains the initialisation values.
    FixVec(int dim, InitType initType, const T* ini)
      : Vector<T>(calcSize(dim))
    {
      TEST_EXIT_DBG(initType == VALUE_LIST)("wrong initType or wrong initializer\n");
      this->setValues(ini);
    }

    /// constructor with default value initialisation. initType must be
    /// DEFAULT_VALUE. All vector entries are set to ini.
    FixVec(int dim, InitType initType, const T& ini)
      : Vector<T>(calcSize(dim))
    {
      TEST_EXIT_DBG(initType == DEFAULT_VALUE)("wrong initType or wrong initializer\n");
      this->set(ini);
    }
    
    FixVec(self const& other)
      : super(other)
    { }

    /// Initialisation for dim.
    inline void init(int dim) 
    {
      this->resize(calcSize(dim));
    }

    /// Initialisation for size
    inline void initSize(int size) 
    {
      this->resize(size);
    }

    /// Returns the \ref size_ of the FixVec.
//     inline int size() const 
//     { 
//       return this->getSize(); 
//     } 

  protected:
    /// Determines needed vector size.
    static int calcSize(int dim) 
    {
      if (dim < 0)
        return Global::getGeo(WORLD);
      else
        return Global::getGeo(d, dim);
    }
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
    VectorOfFixVecs(int d, int s, InitType initType) 
      : size(s),
        dim(d)
    {
      TEST_EXIT_DBG(initType == NO_INIT)("wrong initType or wrong initializer\n");

      vec.resize(size);
      for (int i = 0; i < size; i++)
        vec[i] = new FixVecType(dim, NO_INIT);
    }

    /// constructs a VectorOfFixVecs via an value list.  dim is passed to 
    /// FixVec's constructors. size_ is the number of contained FixVecs. initType
    /// must be VALUE_LIST. ini contains the initialisation values.
    VectorOfFixVecs(int d, int s, InitType initType, FixVecType const* ini)
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
    VectorOfFixVecs(int d, int s, InitType initType, const FixVecType& ini)
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
      if (this != &rhs) {
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
      	vec[i] = new FixVecType(dim, NO_INIT);
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
    MatrixOfFixVecs(int dim, int r, int c, InitType initType)
      : rows(r), 
      	columns(c)
    {
      TEST_EXIT_DBG(initType == NO_INIT)("wrong initType or wrong initializer\n");
      vec = new VectorOfFixVecs<FixVecType>*[rows];
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
    inline VectorOfFixVecs<FixVecType >& operator[](int index)
    {
      TEST_EXIT_DBG(index >= 0 && index < rows)("invalid index\n");
      return *(vec[index]);
    }

    /// assignment operator
    inline const VectorOfFixVecs<FixVecType >& operator[](int index) const
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
    VectorOfFixVecs<FixVecType> **vec;
  };


  /** \ingroup Common
   * \brief
   * A DimVec is a FixVec with dim + 1 entries. It can be used for storing
   * barycentric coordinates or information associated to vertices or
   * parts of an element.
   */
  template<typename T>
  class DimVec : public FixVec<T,PARTS> 
  {
  public:
    DimVec() {}

    /// Calls the corresponding constructor of FixVec
    DimVec(int dim, InitType initType = NO_INIT)
      : FixVec<T,PARTS>(dim, initType)
    {}

    /// Calls the corresponding constructor of FixVec
    DimVec(int dim, InitType initType, T const* ini)
      : FixVec<T,PARTS>(dim, initType, ini)
    {}

    /// Calls the corresponding constructor of FixVec
    DimVec(int dim, InitType initType, const T& ini)
      : FixVec<T,PARTS>(dim, initType, ini)
    {}

    /// Multplication of a matrix with a vector.
    void multMatrixVec(DimMat<T> &m, DimVec<T> &v);
  };


  /** \ingroup Common
   * \brief
   * A DimMat is a VectorOfFixVecs of dim+1 DimVecs. 
   */
  template<typename T>
  class DimMat : public Matrix<T>
  {
  public:
    /// Calls the corresponding constructor of VectorOfFixVecs
    DimMat(int dim, InitType initType = NO_INIT)
      : Matrix<T>(dim + 1, dim + 1)
    {}

    /// Calls the corresponding constructor of VectorOfFixVecs
    DimMat(int dim, InitType initType, const T& ini)
      : Matrix<T>(dim + 1, dim + 1)
    {
      TEST_EXIT_DBG(initType == DEFAULT_VALUE)
        ("wrong initType or wrong initializer\n");    
      this->set(ini);
    }

    /// Calls the corresponding constructor of VectorOfFixVecs
    DimMat(int dim, InitType initType, T const* ini)
      : Matrix<T>(dim + 1, dim + 1)
    {
      TEST_EXIT_DBG(initType == VALUE_LIST)("wrong initType or wrong initializer\n");
      setValues(ini);
    }
  };


  /** \ingroup Common
   * \brief
   * A WorldVector is an AlgoVec with DIM_OF_WORLD entries of type double.
   * Can be used for storing world coordinates.
   */
  template <typename T>
  class WorldVector : public FixVec<T, WORLD>
  {
  public:
    typedef WorldVector       self;
    typedef FixVec<T, WORLD>  super;
  
    /// Calls the corresponding constructor of AlgoVec
    WorldVector() 
      : super(Global::getGeo(WORLD), NO_INIT) 
    {}

    /// Calls the corresponding constructor of AlgoVec
    WorldVector(InitType initType, T const* ini) 
      : super(Global::getGeo(WORLD), initType, ini)
    {}

    /// Calls the corresponding constructor of AlgoVec
    WorldVector(InitType initType, const T& ini)
      : super(Global::getGeo(WORLD), initType, ini)
    {}
//     
    /// Copy constructor for other of same value_type
    WorldVector(self const& other)
      : super(other) 
    {}
    
    /// Copy assignement operator
//     self& operator=(self const& other)
//     {
//       assert( Global::getGeo(WORLD) == other.getSize() );
//       this->setValues(other.getValArray());
//       return *this;
//     }
      
    
    /// Assignement operator
    template <typename S>
    self& operator=(const Vector<S>& other) 
    {
      TEST_EXIT_DBG( other.getSize() == this->size )
		   ("Wrong dimensions in assignment.\n");
      this->setValues(other.getValArray());
      return *this;
    }

    /// Sets all entries to scal
    template <typename S>
    typename enable_if<boost::is_convertible<S, T>, self>::type &
    operator=(S scal)
    {
      this->set(scal);
      return *this;
    }

    /// Sets the arrays value to the geometric midpoint of the points  p1 and p2.
    inline void setMidpoint(const WorldVector<T> &p1, const WorldVector<T> &p2)
    {
      int dow = Global::getGeo(WORLD);

      for (int i = 0; i < dow; i++)
	this->valArray[i] = 0.5 * (p1[i] + p2[i]);
    }

    /// Multplication of a matrix with a vector.
    void multMatrixVec(const WorldMatrix<T> &m, const WorldVector<T> &v);
  };


  /** \ingroup Common
   * \brief
   * A WorldMatrix is a FixVec with DIM_OF_WORLD WorldVector entries.
   * So it can be seen as matrix with DIM_OF_WORLD x DIM_OF_WORLD double 
   * entries.
   * Here it's not necessary to use the VectorOfFixVecs class, because the 
   * FixVec constructor assumes dim = DIM_OF_WORLD, if no dim is spezified, 
   * what is the right assumption in this case.
   */
  template<typename T>
  class WorldMatrix : public Matrix<T>
  {
  public:
    typedef WorldMatrix   self;
    typedef Matrix<T>     super;
    
    /// Calls the corresponding constructor of FixVec
    WorldMatrix()
      : super(Global::getGeo(WORLD), Global::getGeo(WORLD))
    {}

    /// Calls the corresponding constructor of FixVec
    WorldMatrix(InitType initType, T* ini)
      : super(Global::getGeo(WORLD), Global::getGeo(WORLD))
    {
      TEST_EXIT_DBG(initType == VALUE_LIST)("???\n");
      this->setValues(ini);
    }

    /** \brief
     * Calls the corresponding constructor of FixVec and sets all matrix entries
     * to ini
     */
    WorldMatrix(InitType initType, const T& ini)
      : super(Global::getGeo(WORLD), Global::getGeo(WORLD))
    {
      TEST_EXIT_DBG(initType == DEFAULT_VALUE)("wrong initType or wrong initializer\n");
      this->set(ini);
    }
    
    /// Copy constructor 
    WorldMatrix(self const& other)
      : super(other) 
    { }
    
    /// Copy assignment operator
//     self& operator=(self other) 
//     {
//       swap(*this, other);
//       return *this;
//     }
    
    /// Assignment operator
    template <typename S>
    self& operator=(const Matrix<S>& other) 
    {
      TEST_EXIT_DBG( this->size == other.getSize() && 
		     this->rows == other.getNumRows() && 
		     this->cols == other.getNumCols() )
		   ("Wrong dimensions in assignment.\n");
      this->setValues(other.getValArray());
      return *this;
    }
    
    /// Assignment operator for scalars
    template <typename S>
    typename enable_if<boost::is_convertible<S, T>, self>::type &
    operator=(S value) 
    {
      this->set(value);
      return *this;
    }
  
    /// Returns true if the matrix is a diagonal matrix, returns false otherwise.
    bool isDiagMatrix() const;

    /// Returns true if the matrix is symmetric, returns false otherwise.
    bool isSymmetric() const;

    /// Sets the diagonal entries to the given value.
    void setDiag(T value);

    /** \brief
     * Creates a matrix from the "multiplication" of two vectors.
     *
     * 2d-Example:
     *
     * /a\   /c\   /ac ad\
     * | | * | | = |     |
     * \b/   \d/   \bc bd/
     */
    void vecProduct(const WorldVector<T>& v1, const WorldVector<T>& v2);
  };

  /// FixVec operator for stream output
  template<typename T, GeoIndex d>
  std::ostream& operator <<(std::ostream& out, const FixVec<T,d>& fixvec)
  {
    for (int i = 0; i < fixvec.getSize() - 1; i++)
      out << fixvec[i] << " ";
    out << fixvec[fixvec.getSize()-1];
    return out;
  }

  /// creates and inits a VectorOfFixVecs<DimVec<double> >
  VectorOfFixVecs<DimVec<double> > *createAndInit(int dim, int size, ...);

  /// creates and inits and double array
  double *createAndInitArray(int size, ...); 

  template <class T>
  struct GradientType
  {
    typedef WorldVector<T> type;
    
    static DenseVector<double> getValues(type &t) {
      DenseVector<double> result(t.getSize());
      for (size_t i = 0; i < t.getSize(); i++)
	result[i] = static_cast<T>(t[i]);
      return result;
    }
  };

  template <class T>
  struct GradientType<WorldVector<T> >
  {
    typedef WorldVector<WorldVector<T> > type;

    static DenseVector<double> getValues(type &t) {
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
    typedef WorldMatrix<T> type;
    
    static DenseVector<double> getValues(type &t) {
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
  inline void resize(FixVec<T,d>& vec, size_t newSize)
  { }
  
} // end namespace AMDiS


#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>

namespace mtl 
{
  template <typename T>
  struct Collection<AMDiS::WorldVector<AMDiS::WorldVector<T> > >
  {
    typedef T      value_type;
    typedef int    size_type;
  };

  namespace ashape 
  {
    template <typename T> 
    struct ashape_aux<AMDiS::WorldVector<AMDiS::WorldVector<T> > > 
    { 
      typedef mat<typename ashape<T>::type> type;  // nonscal
    };
  } // end namespace ashape  
  
  namespace traits
  {
    template <typename T>
    struct category<AMDiS::WorldVector<AMDiS::WorldVector<T> > > 
    {
	typedef tag::dense2D type;
    };
    
    template <typename T>
    struct is_row_major<AMDiS::WorldVector<AMDiS::WorldVector<T> > >
      : boost::mpl::true_ {};
      
  } // end namespace traits
  
} // end namespace mtl

#include "FixVec.hh"
