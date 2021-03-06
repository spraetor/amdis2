#pragma once

// AMDiS headers
#include "Config.hpp"
#include "Log.hpp"			// TEST_EXIT_DBG
#include "operations/generic_loops.hpp"	// meta::FOR
#include "utility/aligned_alloc.hpp"	// ALIGNED_ALLOC, ALIGNED_FREE, ...

namespace AMDiS
{
  /// Memory base for vector types using static storage
  /** The template parameter \p T describes the value-type of the
   *  data elements. The maximal size (capacity) of the container
   *  ist set by the non-type template parameter \p N. In the case
   *  that the container is a matrix, a second dimension can be
   *  specified by the non-type template parameter \p M. Then the
   *  total size is N*M.
   **/
  template <class T, small_t N, small_t M = 1>
  struct MemoryBaseStatic
  {
    static_assert( N*M > 0 , "Container size must be > 0" );

    using value_type    = T;
    using size_type     = small_t;
    using pointer       = value_type*;
    using const_pointer = value_type const*;

    // static sizes
    static constexpr int _SIZE = N*M;
    static constexpr int _ROWS = N;
    static constexpr int _COLS = M;

  protected:
    static constexpr size_type _size = _SIZE;
    static constexpr size_type _capacity = _SIZE; // TODO: eventuell aufrunden

    AMDIS_ALIGNED(T, _elements, _capacity);   // T _elements[N];

  protected:
    /// default constructor
    explicit MemoryBaseStatic(size_type DBG_VAR( s ) = 0)
    {
      TEST_EXIT_DBG(s == 0 || s == _SIZE)("Size must be equal to capacity!\n");
    }

    /// destructor
    ~MemoryBaseStatic() = default;

  public:
    // use default implementations for copy and move operations
    MemoryBaseStatic(MemoryBaseStatic const&)            = default;
    MemoryBaseStatic(MemoryBaseStatic&&)                 = default;
    MemoryBaseStatic& operator=(MemoryBaseStatic const&) = default;
    MemoryBaseStatic& operator=(MemoryBaseStatic&&)      = default;

  public:
    /// return the \ref _size of the vector.
    size_type getSize() const
    {
      return _size;
    }

    /// return the \ref _capacity of the vector.
    static constexpr size_type getCapacity()
    {
      return _capacity;
    }

    /// return the amount of memory in Bytes allocated by this vector.
    size_type getMemoryUsage() const
    {
      return _capacity*sizeof(T) + sizeof(size_type);
    }

    /// return address of contiguous memory block \ref _elements
    pointer data()
    {
      return _elements;
    }

    /// return address of contiguous memory block \ref _elements (const version)
    const_pointer data() const
    {
      return _elements;
    }

    /// resize the vector. Only possible, if \p s <= _capacity
    void resize(size_type s)
    {
      if (s != _size)
      {
        // Not allowed
        assert( false );
      }
    }

  protected:

    template <class Target, class Source, class Assigner>
    void assign_aux(Target& target, Source const& src, Assigner assigner)
    {
      using meta::FOR;
      FOR<0,_SIZE>::assign(target, src, assigner);
      //       for (size_type i = 0; i < _size; ++i)
      // 	Assigner::apply(target(i), src(i));
    }

    template <class Functor>
    void for_each_aux(Functor f)
    {
      using meta::FOR;
      FOR<0,_SIZE>::for_each(_elements, f);
      //       for (size_type i = 0; i < _size; ++i)
      // 	f(_elements[i]);
    }
  };

  // ===========================================================================

  /// Memory base for vector types using dynamic storage
  /** The template parameter \p T describes the value-type of the
   *  data elements. The memory is allocated on the heap. When
   *  \p aligned is set to true an 16-Byte alignement of the data
   *  is enforced, in order to use vectorization methods.
   **/
  template <class T, bool aligned>
  struct MemoryBaseDynamic
  {
    using Self = MemoryBaseDynamic;

    using value_type    = T;
    using size_type     = small_t;
    using pointer       = value_type*;
    using const_pointer = value_type const*;

    // static sizes (by default -1 := dynamic size)
    static constexpr int _SIZE = -1;
    static constexpr int _ROWS = -1;
    static constexpr int _COLS = -1;

  protected:
    size_type  _size;
    size_type  _capacity;
    T*         _elements;

  protected:
    /// default constructor
    explicit MemoryBaseDynamic(size_type s = 0)
      : _size(s),
        _capacity(s),
        _elements(_size ? (aligned ? AMDIS_ALIGNED_ALLOC(T, s) : new T[s]) : NULL)
    {}

    /// destructor
    ~MemoryBaseDynamic()
    {
      if (_elements)
      {
        if (aligned)
        {
          AMDIS_ALIGNED_FREE(_elements);
        }
        else
        {
          delete [] _elements;
        }
        _elements = NULL;
      }
    }

  public:

    /// copy constructor
    MemoryBaseDynamic(Self const& other)
      : _size(other._size),
        _capacity(other._capacity),
        _elements(_size ? (aligned ? AMDIS_ALIGNED_ALLOC(T, _size) : new T[_size]) : NULL)
    {
      std::copy(other._elements, other._elements + _size, _elements);
    }

    /// copy assignment operator
    MemoryBaseDynamic& operator=(Self const& other)
    {
      if (this != &other) {
        resize(other._size);
        std::copy(other._elements, other._elements + other._size, _elements);
      }
      return *this;
    }

    /// move constructor
    MemoryBaseDynamic(Self&& other)
      : _size(other._size),
        _capacity(other._capacity),
        _elements(other._elements)
    {
      other._elements = NULL;
      other._size = 0;
    }

    // move assignment operator
    MemoryBaseDynamic& operator=(Self&& other)
    {
      if (this != &other)
      {
        if (_elements)
        {
          if (aligned)
          {
            AMDIS_ALIGNED_FREE(_elements);
          }
          else
          {
            delete [] _elements;
          }
          _elements = NULL;
        }

        _size = other._size;
        _capacity = other._capacity;
        _elements = other._elements;

        other._elements = NULL;
        other._size = 0;
      }
      return *this;
    }

  public:
    /// return the \ref _size of the vector.
    size_type getSize() const
    {
      return _size;
    }

    /// return the \ref _capacity of the vector.
    size_type getCapacity() const
    {
      return _capacity;
    }

    /// return the amount of memory in Bytes allocated by this vector.
    size_type getMemoryUsage() const
    {
      return (aligned ? AMDIS_ALIGNED_SIZE(_capacity*sizeof(T), CACHE_LINE)
                      : _capacity*sizeof(T))
             + 2*sizeof(size_type);
    }

    /// return address of contiguous memory block \ref _elements
    pointer data()
    {
      return _elements;
    }

    /// return address of contiguous memory block \ref _elements (const version)
    const_pointer data() const
    {
      return _elements;
    }

    /// resize the vector. If \p s <= \ref _capacity simply set the \ref _size
    /// attribute the s, otherwise realloc memory
    void resize(size_type s)
    {
      if (s <= _capacity)
      {
        _size = s;
        // TODO: (maybe) set unused entries to NULL
      }
      else
      {
        if (_elements)
        {
          if (aligned)
          {
            AMDIS_ALIGNED_FREE(_elements);
          }
          else
          {
            delete [] _elements;
          }
        }
        _elements = aligned ? AMDIS_ALIGNED_ALLOC(T, s) : new T[s];
        _size = s;
      }
    }

  protected:
    template <class Target, class Source, class Assigner>
    void assign_aux(Target& target, Source const& src, Assigner assigner)
    {
      assign_aux(target, src, assigner, bool_<aligned>());
    }

    template <class Target, class Source, class Assigner> // not assume aligned
    void assign_aux(Target& target, Source const& src, Assigner /*assigner*/, false_)
    {
      for (size_type i = 0; i < _size; ++i)
        Assigner::apply(target(i), src(i));
    }

    template <class Target, class Source, class Assigner> // assume aligned
    void assign_aux(Target& /*target*/, Source const& src, Assigner /*assigner*/, true_)
    {
      value_type* var = (value_type*)AMDIS_ASSUME_ALIGNED(_elements);
      for (size_type i = 0; i < _size; ++i)
        Assigner::apply(var[i], src(i));
    }

    template <class Functor>
    void for_each_aux(Functor f)
    {
      for (size_type i = 0; i < _size; ++i)
        f(_elements[i]);
    }
  };

  // ===========================================================================

  /// Memory base for vector types using static storage, with dynamic size
  /** The template parameter \p T describes the value-type of the
   *  data elements. The maximal size (capacity) of the container
   *  ist set by the non-type template parameter \p N. In the case
   *  that the container is a matrix, a second dimension can be
   *  specified by the non-type template parameter \p M. Then the
   *  total size allocated is N*M. The internal size used is set in the
   *  constructor or the \ref resize function.
   **/
  template <class T, small_t N, small_t M=1>
  struct MemoryBaseHybrid
  {
    static_assert( N*M > 0 , "Container size must be > 0" );

    using value_type    = T;
    using size_type     = small_t;
    using pointer       = value_type*;
    using const_pointer = value_type const*;

    // static sizes
    static constexpr int _SIZE = -1;
    static constexpr int _ROWS = -1;
    static constexpr int _COLS = -1;

  protected:
    size_type _size;
    static constexpr size_type _capacity = N*M;

    AMDIS_ALIGNED(T, _elements, _capacity);   // T _elements[N];

  protected:
    /// default constructor
    explicit MemoryBaseHybrid(size_type s = 0)
      : _size(s)
    {
      TEST_EXIT_DBG(s <= _capacity)("Size must be <= capacity!\n");
    }

    /// destructor
    ~MemoryBaseHybrid() = default;

  public:
    // use default implementations for copy and move operations
    MemoryBaseHybrid(MemoryBaseHybrid const&)            = default;
    MemoryBaseHybrid(MemoryBaseHybrid&&)                 = default;
    MemoryBaseHybrid& operator=(MemoryBaseHybrid const&) = default;
    MemoryBaseHybrid& operator=(MemoryBaseHybrid&&)      = default;

  public:
    /// return the \ref _size of the vector.
    size_type getSize() const
    {
      return _size;
    }

    /// return the \ref _capacity of the vector.
    size_type getCapacity() const
    {
      return _capacity;
    }

    /// return the amount of memory in Bytes allocated by this vector.
    size_type getMemoryUsage() const
    {
      return _capacity*sizeof(T) + sizeof(size_type);
    }

    /// return address of contiguous memory block \ref _elements
    pointer data()
    {
      return _elements;
    }

    /// return address of contiguous memory block \ref _elements (const version)
    const_pointer data() const
    {
      return _elements;
    }

    /// resize the vector. If \p s <= \ref _capacity simply set the \ref _size
    /// attribute the s, otherwise realloc memory
    void resize(size_type s)
    {
      if (s <= _capacity)
      {
        _size = s;
      }
      else
      {
        // not supported
        assert( false );
      }
    }

  protected:
    template <class Target, class Source, class Assigner>
    void assign_aux(Target& target, Source const& src, Assigner /*assigner*/)
    {
      for (size_type i = 0; i < _size; ++i)
        Assigner::apply(target(i), src(i));
    }

    template <class Functor>
    void for_each_aux(Functor f)
    {
      for (size_type i = 0; i < _size; ++i)
        f(_elements[i]);
    }
  };

} // end namespace AMDiS
