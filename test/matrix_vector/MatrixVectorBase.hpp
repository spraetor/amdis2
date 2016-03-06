#pragma once

// std c++ headers
#include <algorithm> // std::copy, std::fill
#include <ostream>   // std::basic_ostream

// MTL4 includes
// #include <boost/numeric/mtl/operation/assign_mode.hpp>

// AMDiS includes
#include "Log.hpp"			// TEST_EXIT_DBG
#include "matrix_vector/expr/base_expr.hpp"     // BaseExpr
#include "operations/assign.hpp"
#include "operations/meta.hpp"
#include "traits/basic.hpp"
#include "traits/traits_fwd.hpp"

#define DEFAULT_SIZE 0

namespace AMDiS
{
  /// Base class for matrix and vector types.
  /**
    * This class manages the assignment operators and element access operators.
    */
  template <class Model, class MemoryPolicy>
  struct MatrixVectorBase
    : public MemoryPolicy
  {
    using Self  = MatrixVectorBase;
    using Super = MemoryPolicy;

    using value_type = Value_t<Super>;
    using size_type  = Size_t<Super>;

    using pointer        = value_type*;
    using const_pointer  = value_type const*;
    using iterator       = pointer;
    using const_iterator = const_pointer;

    using Super::_elements;
    using Super::_size;

    // ---------------------------------------------------------------------------
  protected:
    /// default constructor
    explicit MatrixVectorBase(size_type s = 0)
      : Super(s)
    {}

  public:
    /// assignment of an expression
    template <class Expr>
    explicit MatrixVectorBase(BaseExpr<Expr> const& expr)
      : Super(size(expr))
    {
      MSG("MatrixVectorBase->Copy(BaseExpr(" << &expr.sub() << "))");
      this->operator=(expr.sub());
    }

    // use default implementations for copy and move operations
    MatrixVectorBase(Self const&) = default;
    MatrixVectorBase(Self&&)      = default;

    Self& operator=(Self const&)  = default;
    Self& operator=(Self&&)       = default;

    /// fill vector with scalar value
    template <class S>
    Requires_t<std::is_convertible<S, value_type>>
    set(S const& value)
    {
      std::fill(_elements, _elements + _size, value);
    }

    /// fill vector with values from pointer. No check of length is performed!
    void setValues(value_type* values)
    {
      std::copy(values, values + _size, _elements);
    }

    /// Access to the i-th data element.
    value_type& operator()(size_type i)
    {
      return _elements[i];
    }

    /// Access to the i-th data element. (const variant)
    value_type const& operator()(size_type i) const
    {
      return _elements[i];
    }

    /// Returns pointer to the first vector element.
    iterator begin()
    {
      return _elements;
    }

    /// Returns pointer to the first vector element. (const variant)
    const_iterator begin() const
    {
      return _elements;
    }

    /// Returns pointer after the last vector element.
    iterator end()
    {
      return _elements + _size;
    }

    /// Returns pointer after the last vector element. (const variant)
    const_iterator end() const
    {
      return _elements + _size;
    }

    // ---------------------------------------------------------------------------
  public:
    /// assignment of an expression
    template <class Expr>
    Model& operator=(BaseExpr<Expr> const& expr)
    {
      MSG("MatrixVectorBase->Assign(BaseExpr(" << &expr.sub() << "))");
      assign(expr.sub(), assign::assign<value_type>());
      return static_cast<Model&>(*this);
    }

    /// compound plus-assignment of an expression
    template <class Expr>
    Model& operator+=(BaseExpr<Expr> const& expr)
    {
      assign(expr.sub(), assign::plus<value_type>());
      return static_cast<Model&>(*this);
    }

    /// compound minus-assignment of an expression
    template <class Expr>
    Model& operator-=(BaseExpr<Expr> const& expr)
    {
      assign(expr.sub(), assign::minus<value_type>());
      return static_cast<Model&>(*this);
    }

    /// compound times-assignment of an expression
    template <class Expr>
    Model& operator*=(BaseExpr<Expr> const& expr)
    {
      assign(expr.sub(), assign::multiplies<value_type>());
      return static_cast<Model&>(*this);
    }

    /// compound divides-assignment of an expression
    template <class Expr>
    Model& operator/=(BaseExpr<Expr> const& expr)
    {
      assign(expr.sub(), assign::divides<value_type>());
      return static_cast<Model&>(*this);
    }

    // ---------------------------------------------------------------------------
  public:
    /// Assignment operator for scalars
    template <class S>
    Requires_t<std::is_convertible<S, value_type>, Model> &
    operator=(S value)
    {
      for_each(assign::value<value_type, S>(value));
      return static_cast<Model&>(*this);
    }

    /// compound assignment *= of a scalar
    template <class S>
    Requires_t<std::is_convertible<S, value_type>, Model> &
    operator*=(S value)
    {
      for_each(assign::mult_value<value_type, S>(value));
      return static_cast<Model&>(*this);
    }

    /// compound assignment /= of a scalar
    template <class S>
    Requires_t<std::is_convertible<S, value_type>, Model> &
    operator/=(S value)
    {
      for_each(assign::div_value<value_type, S>(value));
      return static_cast<Model&>(*this);
    }

    // ---------------------------------------------------------------------------
  private:

    /// basic assignment for compound operators given by the Assigner
    template <class Expr, class Assigner>
    void assign(Expr const& expr, Assigner assigner)
    {
      MSG("assign()(Expr(" << &expr << "))");
#ifndef NDEBUG
      size_t s = size(expr);
      TEST_EXIT( _size == s,
                 "Sizes do not match! _size = " << _size << ", size(expr) = " << s << "\n");
#endif
      Super::assign_aux(static_cast<Model&>(*this), expr, assigner);
    }

    /// basic assignment for compound operators given by the Assigner
    template <class Functor>
    inline void for_each(Functor&& f)
    {
      Super::for_each_aux(std::forward<Functor>(f));
    }
  };


  // ===========================================================================

  struct DefaultSizePolicy
  {
    /// return argument \param s
    template <class size_type>
    static constexpr size_type eval(size_type s)
    {
      return s;
    }
  };

  template <size_t S>
  struct StaticSizePolicy
  {
    static constexpr size_t value = S;

    /// return static size parameter \p S
    template <class size_type>
    static constexpr size_type eval(size_type)
    {
      return size_type(S);
    }
  };

  // ===========================================================================

  // OUTPUT << MatrixVectorBase
  template <class charT, class traits, class Model, class Memory>
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT,traits>& out,
             MatrixVectorBase<Model, Memory> const& vector)
  {
    auto it = vector.begin();
    out << *it;
    ++it;
    for (; it != vector.end(); ++it)
      out << ' ' << *it;
    return out;
  }

} // end namespace AMDiS
