/** \file MatrixVectorBase.hpp */

#pragma once

#include <algorithm> // std::copy, std::fill
#include <ostream>   // std::basic_ostream

#include <boost/numeric/mtl/operation/assign_mode.hpp>

#include <Log.h>			// TEST_EXIT_DBG

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

#include <operations/meta.hpp>
#include <operations/assign.hpp>

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
    using model = Model;

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
      this->operator=(expr);
    }

    /// destructor
    ~MatrixVectorBase() {}

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
    model& operator=(BaseExpr<Expr> const& expr)
    {
      assign(static_cast<Expr const&>(expr), mtl::assign::assign_sum());
      return static_cast<model&>(*this);
    }

    /// compound plus-assignment of an expression
    template <class Expr>
    model& operator+=(BaseExpr<Expr> const& expr)
    {
      assign(static_cast<Expr const&>(expr), mtl::assign::plus_sum());
      return static_cast<model&>(*this);
    }

    /// compound minus-assignment of an expression
    template <class Expr>
    model& operator-=(BaseExpr<Expr> const& expr)
    {
      assign(static_cast<Expr const&>(expr), mtl::assign::minus_sum());
      return static_cast<model&>(*this);
    }

    /// compound times-assignment of an expression
    template <class Expr>
    model& operator*=(BaseExpr<Expr> const& expr)
    {
      assign(static_cast<Expr const&>(expr), mtl::assign::times_sum());
      return static_cast<model&>(*this);
    }

    /// compound divides-assignment of an expression
    template <class Expr>
    model& operator/=(BaseExpr<Expr> const& expr)
    {
      assign(static_cast<Expr const&>(expr), mtl::assign::divide_sum());
      return static_cast<model&>(*this);
    }

    // ---------------------------------------------------------------------------
  public:
    /// Assignment operator for scalars
    template <class S>
    Requires_t<std::is_convertible<S, value_type>, model> &
    operator=(S value)
    {
      for_each(assign::value<value_type, S>(value));
      return static_cast<model&>(*this);
    }

    /// compound assignment *= of a scalar
    template <class S>
    Requires_t<std::is_convertible<S, value_type>, model> &
    operator*=(S value)
    {
      for_each(assign::mult_value<value_type, S>(value));
      return static_cast<model&>(*this);
    }

    /// compound assignment /= of a scalar
    template <class S>
    Requires_t<std::is_convertible<S, value_type>, model> &
    operator/=(S value)
    {
      for_each(assign::div_value<value_type, S>(value));
      return static_cast<model&>(*this);
    }

    // ---------------------------------------------------------------------------
  private:

    /// basic assignment for compound operators given by the Assigner
    template <class Expr, class Assigner>
    void assign(Expr const& expr, Assigner assigner)
    {
      TEST_EXIT_DBG( _size == size(expr) )("Sizes do not match!\n");
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
    static constexpr size_t eval(size_t s)
    {
      return s;
    }
  };

  template <size_t S>
  struct StaticSizePolicy
  {
    static constexpr size_t value = S;

    /// return static size parameter \p S
    static constexpr size_t eval(size_t)
    {
      return S;
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
