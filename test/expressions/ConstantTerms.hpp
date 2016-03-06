#pragma once

#include <memory> // std::addressof

#include "expressions/BaseTerms.hpp"
#include "expressions/LazyOperatorTerm.hpp"

namespace AMDiS
{
  /// Expression that encapsulates a runtime value
  template <class T>
  struct RTConstant
    : public ShapedTerm_t<Decay_t<T>, RTConstant<T>>,
      public LazyOperatorTermBase
  {
    using Self = RTConstant;
    using value_type = Decay_t<T>;

    RTConstant(value_type const& value_)
      : value(value_)
    {
      MSG("RTConstant()");
    }

//     RTConstant(Self const&) = default;
    RTConstant(Self&&) = default;
    Self& operator=(Self&&) = default;

    value_type operator()(WorldVector<double> /* x */) const
    {
      return value;
    }

    std::string str() const
    {
      return std::to_string(value);
    }

  private:
    value_type value;
  };


  /// Expression that encapsulates a compiletime value
  template <long V>
  struct CTConstant
    : public BaseTerm<CTConstant<V>>,
      public LazyOperatorTermBase
  {
    using Self = CTConstant;
    using value_type = long;

    constexpr CTConstant() {}

    value_type operator()(WorldVector<double> /* x */) const
    {
      return value;
    }

    std::string str() const
    {
      return std::string("[") + std::to_string(value) + "]";
    }

  private:
    static constexpr value_type value = V;
  };


  /// Expression that points to a value given by reference
  template <class T>
  struct Reference
    : public ShapedTerm_t<Decay_t<T>, Reference<T>>,
      public LazyOperatorTermBase
  {
    using Self = Reference;
    using value_type = Decay_t<T>;

    // construct/copy/destroy
    Reference(value_type const& ref) noexcept
      : ptr(std::addressof(ref)) {}

    Reference(value_type&&) = delete;
    Reference(Self const&) noexcept = default;

    value_type const& operator()(WorldVector<double> /* x */) const
    {
      return *ptr;
    }

    std::string str() const
    {
      return std::string("&(") + std::to_string(*ptr) + ")";
    }

  private:
    value_type const* ptr;
  };

} // end namespace AMDiS
