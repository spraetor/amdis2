/** \file ConstantTerms.hpp */

#pragma once

#include <expressions/LazyOperatorTerm.h>
#include "BaseTerms.hpp"

namespace AMDiS
{
  /// Expression that encapsulates a runtime value
  template <class T>
  struct RTConstant : public BaseTerm<RTConstant<T>>,
        public LazyOperatorTermBase
  {
    using Self = RTConstant;
    using value_type = Decay_t<T>;

    constexpr RTConstant(value_type const& value_)
      : value{value_}
    {}

    value_type operator[](int) const
    {
      return value;
    }

    value_type operator()(WorldVector<double>) const
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
  template <int V>
  struct CTConstant : public BaseTerm<CTConstant<V>>,
        public LazyOperatorTermBase
  {
    using Self = CTConstant;
    using value_type = int;

    constexpr CTConstant() {}

    value_type operator[](int) const
    {
      return value;
    }

    value_type operator()(WorldVector<double>) const
    {
      return value;
    }

    std::string str() const
    {
      return std::string("[") + std::to_string(value) + "]";
    }

  private:
    constexpr static value_type value = V;
  };


  /// Expression that points to a value given by reference
  template <class T>
  struct Reference : public BaseTerm<Reference<T>>,
        public LazyOperatorTermBase
  {
    using Self = Reference;
    using value_type = Decay_t<T>;

    constexpr Reference(value_type const& value_)
      : value{value_} {}

    value_type operator[](int) const
    {
      return value;
    }

    value_type operator()(WorldVector<double>) const
    {
      return value;
    }

    std::string str() const
    {
      return std::string("&(") + std::to_string(value) + ")";
    }

  private:
    value_type const& value;
  };

} // end namespace AMDiS
