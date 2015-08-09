/** \file ConstantTerms.hpp */

#pragma once

#include "BaseTerms.hpp"
#include "LazyOperatorTerm.h"

namespace AMDiS 
{
  /// Expression that encapsulates a runtime value
  template <class T>
  struct RTConstant : public BaseTerm<RTConstant<T>>,
                      public LazyOperatorTermBase
  {
    using Self = RTConstant;
    using value_type = Decay_t<T>;
    value_type value;
    
    constexpr RTConstant(value_type const& value_)
      : value{value_}
    {}

    value_type operator()(int) const { return value; }
    
    std::string str() const { return std::to_string(value); }
  };


  /// Expression that encapsulates a compiletime value
  template <int V>
  struct CTConstant : public BaseTerm<CTConstant<V>>,
                      public LazyOperatorTermBase
  {
    using Self = CTConstant;
    using value_type = int;
    static constexpr value_type value = V;
    
    constexpr CTConstant() {}
    
    value_type operator()(int) const { return value; }
    
    std::string str() const { return std::string("[") + std::to_string(value) + "]"; }
  };
  
  
  /// Expression that points to a value given by reference
  template <class T>
  struct Reference : public BaseTerm<Reference<T>>,
                     public LazyOperatorTermBase
  {
    using Self = Reference;
    using value_type = Decay_t<T>;
    const value_type& value;
    
    constexpr Reference(value_type const& value_) 
      : value{value_} {}

    value_type operator()(int) const { return value; }
    
    std::string str() const { return std::string("&(") + std::to_string(value) + ")"; }
  };

} // end namespace AMDiS
