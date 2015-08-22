/** \file TermGenerator.hpp */

#pragma once

#include <traits/basic.hpp>

#include "BaseTerms.hpp"
#include "ConstantTerms.hpp"
#include "CoordsTerm.hpp"
#include "FunctorTerm.hpp"
#include "TermConcepts.hpp"

namespace AMDiS
{
  namespace detail
  {
    template <class T, class = void>
    struct ToTerm 
    {
      using type = RTConstant<T>;
      
      template <class T_>
      constexpr static type get(T_&& t) { return {t}; }
    };
    
    template <class T>
    struct ToTerm<T, Requires_t<concepts::Term<T>> >
    {
      using type = T;
      
      template <class T_>
      constexpr static T_&& get(T_&& t) { return std::forward<T_>(t); }
    };
    
    template <class F>
    struct ToTerm<F, Requires_t<concepts::CoordsFunctor<F>> >
    {
      using type = FunctorTerm<F, CoordsTerm<>>;
      
      template <class F_>
      static type get(F_&& f) { return {std::forward<F_>(f), CoordsTerm<>()}; }
    };
    
  } // end namespace detail
  
  template <class T>
  using ToTerm = detail::ToTerm<T>;
  
  template <class T>
  using ToTerm_t = typename ToTerm<T>::type;
  
  template <class T>
  inline auto toTerm(T&& t) RETURNS( ToTerm<T>::get(std::forward<T>(t)) )

} // end namespace AMDiS
