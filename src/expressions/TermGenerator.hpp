/** \file TermGenerator.hpp */

#pragma once

#include "BaseTerms.hpp"
#include "ConstantTerms.hpp"
#include "CoordsTerm.hpp"
#include "FunctorTerm.hpp"
#include <traits/basic.hpp>
#include <traits/concepts.hpp>

namespace AMDiS
{
  namespace detail
  {
    template <class T, class = void>
    struct ToTerm 
    {
      using type = RTConstant<T>;
      
      template <class T_>
      static constexpr type get(T_&& t) { return {t}; }
    };
    
    template <class T>
    struct ToTerm<T, typename enable_if< concepts::Term<T> >::type>
    {
      using type = T;
      
      template <class T_>
      static constexpr T_&& get(T_&& t) { return std::forward<T_>(t); }
    };
    
    template <class F>
    struct ToTerm<F, typename enable_if< concepts::CoordsFunctor<F> >::type >
    {
      using type = FunctorTerm<F, CoordsTerm>;
      
      template <class F_>
      static type get(F_&& f) { return {std::forward<F_>(f), CoordsTerm()}; }
    };
    
  } // end namespace detail
  
  template <class T>
  using ToTerm = detail::ToTerm<T>;
  
  template <class T>
  using ToTerm_t = typename ToTerm<T>::type;
  
  template <class T>
  inline auto toTerm(T&& t) RETURNS( ToTerm<T>::get(std::forward<T>(t)) )

} // end namespace AMDiS
