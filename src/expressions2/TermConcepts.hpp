/** \file TermConcepts.hpp */

#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// AMDiS headers
#include <traits/basic.hpp>
#include <traits/meta_basic.hpp>
#include <traits/concepts_base.hpp>

#include <MatrixVector_fwd.h> // WorldVector<T>

namespace AMDiS
{
  // forward declaration
  class SubAssembler;
  class ElInfo;
  class Quadrature;
  class BasisFunction;
  
  namespace concepts
  {    
    HAS_MEMBER_GENERATE(insertFeSpaces)
    HAS_MEMBER_GENERATE(initElement)
    HAS_MEMBER_GENERATE(getDegree)
    
    // expressions concepts
    template <class... Ts>
    struct Term
    {
      constexpr static bool value = and_<Term<Ts>...>::value;
    };
	
    template <class T>
    struct Term<T>
    {
    private:
      using Type = typename std::decay<T>::type;
      
      template <class T1> static Value_t<T1> test(int);
      template <class> static void test(...);
      
      constexpr static bool value0 = HAS_MEMBER(insertFeSpaces)<Type, void(int&)>::value;
      constexpr static bool value1 = HAS_MEMBER(initElement)<Type, void(ElInfo const*, SubAssembler*, Quadrature*, BasisFunction const*)>::value;
      constexpr static bool value2 = !std::is_void<decltype(test<Type>(0))>::value;
      
      template <bool, class T1>
      struct check {
        constexpr static bool value = check_signature_weak<T1, Value_t<T1>(int)>::value;
      };
      template <class T1> struct check<false, T1> : false_ {};
      
    public:
      constexpr static bool value = value0 && value1 && check<value2, Type>::value;
    };
        
    
    // check wether a functor has a member getDegree, with N arguments of type int
    template <class F, int N>
    struct TermFunctor
    {
      template <int I, class... Ints>
      constexpr static bool has_member(int_<I>, Ints...)
      {
        return has_member(int_<I-1>(), 0, Ints(0)...);
      }
      
      template <class... Ints>
      constexpr static bool has_member(int_<0>, Ints...)
      {
        return HAS_MEMBER(getDegree)<F, int(Ints...)>::value;
      }
      
      constexpr static bool value = has_member(int_<N>());
    };
    
    
    template <class F>
    using CoordsFunctor = check_signature<F, double(WorldVector<double>)>;
    
  } // end namespace concepts
    
    
  namespace requires
  {
    // test whether one of the arguments is term
    template <class Result, class... Args>
    using Term = Requires_t< or_<concepts::Term<Args>...>, Result >;
    
    template <class F, int N>
    using TermFunctor = Requires_t< concepts::TermFunctor<F,N> >;
    
    template <class Result, class F>
    using CoordsFunctor = Requires_t< concepts::CoordsFunctor<F>, Result >;
    
  } // end namespace require
} // end namespace AMDiS
