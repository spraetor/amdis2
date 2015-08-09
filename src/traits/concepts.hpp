/** \file concepts.hpp */

#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// AMDiS headers
#include <traits/meta_basic.hpp>
#include <AMDiS_fwd.h>
#include <MatrixVector_fwd.h>


#define HAS_MEMBER_GENERATE(name) \
  template <class,class,class = void> struct has_member_ ## name : false_ {}; \
  template <class F, class Return, class... Args>                             \
  struct has_member_ ## name <F, Return(Args...),                             \
    typename enable_if<                                                       \
      std::is_same<                                                           \
        decltype(std::declval<F>().name (std::declval<Args>()...)),           \
        Return                                                                \
      >, void>::type>                                                         \
  : true_ {};
  
#define HAS_MEMBER(name) has_member_ ## name

namespace AMDiS
{
  // forward declaration
  class SubAssembler;
  
  namespace concepts
  {
    // source: http://stackoverflow.com/questions/25603240/checking-callable-template-parameter-types
    template <class, class, class = void>
    struct check_signature : false_ {};

    template <class Func, class Return, class... Args>
    struct check_signature<Func, Return(Args...),
      typename enable_if<
          std::is_same<
            decltype(std::declval<Func>()(std::declval<Args>()...)),
            Return
          >, void>::type>
      : true_ {};
	
    template <class, class, class = void>
    struct check_signature_weak : false_ {};

    template <class Func, class Return, class... Args>
    struct check_signature_weak<Func, Return(Args...),
      typename enable_if<
          std::is_convertible<
            decltype(std::declval<Func>()(std::declval<Args>()...)),
            Return
          >, void>::type>
      : true_ {};
	
    
    HAS_MEMBER_GENERATE(insertFeSpaces)
    HAS_MEMBER_GENERATE(initElement)
    HAS_MEMBER_GENERATE(getDegree)
    
    /* ---------------------------------------------------------------------- */
    
    // expressions concepts
    template <class... Ts>
    struct Term
    {
      static constexpr bool value = and_<Term<Ts>...>::value;
    };
	
    template <class T>
    struct Term<T>
    {
    private:
      using Type = typename std::decay<T>::type;
      
      template <class T1>
      static decltype(/*Size_t<T1>(), */Value_t<T1>(), void(), 0) test(int);
      template <class>
      static void test(...);
      
      static constexpr bool value0 = HAS_MEMBER(insertFeSpaces)<Type, void(int&)>::value;
      static constexpr bool value1 = HAS_MEMBER(initElement)<Type, void(ElInfo const*, SubAssembler*, Quadrature*, BasisFunction const*)>::value;
      static constexpr bool value2 = !std::is_void<decltype(test<Type>(0))>::value;
      
      template <bool, class T1>
      struct check {
        static constexpr bool value = check_signature_weak<T1, Value_t<T1>(int)>::value;
      };
      template <class T1> struct check<false, T1> : false_ {};
      
    public:
      static constexpr bool value = value0 && value1 && check<value2, Type>::value;
    };
    
    // check wether a functor has a member getDegree, with N arguments of type int
    template <class F, int N>
    struct TermFunctor
    {
      template <int I, class... Ints>
      static constexpr bool has_member(int_<I>, Ints...)
      {
        return has_member(int_<I-1>(), 0, Ints(0)...);
      }
      
      template <class... Ints>
      static constexpr bool has_member(int_<0>, Ints...)
      {
        return HAS_MEMBER(getDegree)<F, int(Ints...)>::value;
      }
      
      static constexpr bool value = has_member(int_<N>());
    };
    
    template <class F>
    using CoordsFunctor = check_signature<F, double(WorldVector<double>)>;
    
  } // end namespace concepts
    
  namespace requires
  {
    // test whether one of the arguments is term
    template <class Result, class... Args>
    using Term = typename enable_if< or_<concepts::Term<Args>...>, Result >::type;
    
    template <class F, int N>
    using TermFunctor = typename enable_if< concepts::TermFunctor<F,N> >::type;
    
    template <class Result, class F>
    using CoordsFunctor = typename enable_if< concepts::CoordsFunctor<F>, Result >::type;
    
  } // end namespace require
} // end namespace AMDiS
