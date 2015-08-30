/** \file concepts_base.hpp */

#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// AMDiS headers
#include <traits/basic.hpp>
#include <traits/meta_basic.hpp>

// macro to generate concept-checks
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
  namespace traits
  {
    template <class T>
    struct HasValueType
    {
    private:
      template <class T1> static typename T1::value_type test(int);
      template <class> static void test(...);

    public:
      constexpr static bool value = !std::is_void<decltype(test<T>(0))>::value;
    };

    template <class T>
    struct HasSizeType
    {
    private:
      template <class T1> static typename T1::size_type test(int);
      template <class> static void test(...);

    public:
      constexpr static bool value = !std::is_void<decltype(test<T>(0))>::value;
    };

    template <class T>
    struct HasResultType
    {
    private:
      template <class T1> static typename T1::result_type test(int);
      template <class> static void test(...);

    public:
      constexpr static bool value = !std::is_void<decltype(test<T>(0))>::value;
    };
  }

  /// Namespace that implements basic type-traits
  namespace concepts
  {
    /// @defgroup concepts Concepts definitions, using constexpr and other c++11 features.

    // source: http://stackoverflow.com/questions/25603240/checking-callable-template-parameter-types
    template <class, class, class = void>
    struct check_functor : false_ {};

    template <class Func, class Return, class... Args>
    struct check_functor<Func, Return(Args...),
           typename enable_if<
           std::is_same<
           decltype(std::declval<Func>()(std::declval<Args>()...)),
           Return
           >, void>::type>
       : true_ {};

    template <class, class, class = void>
    struct check_functor_weak : false_ {};

    template <class Func, class Return, class... Args>
    struct check_functor_weak<Func, Return(Args...),
           typename enable_if<
           std::is_convertible<
           decltype(std::declval<Func>()(std::declval<Args>()...)),
           Return
           >, void>::type>
       : true_ {};

    // -------------------------------------------------------------------------
#if 0
    template <class, class, class = void>
    struct check_vector : false_ {};

    template <class Func, class Return, class... Args>
    struct check_vector<Func, Return(Args...),
           typename enable_if<
           std::is_same<
           decltype(std::declval<Func>()[std::declval<Args>()...]),
           Return
           >, void>::type>
       : true_ {};

    template <class, class, class = void>
    struct check_vector_weak : false_ {};

    template <class Func, class Return, class... Args>
    struct check_vector_weak<Func, Return(Args...),
           typename enable_if<
           std::is_convertible<
           decltype(std::declval<Func>()[std::declval<Args>()...]),
           Return
           >, void>::type>
       : true_ {};
#endif
  } // end namespace concepts


  /// Namespace that implements basic concepts-checks using enable_if and
  /// type-traits from namespace \ref concepts
  namespace requires {}

} // end namespace AMDiS
