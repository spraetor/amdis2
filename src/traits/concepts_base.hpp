/** \file concepts_base.hpp */

#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// AMDiS headers
#include "traits/basic.hpp"
#include "traits/meta_basic.hpp"

// macro to generate concept-checks
#define HAS_MEMBER_GENERATE(name)                                           \
  template <class F, class... Args>                                         \
  inline auto invoke_member_ ## name(F&& f, Args&&... args) ->              \
    decltype(std::forward<F>(f).name (std::forward<Args>(args)...)) {       \
    return std::forward<F>(f).name (std::forward<Args>(args)...);           \
  }                                                                         \
  inline no_valid_type invoke_member_ ## name(...);                         \
  template <class, class> struct has_member_ ## name : false_ {};           \
  template <class F, class Return, class... Args>                           \
  struct has_member_ ## name <F, Return(Args...)>                           \
    : std::is_convertible<                                                  \
        decltype(invoke_member_ ## name(std::declval<F>(), std::declval<Args>()...)), \
        Return > {}; 

#define HAS_MEMBER(name) has_member_ ## name

namespace AMDiS
{
  namespace traits
  {
    // Some type traits to test for type-attributes of a class
    // -------------------------------------------------------------------------

    namespace detail
    {
      template <class T, class = typename T::value_type >
      true_ HasValueTypeImpl(int);

      template <class T> false_ HasValueTypeImpl(...);


      template <class T, class = typename T::size_type >
      true_ HasSizeTypeImpl(int);

      template <class T> false_ HasSizeTypeImpl(...);


      template <class T, class = typename T::result_type >
      true_ HasResultTypeImpl(int);

      template <class T> false_ HasResultTypeImpl(...);

    } // end namespace detail


    template <class T>
    using HasValueType = decltype(detail::HasValueTypeImpl<T>(int{}));

    template <class T>
    using HasSizeType = decltype(detail::HasSizeTypeImpl<T>(int{}));

    template <class T>
    using HasResultType = decltype(detail::HasResultTypeImpl<T>(int{}));


    // Type traits to test whether a class is a functor, i.e. has a operator()
    // -------------------------------------------------------------------------
    
    namespace detail 
    {
        template <class F, class... Args>                                         
        inline auto invoke_functor(F&& f, Args&&... args) RETURNS
        (      
            std::forward<F>(f)(std::forward<Args>(args)...)
        )

        inline no_valid_type invoke_functor(...);
    }
    
    template <class, class>
    struct IsFunctor : false_ {};

    template <class F, class Return, class... Args>
    struct IsFunctor<F, Return(Args...)>
        : std::is_same< decltype(detail::invoke_functor(std::declval<F>(), std::declval<Args>()...)),
                        Return > {};

    template <class, class>
    struct IsFunctorWeak : false_ {};

    template <class F, class Return, class... Args>
    struct IsFunctorWeak<F, Return(Args...)>
        : std::is_convertible< decltype(detail::invoke_functor(std::declval<F>(), std::declval<Args>()...)),
                               Return > {};

    // -------------------------------------------------------------------------

    template <class>
    struct IsCallable : false_ {};

    template <class F, class... Args>
    struct IsCallable<F(Args...)>
        : not_< std::is_same< decltype(detail::invoke_functor(std::declval<F>(), std::declval<Args>()...)),
                              no_valid_type > > {};
                 
    // -------------------------------------------------------------------------
                 
    namespace detail 
    {
        template <class F, class Arg>                                         
        inline auto invoke_vector(F&& f, Arg&& arg) RETURNS
        (      
            std::forward<F>(f)[std::forward<Arg>(arg)]
        )

        inline no_valid_type invoke_vector(...);
    }
    
    template <class, class>
    struct IsVector : false_ {};

    template <class F, class Return, class Arg>
    struct IsVector<F, Return(Arg)>
        : std::is_same< decltype(detail::invoke_vector(std::declval<F>(), std::declval<Arg>())),
                        Return > {};

    template <class, class>
    struct IsVectorWeak : false_ {};

    template <class F, class Return, class Arg>
    struct IsVectorWeak<F, Return(Arg)>
        : std::is_convertible< decltype(detail::invoke_vector(std::declval<F>(), std::declval<Arg>())),
                               Return > {}; 
                               
                               
    namespace detail
    {
      template <class T, class = decltype(&T::push_back())>
      true_ PushBackImpl(int);

      template <class T> false_ PushBackImpl(...);

    } // end namespace detail
    
    template <class T> 
    using HasPushBack = decltype(detail::PushBackImpl<T>(int{}));
    
  } // end namespace traits


  /// Namespace that implements basic concept checks
  namespace concepts
  {
    /// @defgroup concepts Concepts definitions, using decltype and other c++11 features.

//     namespace detail
//     {
//       template <class T, class = decltype(&T::operator[])>
//       true_ VectorImpl(int);
// 
//       template <class T> false_ VectorImpl(...);
// 
//     } // end namespace detail

    // template <class T>
    // using Vector = decltype(detail::VectorImpl<T>(int{}));

    template <class T>
    using Vector = traits::IsVectorWeak<T, Value_t<T>(size_t)>;

    template <class T>
    using Matrix = traits::IsFunctorWeak<T, Value_t<T>(size_t, size_t)>;

  } // end namespace concepts


  /// Namespace that implements basic concepts-checks using enable_if and
  /// type-traits from namespace \ref concepts
  namespace requires {}

} // end namespace AMDiS
