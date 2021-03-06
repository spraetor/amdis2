/** \file basic.hpp */

#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// a helper macro to reduce typing
#define RETURNS(...) \
  -> decltype(__VA_ARGS__) { return (__VA_ARGS__); }

#define RETURNS_REF(...) \
  -> decltype(__VA_ARGS__) & { return (__VA_ARGS__); }

#define RETURNS_RREF(...) \
  -> decltype(__VA_ARGS__) && { return std::move(__VA_ARGS__); }

#define RETURNS_CONST_REF(...) \
  -> decltype(__VA_ARGS__) const & { return (__VA_ARGS__); }

namespace AMDiS
{
    namespace detail
    {
        // workaround for MSVC (problems with alias templates in pack expansion)
        template <class, class T>
        struct InvokeType { using type = T; };

        template <class, class, class T>
        struct InvokeType2 { using type = T; };
    }

  template <class T>
  using Value_t = typename detail::InvokeType<T, typename T::value_type>::type;

  template <class T>
  using Size_t = typename detail::InvokeType<T, typename T::size_type>::type;

  template <class T>
  using Result_t = typename detail::InvokeType<T, typename T::result_type>::type;

  template <class T>
  using Decay_t = typename detail::InvokeType<T, typename std::decay<T>::type>::type;

  template <class T1, class T2>
  using Common_t = typename detail::InvokeType2<T1, T2, typename std::common_type<T1,T2>::type>::type;

  namespace detail
  {
    template <class T, class = void>
    struct assign_type
    {
      using type = T;
    };

  } // end namespace detail

  template <class T>
  using Assign_t = typename detail::assign_type<T>::type;


  // ---------------------------------------------------------------------------


  template <class... Ts>
  struct Types {};

  template <class... Ts>
  using Types_t = Types<Decay_t<Ts>...>;

  template <int... Is>
  struct Ints {};


  // ---------------------------------------------------------------------------


  template <class C, class T = void>
  using enable_if = std::enable_if<C::value, T>;

  template <bool C, class T = void>
  using enable_if_c = std::enable_if<C, T>;

  template <class C, class T = void>
  using disable_if = std::enable_if<!C::value, T>;

  template <bool C, class T = void>
  using disable_if_c = std::enable_if<!C, T>;


  // alias for enable_if
  template <class C, class T = void>
  using Requires_t = typename enable_if<C,T>::type;

  template <bool C, class T = void>
  using Requires_c = typename enable_if_c<C,T>::type;


  // ---------------------------------------------------------------------------


  template <class T>
  struct Id
  {
    using type = T;
  };

  template <class T>
  using Id_t = typename Id<T>::type;


  // ---------------------------------------------------------------------------


  // dummy type
  class no_valid_type {};

  namespace detail
  {
    template <bool Valid, class Type>
    struct Enabler
    {
      using type = typename Type::type;
    };

    template <class Type>
    struct Enabler<false, Type>
    {
      using type = no_valid_type;
    };

    template <bool Valid, class Type>
    using Enabler_t = typename Enabler<Valid, Type>::type;

  } // end namespace detail



  template <class T>
  void print_type(T const& t) {}

} // end namespace AMDiS
