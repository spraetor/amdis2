/** \file basic.hpp */

#pragma once

#include <boost/version.hpp>

#include "boost/numeric/ublas/detail/returntype_deduction.hpp"
#if BOOST_VERSION >= 105600
#include <boost/core/enable_if.hpp>
#else
#include <boost/utility/enable_if.hpp>
#endif

#include <utility>
#include <type_traits>

#include "meta_basic.hpp"


// a helper mycro to reduce typing
#define RETURNS(...) \
  noexcept(noexcept(__VA_ARGS__)) \
    -> decltype(__VA_ARGS__) { return (__VA_ARGS__); }


namespace AMDiS 
{
  
  template <class T>
  using Value_t = typename T::value_type;
  
  template <class T>
  using Size_t = typename T::size_type;
  
  template <class T>
  using Result_t = typename T::result_type;

  // introduce some shortcuts for boost::mpl
  // ---------------------------------------
  using boost::enable_if;
  using boost::enable_if_c;
  using boost::disable_if;
  using boost::disable_if_c;
  
#if 0
  template <class... Iters>
  struct MultiIterator
  {
    template <class... Iters_>
    MultiIterator(Iters_&&... iters_) 
      : iters(std::forward<Iters_>(iters_)...) { }
    
    auto begin() {
      using Indices = std::make_index_sequence<sizeof...(Iters)>;
      return begin_impl(Indices());
    }
    
    auto end() {
      using Indices = std::make_index_sequence<sizeof...(Iters)>;
      return end_impl(Indices());
    }
    
  private:    
    template <size_t ... I>
    auto begin_impl(std::index_sequence<I...>) {
      return std::make_tuple(std::get<I>(iters).begin()...);
    }
    template <size_t ... I>
    auto end_impl(std::index_sequence<I...>) {
      return std::make_tuple(std::get<I>(iters).end()...);
    }
    
    using TupleType = std::tuple<Iters...>;
    TupleType iters;
  };
  
  
  template <class... Iters>
  auto make_iter(Iters&&... iters)
  {
    return MultiIterator<Iters...>(std::forward<Iters>(iters)...);
  }
#endif
  
  
  // some traits to test for binary operators on types
  namespace traits 
  {
    namespace detail
    {
      template <class...>
      struct voider { using type = void; };
      
      template <class...Ts>
      using void_t = typename voider<Ts...>::type;
  
      template <class T, class U, class BinOp, class Enable = void>
      struct IsBinopAbleImpl : false_ {};
      
      template <class T, class U, class BinOp>
      struct IsBinopAbleImpl<T, U, BinOp, void_t<
        decltype(std::declval<BinOp>()(std::declval<T>(), 
                                       std::declval<U>() )) > > : true_ {};

    } // end namespace detail
    
    template <class T, class U, class BinOp>
    using IsBinopAble = typename detail::IsBinopAbleImpl<T, U, BinOp>::type;
  
    namespace detail
    {
      struct PlusOp {
        template <class T, class U>
        auto operator()(T t, U u) RETURNS(t + u)
      };
      
      struct MultipliesOp {
        template <class T, class U>
        auto operator()(T t, U u) RETURNS(t * u)
      };

    } // end namespace detail
    
    template <class T, class U = T>
    using is_addable = IsBinopAble<T, U, detail::PlusOp>;
    
    template <class T, class U = T>
    using is_multiplicable = IsBinopAble<T, U, detail::MultipliesOp>;
      
    template <class T>
    using is_trivially_copyable = std::is_trivially_copyable<T>;

    // larger types
    template <class... Ts>
    struct larger_type;
    
    template <class T1, class T2, class... Ts>
    struct larger_type<T1, T2, Ts...>
    {
      using type = typename if_then_else< (sizeof(T1) > sizeof(T2)), 
					  larger_type<T1,Ts...>, 
					  larger_type<T2,Ts...> >::type;
    };
    
    template <class T1, class T2>
    struct larger_type<T1, T2>
    {
      using type = if_then_else<(sizeof(T1) > sizeof(T2)), T1, T2>;
    };
    
    template <class T> struct larger_type<T> { using type = T; };
    
    // maximal size type
    template <class... Es>
    using max_size_type = typename larger_type<Size_t<Es>...>::type;
  } // end namespace traits
  
} // end namespace AMDiS
