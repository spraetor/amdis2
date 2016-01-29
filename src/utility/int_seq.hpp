/** \file int_seq.hpp */

#pragma once

// from http://stackoverflow.com/questions/17424477/implementation-c14-make-integer-sequence/17426611#17426611
namespace AMDiS
{
  template <int...> struct Seq {};

  namespace detail
  {
    template <int s, class S>
    struct Concat;

    template <int s, int ... i>
    struct Concat<s, Seq<i... >>
    {
      using type = Seq<i..., (s + i)... >;
    };

    template <bool, class S>
    struct IncSeq_if 
    { 
      using type = S; 
    };

    template <int... Is>
    struct IncSeq_if<true, Seq<Is...>>
    {
      using type = Seq<Is..., sizeof...(Is)>;
    };
    
  } // end namespace detail
  
  template <int N>
  struct MakeSeq 
  {
    using type = typename detail::IncSeq_if< (N % 2 != 0), 
      typename detail::Concat<N/2, typename MakeSeq<N/2>::type>::type >::type;
  };

  // break condition
  template <> struct MakeSeq<0> { using type = Seq<>; };

  // alias template
  template <int N>
  using MakeSeq_t = typename MakeSeq<N>::type;
  
} // end namespace AMDiS
