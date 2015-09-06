/** \file int_seq.hpp */

#pragma once

// from http://stackoverflow.com/questions/17424477/implementation-c14-make-integer-sequence/17426611#17426611
namespace AMDiS
{
  namespace detail
  {
    template <unsigned...> struct Seq { using type = Seq; };

    template <class S1, class S2> struct Concat;

    template <unsigned... I1, unsigned... I2>
    struct Concat<Seq<I1...>, Seq<I2...>>
      : Seq<I1..., (sizeof...(I1)+I2)...> {};

    template <class S1, class S2>
    using Concat_t = typename Concat<S1, S2>::type;
    
  } // end namespace detail

  template <unsigned N> struct IntSeq;
  template <unsigned N> using IntSeq_t = typename IntSeq<N>::type;

  template <unsigned N>
  struct IntSeq 
    : Concat_t<IntSeq_t<N/2>, IntSeq_t<N - N/2>> {};

  template <> struct IntSeq<0> : Seq<> {};
  template <> struct IntSeq<1> : Seq<0> {};
  
} // end namespace AMDiS
