#pragma once

// std c++ headers
#include <tuple>
#include <utility>

// AMDiS headers
#include "expressions/BaseTerms.hpp"
#include "expressions/LazyOperatorTerm.hpp"
#include "expressions/TermConcepts.hpp"
// #include "expressions/ComponentView.hpp"
#include "traits/basic.hpp"
#include "traits/traits_fwd.hpp"
#include "traits/traits.hpp"
#include "traits/meta_basic.hpp"
#include "utility/int_seq.hpp"

namespace AMDiS
{
  namespace traits
  {
    /// get the degree of a functor by combining the degrees of the arguments
    template <class F, int N, class = void>
    struct functor_degree
    {
      // if no getDegree method is available, use default value 0
      template <class... Int>
      static constexpr int eval(F, Int...)
      {
        return 0;
      }
    };

    template <class F, int N>
    struct functor_degree<F, N, requires::TermFunctor<F,N>>
    {
      template <class... Int>
      static constexpr int eval(F const& f, Int... d)
      {
        return f.getDegree(d...);
      }
    };

  } // end namespace traits

  // forward declaration
  template <class F, class... T>
  struct FunctorValue
  {
    using type = typename std::result_of<F(typename T::value_type...)>::type;
  };

  template <class F, class... Terms>
  using FunctorValue_t = typename FunctorValue<F, Terms...>::type;

  // the expressions
  // ___________________________________________________________________________

  /// Functor that takes arbitrary number of arguments
  template <class F, class... Terms>
  class FunctorTerm
    : public ShapedTerm_t<FunctorValue_t<F, Terms...>, FunctorTerm<F, Terms...>>,
      public LazyOperatorTerms<Terms...>
  {
    static constexpr int N = sizeof...(Terms);
    static_assert( N > 0, "Nullary Functor not allowed here!" );

    using Self  = FunctorTerm;
    using Super = LazyOperatorTerms<Terms...>;

  public:
    using value_type = FunctorValue_t<F, Terms...>;

    template <class... Terms_,
      class = Requires_t<traits::IsCompatible<Types<Terms...>, Types<Terms_...>>> >
    FunctorTerm(Terms_&&... terms_)
      : Super(std::forward<Terms_>(terms_)...),
        fct{}
    {
      MSG("FunctorTerm()");
    }

    template <class F_, class... Terms_,
      class = Requires_t<traits::IsCompatible<Types<F,Terms...>, Types<F_,Terms_...>>> >
    FunctorTerm(F_&& f_, Terms&&... terms_)
      : Super(std::forward<Terms_>(terms_)...),
        fct(std::forward<F_>(f_))
    {
      MSG("FunctorTerm()");
    }

    ~FunctorTerm() { MSG("~FunctorTerm()"); }

    /// return the required quadrature degree to integrate term
    int getDegree() const
    {
      return getDegreeImpl(IndexSeq{});
    }

    /// eval at point with coordinate x
    value_type operator()(WorldVector<double> const& x) const
    {
      MSG("FunctorTerm::operator()()");
      return evalImpl(x, IndexSeq{});
    }

  private:
    using IndexSeq  = MakeSeq_t<N>; // AMDiS::Seq<0,1,2,...,N-1>
    using FctDegree = traits::functor_degree<F,N>;

    // fct.getDegree(t1.getDegre(), t2.getDegree(), ...)
    template <int... I>
    int getDegreeImpl(AMDiS::Seq<I...>) const
    {
      return FctDegree::eval(fct, Super::getDegree(int_<I>{})...);
    }

    // fct(t1(x), t2(x), t3(x),...)
    template <int... I>
    value_type evalImpl(WorldVector<double> const& x, AMDiS::Seq<I...>) const
    {
      MSG("FunctorTerm::evalImpl(" << this << ")");
      return fct(Super::evalTerm(x, int_<I>{})...);
    }

  private:
    F fct; ///< the functor
  };


  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    template <class F, class... Terms>
    struct category<FunctorTerm<F, Terms...>>
    {
      using value_type = FunctorValue_t<F, Terms...>;
      using size_type  = int;
      using tag        = typename category<value_type>::tag;
    };
    /// \endcond

  } // end namespace traits

} // end namespace AMDiS
