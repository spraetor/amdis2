#pragma once

#include "expressions/LazyOperatorTermBase.hpp"
#include "expressions/TermConcepts.hpp"
#include "traits/basic.hpp"
#include "traits/traits.hpp"
#include "utility/foreach.hpp"

namespace AMDiS
{
  /// Operator term with arbitrary number of sub-term (expressions)
  template <class... Terms>
  class LazyOperatorTerms
    : public LazyOperatorTermBase
  {
    using Self      = LazyOperatorTerms;
    using TermTuple = std::tuple<Terms...>;

    template <int N>
    using Term_t = typename std::tuple_element<size_t(N), TermTuple>::type;

  public:
    template <class... Terms_,
        class = Requires_t< concepts::Term<Terms_...> > >
    LazyOperatorTerms(Terms_&&... terms_)
      : terms(std::forward<Terms_>(terms_)...)
    {
      MSG("LazyOperatorTerms(" << this << ")");
    }

    ~LazyOperatorTerms() { MSG("~LazyOperatorTerms(" << this << ")"); }

    template <int N>
    int getDegree(int_<N>) const
    {
      return std::get<N>(terms).getDegree();
    }

    template <int N>
    Term_t<N>& getTerm(int_<N>)
    {
      return std::get<N>(terms);
    }

    template <int N>
    Term_t<N> const& getTerm(int_<N>) const
    {
      return std::get<N>(terms);
    }

  protected:
    template <int N>
    Value_t<Term_t<N>> evalTerm(WorldVector<double> const& x, int_<N>) const
    {
      MSG("LazyOperatorTerms::evalTerm(" << &std::get<N>(terms) << " => " << typeid(std::get<N>(terms)).name() << ")");
      return std::get<N>(terms)(x);
    }

  private:
    TermTuple terms;
  };

} // end namespace AMDiS
