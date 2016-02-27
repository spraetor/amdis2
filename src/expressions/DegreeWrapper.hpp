#pragma once

#include "operations/functors.hpp"
#include "traits/basic.hpp"

namespace AMDiS
{
  /// class is instatiated be generator function \ref deg
  template <int D, class F>
  struct DegreeWrapper : public FunctorBase
  {
    template <class F_,
      class = Requires_t<traits::IsCompatible<F, F_>>>
    DegreeWrapper(F_&& fct_)
      : fct{fct_}
    {}

    template <class... Int>
    constexpr int getDegree(Int... degrees) const
    {
      return D;
    }

    template <class... Args>
    auto operator()(Args&&... args) const RETURNS
    (
      fct(std::forward<Args>(args)...)
    )

  protected:
    F fct;
  };


  /// class is instatiated be generator function \ref deg
  template <class F, class DegF>
  struct DegreeWrapper2 : public FunctorBase
  {
    template <class F_, class DegF_,
      class = Requires_t<traits::IsCompatible<Types<F,DegF>, Types<F_,DegF_>>>>
    DegreeWrapper2(F_&& fct_, DegF_&& degfct_)
      : fct{fct_},
	    degfct{degfct_}
    {}

    template <class... Int>
    constexpr int getDegree(Int... degrees) const
    {
      return degfct(degrees...);
    }

    template <class... Args>
    auto operator()(Args&&... args) const RETURNS
    (
      fct(std::forward<Args>(args)...)
    )

  protected:
    F fct;
    DegF degfct;
  };


  // add a constant quadrature degree to the functor
  template <int Degree, class F>
  inline DegreeWrapper<Degree, F> deg(F const& fct)
  {
    return {fct};
  }

  // add a quadrature degree to a functor by providing a
  // degree-functor, that takes the degree of the arguments passed
  // to the functor and returns the combined degree.
  template <class F, class DegF>
  inline DegreeWrapper2<F, DegF> deg(F const& fct, DegF const& degfct)
  {
    return {fct, degfct};
  }
  
} // end namespace AMDiS
