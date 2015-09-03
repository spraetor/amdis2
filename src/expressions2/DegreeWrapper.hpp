/** \file DegreeWrapper.hpp */

#pragma once

#include <operations/functors.hpp>
#include <traits/basic.hpp>

namespace AMDiS
{
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
    auto operator()(Args&& ... args) const RETURNS
    (
      fct(std::forward<Args>(args)...)
         )

  protected:
    F fct;
  };


  template <class F, class DegF>
struct DegreeWrapper2 : public FunctorBase
  {
    template <class F_, class DegF_,
              class = Requires_t<traits::IsCompatible<Types<F,DegF>, Types<F_,DegF_>>>>
                                 DegreeWrapper2(F&& fct_, DegF&& degfct_)
                                 : fct{fct_},
                                 degfct{degfct_}
  {}

  template <class... Int>
  constexpr int getDegree(Int... degrees) const
  {
    return degfct(degrees...);
  }

  template <class... Args>
  auto operator()(Args&& ... args) const RETURNS
  (
    fct(std::forward<Args>(args)...)
       )

  protected:
  F fct;
  DegF degfct;
  };


  template <int Degree, class F>
  inline DegreeWrapper<Degree, F> deg(F&& fct)
  {
    return {std::forward<F>(fct)};
  }

  template <class F, class DegF>
  inline DegreeWrapper2<F, DegF> deg(F&& fct, DegF&& degfct)
  {
    return {std::forward<F>(fct), std::forward<DegF>(degfct)};
  }
}
