#pragma once

// std c++ headers
#include <memory>
#include <string>

// AMDiS headers
#include "MatrixVector.hpp"
#include "expressions/BaseTerms.hpp"
#include "expressions/LazyOperatorTerm.hpp"

namespace AMDiS
{
  namespace detail
  {
    template <int I, class M>
    using ShapeByIdx = if_then_else< I == -1, VectorTerm<M>, BaseTerm<M> >;
  }

  /// Expression that representy the coordinate vector
  // Comp == -1: WorldVector, Comp == -2: component given in constructor
  // Comp != {-1, -2}: component given as template parameter
  template <int Comp = -1>
  struct CoordsTerm
    : public detail::ShapeByIdx<Comp, CoordsTerm<Comp>>,
      public LazyOperatorTermBase
  {
    using coords_type = WorldVector<double>;
    using value_type  = if_then_else< Comp == -1, coords_type, double >;

    CoordsTerm() : C_{Comp} {MSG("CoordsTerm()");}
    CoordsTerm(int C_) : C_{C_}
    {
      MSG("CoordsTerm(c)");
      STATIC_ASSERT( Comp == -2 );
    }

//     CoordsTerm(CoordsTerm const&) = default;
    CoordsTerm(CoordsTerm&&) = default;
    CoordsTerm& operator=(CoordsTerm&&) = default;

    ~CoordsTerm() { MSG("~CoordsTerm()"); }

    constexpr static int getDegree()
    {
      return 1;
    }

    value_type operator()(WorldVector<double> const& x) const
    {
      return evalImpl(x, int_<Comp>{});
    }

    std::string str() const
    {
      return std::string("X") + (Comp == -1 ? "" : std::string("<") + std::to_string(Comp) + ">");
    }

  protected:

    template <int C>
    value_type const& evalImpl(WorldVector<double> const& x, int_<C>) const
    {
      return x[C_];
    }

    value_type const& evalImpl(WorldVector<double> const& x, int_<-1>) const
    {
      return x;
    }

  private:
    int C_;
  };

} // end namespace AMDiS
