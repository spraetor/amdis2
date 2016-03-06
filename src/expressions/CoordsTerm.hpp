#pragma once

// std c++ headers
#include <memory>
#include <string>

// AMDiS headers
#include "BasisFunction.hpp"
#include "ElInfo.hpp"
#include "MatrixVector.hpp"
#include "Quadrature.hpp"
#include "SubAssembler.hpp"
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

    CoordsTerm() : X(NULL), C_{Comp} {}
    CoordsTerm(int C_) : X(NULL), C_{C_}
    {
      STATIC_ASSERT( Comp == -2 );
    }

    CoordsTerm(CoordsTerm const&) = default;

    template <class List>
    void insertFeSpaces(List& /*feSpaces*/) const {}

    constexpr static int getDegree()
    {
      return 1;
    }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler,
                     Quadrature* quad,
                     BasisFunction const* basisFct = NULL)
    {
      if (!X)
        X = new DataType;

      if (subAssembler)
        subAssembler->getCoordsAtQPs(elInfo, quad, *X);
      else if (quad)
      {
        const int nPoints = quad->getNumPoints();

        X->change_dim(nPoints);
        for (int i = 0; i < nPoints; i++)
          elInfo->coordToWorld(quad->getLambda(i), (*X)[i]);
      }
      else if (basisFct)
      {
        const int nBasisFct = basisFct->getNumber();

        X->change_dim(nBasisFct);
        for (int i = 0; i < nBasisFct; i++)
          elInfo->coordToWorld(*basisFct->getCoords(i), (*X)[i]);
      }
    }

    value_type evalAtIdx(int iq) const
    {
      return evalImpl(iq, int_<Comp>{});
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
    // return component of WorldVector
    template <int C>
    value_type const& evalImpl(int iq, int_<C>) const
    {
      return (*X)[iq][C_];
    }

    template <int C>
    value_type const& evalImpl(WorldVector<double> const& x, int_<C>) const
    {
      return x[C_];
    }

    // return WorldVector
    value_type const& evalImpl(int iq, int_<-1>) const
    {
      return (*X)[iq];
    }

    value_type const& evalImpl(WorldVector<double> const& x, int_<-1>) const
    {
      return x;
    }

  private:
    using DataType = DenseVector<coords_type>;
    DataType* X;

    int C_;
  };

} // end namespace AMDiS
