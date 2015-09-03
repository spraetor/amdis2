/** \file coords_expr.hpp */

#pragma once

// std c++ headers
#include <string>
#include <memory>

// AMDiS headers
#include <BasisFunction.h>
#include <ElInfo.h>
#include <MatrixVector.h>
#include <Quadrature.h>
#include <SubAssembler.h>
#include <expressions/LazyOperatorTerm.h>

#include "BaseTerms.hpp"

namespace AMDiS
{
  /// Expression that representy the coordinate vector
  // Comp == -1: WorldVector, Comp == -2: component given in constructor
  // Comp != {-1, -2}: component given as template parameter
  template <int Comp = -1>
  struct CoordsTerm : public VectorTerm<CoordsTerm<Comp>>,
        public LazyOperatorTermBase
  {
    using coords_type = WorldVector<double>;
    using value_type  = if_then_else< Comp == -1, coords_type, double >;

    CoordsTerm() : C_{Comp} {}
    CoordsTerm(int C_) : C_{C_}
    {
      STATIC_ASSERT( Comp == -2 );
    }

    template <class List>
    void insertFeSpaces(List& feSpaces) const {}

    constexpr static int getDegree()
    {
      return 1;
    }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler,
                     Quadrature* quad,
                     BasisFunction const* basisFct = NULL)
    {
      if (subAssembler)
        subAssembler->getCoordsAtQPs(elInfo, quad, X);
      else if (quad)
      {
        const int nPoints = quad->getNumPoints();

        X.change_dim(nPoints);
        for (int i = 0; i < nPoints; i++)
          elInfo->coordToWorld(quad->getLambda(i), X[i]);
      }
      else if (basisFct)
      {
        const int nBasisFct = basisFct->getNumber();

        X.change_dim(nBasisFct);
        for (int i = 0; i < nBasisFct; i++)
          elInfo->coordToWorld(*basisFct->getCoords(i), X[i]);
      }
    }

    value_type operator[](int iq) const
    {
      return eval(int_<Comp>(), iq);
    }

    value_type operator()(WorldVector<double> const& x) const
    {
      return eval(int_<Comp>(), x);
    }


    std::string str() const
    {
      return std::string("X") + (Comp == -1 ? "" : std::string("<") + std::to_string(Comp) + ">");
    }

  protected:
    // return component of WorldVector
    template <int C>
    value_type const& eval(int_<C>, int iq) const
    {
      return X[iq][C_];
    }

    template <int C>
    value_type const& eval(int_<C>, WorldVector<double> const& x) const
    {
      return x[C_];
    }

    // return WorldVector
    value_type const& eval(int_<-1>, int iq) const
    {
      return X[iq];
    }

    value_type const& eval(int_<-1>, WorldVector<double> const& x) const
    {
      return x;
    }

  private:
    using DataType = DenseVector<coords_type>;
    DataType X;

    int C_;
  };

} // end namespace AMDiS
