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

#include "BaseTerms.hpp"
#include "LazyOperatorTerm.h"

namespace AMDiS 
{    
  /// Expression that representy the coordinate vector
  struct CoordsTerm
    : public VectorTerm<CoordsTerm>,
      public LazyOperatorTermBase
  {
    using Self       = CoordsTerm;
    using Super      = LazyOperatorTermBase;    
    using value_type = WorldVector<double>;

    constexpr CoordsTerm() {}

    template <class List>
    void insertFeSpaces(List& feSpaces) const {}
    
    static constexpr int getDegree() { return 1; }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler, 
                     Quadrature* quad, 
                     BasisFunction const* basisFct = NULL)
    {
      x = DataType(new DenseVector<value_type>());
      if (subAssembler)
        subAssembler->getCoordsAtQPs(elInfo, quad, *x); 
      else if (quad) {
        const int nPoints = quad->getNumPoints();
      
        x->change_dim(nPoints);
        for (int i = 0; i < nPoints; i++)
          elInfo->coordToWorld(quad->getLambda(i), (*x)[i]);
      }
      else if (basisFct) {
        const int nBasisFct = basisFct->getNumber();
      
        x->change_dim(nBasisFct);
        for (int i = 0; i < nBasisFct; i++)
          elInfo->coordToWorld(*basisFct->getCoords(i), (*x)[i]);
      }
    }

    value_type operator()(int iq) const { return (*x)[iq]; }
    
    std::string str() const { return "X"; }
    
  private:
    using DataType = std::shared_ptr<DenseVector<value_type>>;
    DataType x;
  };
    
} // end namespace AMDiS
