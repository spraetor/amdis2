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

    constexpr CoordsTerm() : x(NULL), C_{Comp} {}    
    constexpr CoordsTerm(int C_) : x(NULL), C_{C_} { STATIC_ASSERT( Comp == -2 ); }

    template <class List>
    void insertFeSpaces(List& feSpaces) const {}
    
    constexpr static int getDegree() { return 1; }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler, 
                     Quadrature* quad, 
                     BasisFunction const* basisFct = NULL)
    {
      if (!x)
        x = new DenseVector<coords_type>();
      
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
      
      for (size_t i = 0; i < size(*x); ++i)
        std::cout << "x[iq=" << std::to_string(i) << "] = " << (*x)[i] << "\n";
      
      value_type result = (*x)[0];
    }

    value_type operator()(int iq) const { return eval(int_<Comp>(), iq); }
    
    std::string str() const { return std::string("X") + (Comp == -1 ? "" : std::string("<") + std::to_string(Comp) + ">"); }
    
  protected:
    // return component of WorldVector
    template <int C>
    value_type eval(int_<C>, int iq) const { 
      TEST_EXIT_DBG( x )("Coords-vector not initialized!\n");
      return (*x)[iq][C_]; 
    }
    
    // return WorldVector
    value_type eval(int_<-1>, int iq) const { 
      TEST_EXIT_DBG( x && size(*x)>0 && size((*x)[0])>0 )("Coords-vector not initialized!\n");
      for (size_t i = 0; i < size(*x); ++i)
        std::cout << ">>> x[iq=" << std::to_string(i) << "] = " << (*x)[i] << "\n";
      std::cout << "iq = " << iq << "\n";
      return (*x)[iq]; 
    }
    
  private:
    using DataType = DenseVector<coords_type>;
    DataType* x;
    
    int C_;
  };
    
} // end namespace AMDiS
