/** \file DOFVectorTerms.hpp */

#pragma once

#include <AMDiS_fwd.h>
#include <SubAssembler.h>
#include <BasisFunction.h>
#include <DOFVectorBase.h>
#include <expressions/LazyOperatorTerm.h>

#include "BaseTerms.hpp"

namespace AMDiS 
{
  struct _unknown {};

  /// Expressions that extracts the values of a DOFVector at QPs
  template <class Vector, class Name, class = void>
  struct ValueOf : public LazyOperatorTermBase {};
  
  template <class T, class Name>
  struct ValueOf<DOFVector<T>, Name> 
    : public ShapedTerm_t<T, ValueOf<DOFVector<T>, Name> >,
      public LazyOperatorTermBase
  {
    using value_type = T;
    using id         = Name;

    constexpr ValueOf(DOFVectorBase<T> const& vector) : vecDV(&vector) {}
    constexpr ValueOf(DOFVectorBase<T> const* vector) : vecDV(vector) {}

    template <class List>
    void insertFeSpaces(List& feSpaces) const
    {
      feSpaces.insert(vecDV->getFeSpace());
    }
    
    int getDegree() const
    {
      return vecDV->getFeSpace()->getBasisFcts()->getDegree();
    }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler, 
                     Quadrature* quad, 
                     BasisFunction const* basisFct = NULL)
    {
      if (subAssembler)
        subAssembler->getVectorAtQPs(vecDV, elInfo, quad, values);
      else if (quad)
        vecDV->getVecAtQPs(elInfo, quad, NULL, values);
      else if (basisFct) {
        const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
        
        // get coefficients of DOFVector
        DenseVector<T> coeff{(size_t)localBasisFct->getNumber()};
        vecDV->getLocalVector(elInfo->getElement(), coeff);
        
        // eval basisfunctions of DOFVector at coords of given basisFct
        const int nBasisFct = basisFct->getNumber();
        values.change_dim(nBasisFct);	
        for (int i = 0; i < nBasisFct; i++)
          values[i] = localBasisFct->evalUh(*basisFct->getCoords(i), coeff);
      }
    }
    
    value_type operator()(int iq) const { return values[iq]; }
    
    std::string str() const { return std::string("value(") + vecDV->getName() + ")"; }
    
  private:
    DOFVectorBase<T> const* vecDV;
    
    DenseVector<T> values;
  };
  
  
  // ---------------------------------------------------------------------------

  
  /// Expressions that extracts the values of a DOFVector at QPs
  template <class Vector, class Name, class = void>
  struct GradientOf : public LazyOperatorTermBase {};
  
  template <class T, class Name>
  struct GradientOf<DOFVector<T>, Name> 
    : public ShapedTerm_t<Gradient_t<T>, GradientOf<DOFVector<T>, Name> >,
      public LazyOperatorTermBase
  {
    using value_type = Gradient_t<T>;
    using id         = Name;

    constexpr GradientOf(DOFVectorBase<T> const& vector) : vecDV(&vector) {}
    constexpr GradientOf(DOFVectorBase<T> const* vector) : vecDV(vector) {}

    template <class List>
    void insertFeSpaces(List& feSpaces) const
    {
      feSpaces.insert(vecDV->getFeSpace());
    }
    
    int getDegree() const
    {
      return vecDV->getFeSpace()->getBasisFcts()->getDegree();
    }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler, 
                     Quadrature* quad, 
                     BasisFunction const* basisFct = NULL)
    {
      if (subAssembler)
        subAssembler->getGradientsAtQPs(vecDV, elInfo, quad, values);
      else if (quad)
        vecDV->getGrdAtQPs(elInfo, quad, NULL, values);
      else if (basisFct) {
        const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
        
        // get coefficients of DOFVector
        DenseVector<T> coeff{(size_t)localBasisFct->getNumber()};
        vecDV->getLocalVector(elInfo->getElement(), coeff);
        
        // eval basisfunctions of DOFVector at coords of given basisFct
        const int nBasisFct = basisFct->getNumber();
        values.change_dim(nBasisFct); 
        auto& grdLambda = elInfo->getGrdLambda();
        for (int i = 0; i < nBasisFct; i++)
          localBasisFct->evalGrdUh(*basisFct->getCoords(i), grdLambda, coeff, values[i]);
      }
    }
    
    value_type operator()(int iq) const { return values[iq]; }
    
    std::string str() const { return std::string("grad(") + vecDV->getName() + ")"; }
    
  private:
    DOFVectorBase<T> const* vecDV;
    
    DenseVector<value_type> values;
  };
  
} // end namespace AMDiS
