/** \file hessianOf.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
#include "DOFVector.h"
#include "traits/category.hpp"

namespace AMDiS 
{  
  namespace expressions 
  {  
    /// Expressions that extracts the gradient of a DOFVector at QPs
    template <class Vector, class Name>
    struct HessianOf : public LazyOperatorTermBase
    {
      using T = Value_t<traits::category<Vector> >;
      using value_type = typename D2Type<T>::type;
      using id = Name;

      DOFVector<T>*                          vecDV;
      mutable mtl::dense_vector<value_type>  vec;
      mutable mtl::dense_vector<T>           coeff;

      HessianOf(Vector& vector) : vecDV(&vector) {}
      HessianOf(Vector* vector) : vecDV(vector) {}

      template <class List>
      void insertFeSpaces(List& feSpaces) const
      {
        feSpaces.insert(vecDV->getFeSpace());
      }
      
      int getDegree() const
      {
        return vecDV->getFeSpace()->getBasisFcts()->getDegree() /* -1 */;
      }

      void initElement(const ElInfo* elInfo,
            		       SubAssembler* subAssembler, Quadrature *quad, 
            		       const BasisFunction *basisFct = NULL)
      { FUNCNAME("HessianOf::initElement");
      
      	if (subAssembler) {
      	  ERROR_EXIT("Hessian expression not yet implemented for operator-terms!\n");
      	} 
        else if (quad)
      	  vecDV->getD2AtQPs(elInfo, quad, NULL, vec);
      	else if (basisFct) {
      	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
      	  
      	  // get coefficients of DOFVector
      	  coeff.change_dim(localBasisFct->getNumber());
      	  vecDV->getLocalVector(elInfo->getElement(), coeff);
      	  
      	  // eval basisfunctions of DOFVector at coords of given basisFct
      	  size_t nBasisFct = basisFct->getNumber();
      	  vec.change_dim(nBasisFct);	
      	  
      	  const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
      	  for (size_t i = 0; i < nBasisFct; i++)
      	    localBasisFct->evalD2Uh(*basisFct->getCoords(i), grdLambda, coeff, &vec[i]);
      	}
      }

      value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("hessian(") + vecDV->getName() + ")"; }
    };
      
    
    /// Expressions that extracts the partial derivative of a DOFVector at QPs
    template <class Vector, class Name>
    struct LaplacianOf : public LazyOperatorTermBase
    {
      using T = Value_t<traits::category<Vector> >;
      using value_type = T;
      using id = Name;

      DOFVector<T>*                                        vecDV;
      mutable mtl::dense_vector<typename D2Type<T>::type>  vec_tmp;
      mutable mtl::dense_vector<value_type>                vec;
      mutable mtl::dense_vector<T>                         coeff;
      int comp;

      LaplacianOf(Vector& vector) : vecDV(&vector) { }
      LaplacianOf(Vector* vector) : vecDV(vector) { }

      template<typename List>
      void insertFeSpaces(List& feSpaces) const
      {
        feSpaces.insert(vecDV->getFeSpace());
      }
      
      int getDegree() const
      {
        return vecDV->getFeSpace()->getBasisFcts()->getDegree() /* -1 */;
      }

      void initElement(const ElInfo* elInfo,
            		       SubAssembler* subAssembler, Quadrature *quad, 
            		       const BasisFunction *basisFct = NULL)
      { FUNCNAME("LaplacianOf::initElement");
            
      	if (subAssembler) {
      	  ERROR_EXIT("Laplacian expression not yet implemented for operator-terms!\n");
      	} else if (quad) {
      	  vecDV->getD2AtQPs(elInfo, quad, NULL, vec_tmp); 
      	  for (size_t i = 0; i < size(vec_tmp); i++) {
      	    vec[i] = vec_tmp[i][0][0];
      	    for (int j = 1; j < Global::getGeo(WORLD); j++)
      	      vec[i] += vec_tmp[i][j][j];
      	  }
      	} else if (basisFct) {
      	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
      	  
      	  // get coefficients of DOFVector
      	  coeff.change_dim(localBasisFct->getNumber());
      	  vecDV->getLocalVector(elInfo->getElement(), coeff);
      	  
      	  // eval basisfunctions of DOFVector at coords of given basisFct
      	  size_t nBasisFct = basisFct->getNumber();
      	  vec.change_dim(nBasisFct);	
      	  WorldMatrix<double> hessian;
      	  
      	  const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
      	  for (size_t i = 0; i < nBasisFct; i++) {
      	    localBasisFct->evalD2Uh(*basisFct->getCoords(i), grdLambda, coeff, &hessian);
      	    vec[i] = hessian[0][0];
      	    for (int j = 1; j < Global::getGeo(WORLD); j++)
      	      vec[i] += hessian[j][j];
      	  }
      	}
      }

      value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("laplace(") + vecDV->getName() + ")"; }
    };
    
  } // end namespace expressions

  
  // gradient of a DOFVector
  // _____________________________________________________________________________

  template <class Name = _unknown, class T>
  expressions::HessianOf<DOFVector<T>, Name > 
  hessianOf(DOFVector<T>& vector) { return {vector}; }

  template <class Name = _unknown, class T>
  expressions::HessianOf<DOFVector<T>, Name > 
  hessianOf(DOFVector<T>* vector) { return {vector}; }


  // Partial derivative of a DOFVector
  // _____________________________________________________________________________

  // with Name
  template <class Name = _unknown, class T>
  expressions::LaplacianOf<DOFVector<T>, Name > 
  laplacianOf(DOFVector<T>& vector) { return {vector}; }

  template <class Name = _unknown, class T>
  expressions::LaplacianOf<DOFVector<T>, Name > 
  laplacianOf(DOFVector<T>* vector) { return {vector}; }

} // end namespace AMDiS
