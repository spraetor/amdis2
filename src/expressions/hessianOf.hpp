/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file hessianOf.hpp */

#ifndef AMDIS_HESSIAN_OF_HPP
#define AMDIS_HESSIAN_OF_HPP

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
#include "DOFVector.h"
#include "traits/category.hpp"

namespace AMDiS 
{  
  namespace expressions 
  {  
    /// Expressions that extracts the gradient of a DOFVector at QPs
    template<typename Vector, typename Name>
    struct HessianOf : public LazyOperatorTermBase
    {
      typedef typename traits::category<Vector>::value_type  T;
      typedef typename D2Type<T>::type                       value_type;
      typedef Name                                           id;

      DOFVector<T>*                          vecDV;
      mutable mtl::dense_vector<value_type>  vec;
      mutable mtl::dense_vector<T>           coeff;

      HessianOf(Vector& vector) : vecDV(&vector) {}
      HessianOf(Vector* vector) : vecDV(vector) {}

      template<typename List>
      void insertFeSpaces(List& feSpaces) const
      {
	feSpaces.insert(vecDV->getFeSpace());
      }
      
      int getDegree() const
      {
	return vecDV->getFeSpace()->getBasisFcts()->getDegree() /* -1 */;
      }

      template<typename OT>
      void initElement(OT* ot, const ElInfo* elInfo,
		      SubAssembler* subAssembler, Quadrature *quad, 
		      const BasisFunction *basisFct = NULL)
      { FUNCNAME("HessianOf::initElement");
      
	if (ot && subAssembler) {
	  ERROR_EXIT("Hessian expression not yet implemented for operator-terms!\n");
	} else if (quad)
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


      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      { FUNCNAME("HessianOf::initElement");
      
	ERROR_EXIT("Hessian expression not yet implemented for Dual-mesh!\n");
      }

      value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("hessian(") + vecDV->getName() + ")"; }
    };
      
    
    /// Expressions that extracts the partial derivative of a DOFVector at QPs
    template<typename Vector, typename Name>
    struct LaplacianOf : public LazyOperatorTermBase
    {
      typedef typename traits::category<Vector>::value_type   T;
      typedef T                                               value_type;
      typedef Name                                            id;

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

      template<typename OT>
      void initElement(OT* ot, const ElInfo* elInfo,
		      SubAssembler* subAssembler, Quadrature *quad, 
		      const BasisFunction *basisFct = NULL)
      { FUNCNAME("LaplacianOf::initElement");
      
	if (ot && subAssembler) {
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

      template<typename OT>
      void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
		      SubAssembler* subAssembler, Quadrature *quad, 
		      const BasisFunction *basisFct = NULL)
      { FUNCNAME("LaplacianOf::initElement");
      
	ERROR_EXIT("Laplacian expression not yet implemented for Dual-mesh!\n");
      }

      value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("laplace(") + vecDV->getName() + ")"; }
    };
    
  } // end namespace expressions

  
  // gradient of a DOFVector
  // _____________________________________________________________________________

  // with Name
  template<typename Name, typename T>
  expressions::HessianOf<DOFVector<T>, Name > hessianOf(DOFVector<T>& vector) 
  { return expressions::HessianOf<DOFVector<T>, Name >(vector); }

  template<typename Name, typename T>
  expressions::HessianOf<DOFVector<T>, Name > hessianOf(DOFVector<T>* vector) 
  { return expressions::HessianOf<DOFVector<T>, Name >(vector); }

  // without Name
  template<typename T>
  expressions::HessianOf<DOFVector<T>, _unknown > hessianOf(DOFVector<T>& vector) 
  { return expressions::HessianOf<DOFVector<T>, _unknown >(vector); }

  template<typename T>
  expressions::HessianOf<DOFVector<T>, _unknown > hessianOf(DOFVector<T>* vector) 
  { return expressions::HessianOf<DOFVector<T>, _unknown >(vector); }


  // Partial derivative of a DOFVector
  // _____________________________________________________________________________

  // with Name
  template<typename Name, typename T>
  expressions::LaplacianOf<DOFVector<T>, Name > laplacianOf(DOFVector<T>& vector) 
  { return expressions::LaplacianOf<DOFVector<T>, Name >(vector); }

  template<typename Name, typename T>
  expressions::LaplacianOf<DOFVector<T>, Name > laplacianOf(DOFVector<T>* vector) 
  { return expressions::LaplacianOf<DOFVector<T>, Name >(vector); }

  // without Name
  template<typename T>
  expressions::LaplacianOf<DOFVector<T>, _unknown > laplacianOf(DOFVector<T>& vector) 
  { return expressions::LaplacianOf<DOFVector<T>, _unknown >(vector); }

  template<typename T>
  expressions::LaplacianOf<DOFVector<T>, _unknown > laplacianOf(DOFVector<T>* vector) 
  { return expressions::LaplacianOf<DOFVector<T>, _unknown >(vector); }

} // end namespace AMDiS


#endif // AMDIS_HESSIAN_OF_HPP
