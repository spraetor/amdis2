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



/** \file gradientOf.hpp */

#ifndef AMDIS_GRADIENT_OF_HPP
#define AMDIS_GRADIENT_OF_HPP

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
    struct GradientOf : public LazyOperatorTermBase
    {
      typedef typename traits::category<Vector>::value_type T;
      typedef typename GradientType<T>::type value_type;
      typedef Name  id;

      DOFVector<T>* vecDV;
      mutable mtl::dense_vector<typename GradientType<T>::type> vec;
      mutable mtl::dense_vector<T> coeff;

      GradientOf(Vector& vector) : vecDV(&vector) {}
      GradientOf(Vector* vector) : vecDV(vector) {}

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
      {      
	if (ot && subAssembler)
	  ot->getGradientsAtQPs(vecDV, elInfo, subAssembler, quad, vec);
	else if (quad)
	  vecDV->getGrdAtQPs(elInfo, quad, NULL, vec);
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
	    localBasisFct->evalGrdUh(*basisFct->getCoords(i), grdLambda, coeff, vec[i]);
	}
      }


      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {
	if (ot && subAssembler)
	  ot->getGradientsAtQPs(vecDV, smallElInfo, largeElInfo, subAssembler, quad, vec);
	else if (quad)
	  vecDV->getGrdAtQPs(smallElInfo, largeElInfo, quad, NULL, vec);
	else if (basisFct) {
	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
	  
	  // get coefficients of DOFVector
	  coeff.change_dim(localBasisFct->getNumber());
	  vecDV->getLocalVector(smallElInfo->getElement(), coeff);
	  
	  // eval basisfunctions of DOFVector at coords of given basisFct
	  size_t nBasisFct = basisFct->getNumber();
	  vec.change_dim(nBasisFct);	
	  
	  const DimVec<WorldVector<double> > &grdLambda = smallElInfo->getGrdLambda();
	  for (size_t i = 0; i < nBasisFct; i++)
	    localBasisFct->evalGrdUh(*basisFct->getCoords(i), grdLambda, coeff, vec[i]);
	}
      }

      value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("grad(") + vecDV->getName() + ")"; }
    };
      
    
    /// Expressions that extracts the partial derivative of a DOFVector at QPs
    template<int I, typename Vector, typename Name>
    struct DerivativeOf : public LazyOperatorTermBase
    {
      typedef typename traits::category<Vector>::value_type T;
      typedef T     value_type;
      typedef Name  id;

      DOFVector<T>* vecDV;
      mutable mtl::dense_vector<typename GradientType<T>::type> vec;
  //     mutable mtl::dense_vector<T> vec;
      mutable mtl::dense_vector<T> coeff;
      int comp;

      DerivativeOf(Vector& vector) : vecDV(&vector), comp(I) {}
      DerivativeOf(Vector* vector) : vecDV(vector), comp(I) {}
      
      DerivativeOf(Vector& vector, int I0) : vecDV(&vector), comp(I0) 
      {
	TEST_EXIT_DBG( I < 0 && I0 >= 0 ) 
	  ("You yould specify eather template<int I>, or constructor(int I0)\n");
      }
      DerivativeOf(Vector* vector, int I0) : vecDV(vector), comp(I0) 
      {
	TEST_EXIT_DBG( I < 0 && I0 >= 0 )
	  ("You yould specify eather template<int I>, or constructor(int I0)\n");
      }

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
      {      
	if (ot && subAssembler)
	  ot->getGradientsAtQPs(vecDV, elInfo, subAssembler, quad, vec); //subAssembler->getDerivativeAtQPs(vecDV, elInfo, quad, comp, vec);
	else if (quad)
	  vecDV->getGrdAtQPs(elInfo, quad, NULL, vec); //vecDV->getDerivativeAtQPs(elInfo, quad, NULL, comp, vec);
	else if (basisFct) {
	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
	  
	  // get coefficients of DOFVector
	  coeff.change_dim(localBasisFct->getNumber());
	  vecDV->getLocalVector(elInfo->getElement(), coeff);
	  
	  // eval basisfunctions of DOFVector at coords of given basisFct
	  size_t nBasisFct = basisFct->getNumber();
  // 	mtl::dense_vector<typename GradientType<T>::type> helper(nBasisFct);	
	  
	  const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
	  vec.change_dim(nBasisFct);
	  for (size_t i = 0; i < nBasisFct; i++)
	    localBasisFct->evalGrdUh(*basisFct->getCoords(i), grdLambda, coeff, vec[i]); //helper[i]);
	  
  // 	for (size_t i = 0; i < num_rows(helper); i++)
  // 	  vec[i] = helper[i][comp];
	}
      }

      template<typename OT>
      void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
		      SubAssembler* subAssembler, Quadrature *quad, 
		      const BasisFunction *basisFct = NULL)
      {
  //       if (op && subAssembler)
  // 	ot->getGradientsAtQPs(vecDV, smallElInfo, largeElInfo, subAssembler, quad, vec);
  //       else
  // 	vecDV->getGrdAtQPs(smallElInfo, largeElInfo, localQuad, NULL, vec);
	
	if (ot && subAssembler)
	  ot->getGradientsAtQPs(vecDV, smallElInfo, largeElInfo, subAssembler, quad, vec); //subAssembler->getDerivativeAtQPs(vecDV, smallElInfo, largeElInfo, quad, comp, vec);
	else if (quad)
	  vecDV->getGrdAtQPs(smallElInfo, largeElInfo, quad, NULL, vec); // vecDV->getDerivativeAtQPs(smallElInfo, largeElInfo, quad, NULL, comp, vec);
	else if (basisFct) {
	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
	  
	  // get coefficients of DOFVector
	  coeff.change_dim(localBasisFct->getNumber());
	  vecDV->getLocalVector(smallElInfo->getElement(), coeff);
	  
	  // eval basisfunctions of DOFVector at coords of given basisFct
	  size_t nBasisFct = basisFct->getNumber();
  // 	mtl::dense_vector<typename GradientType<T>::type> helper(nBasisFct);	
	  
	  const DimVec<WorldVector<double> > &grdLambda = smallElInfo->getGrdLambda();
	  vec.change_dim(nBasisFct);
	  for (size_t i = 0; i < nBasisFct; i++)
	    localBasisFct->evalGrdUh(*basisFct->getCoords(i), grdLambda, coeff, vec[i]); //helper[i]);
	  
  // 	for (size_t i = 0; i < num_rows(helper); i++)
  // 	  vec[i] = helper[i][comp];
	}
      }

      value_type operator()(const int& iq) const { return vec[iq][comp]; }
  //     value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("deriv<") + boost::lexical_cast<std::string>(I) + ">(" + vecDV->getName() + ")"; }
    };
    
    
  #if 0
    /// Expressions that extracts the divergence of a DOFVector<Vector> at QPs
    template<typename Vector>
    struct DivergenceOf : public LazyOperatorTermBase
    {
      typedef typename traits::ValueType<Vector>::type T;      // e.g. WorldVector<double>
      typedef typename traits::ValueType<T>::type value_type;  // => double

      DOFVector<T>* vecDV;
      mutable mtl::dense_vector<typename GradientType<T>::type> vec;
      mutable mtl::dense_vector<T> coeff;

      GradientOf(Vector& vector) : vecDV(&vector) {}
      GradientOf(Vector* vector) : vecDV(vector) {}

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
      {      
	if (ot && subAssembler)
	  ot->getGradientsAtQPs(vecDV, elInfo, subAssembler, quad, vec);
	else if (quad)
	  vecDV->getGrdAtQPs(elInfo, quad, NULL, vec);
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
	    localBasisFct->evalGrdUh(*basisFct->getCoords(i), grdLambda, coeff, vec[i]);
	}
      }


      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {
	if (ot && subAssembler)
	  ot->getGradientsAtQPs(vecDV, smallElInfo, largeElInfo, subAssembler, quad, vec);
	else if (quad)
	  vecDV->getGrdAtQPs(smallElInfo, largeElInfo, quad, NULL, vec);
	else if (basisFct) {
	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
	  
	  // get coefficients of DOFVector
	  coeff.change_dim(localBasisFct->getNumber());
	  vecDV->getLocalVector(smallElInfo->getElement(), coeff);
	  
	  // eval basisfunctions of DOFVector at coords of given basisFct
	  size_t nBasisFct = basisFct->getNumber();
	  vec.change_dim(nBasisFct);	
	  
	  const DimVec<WorldVector<double> > &grdLambda = smallElInfo->getGrdLambda();
	  for (size_t i = 0; i < nBasisFct; i++)
	    localBasisFct->evalGrdUh(*basisFct->getCoords(i), grdLambda, coeff, vec[i]);
	}
      }

      value_type operator()(const int& iq) const { return vec[iq]; }
    };
      
  #endif
    
  } // end namespace expressions

  
  // gradient of a DOFVector
  // _____________________________________________________________________________

  // with Name
  template<typename Name, typename T>
  expressions::GradientOf<DOFVector<T>, Name > gradientOf(DOFVector<T>& vector) 
  { return expressions::GradientOf<DOFVector<T>, Name >(vector); }

  template<typename Name, typename T>
  expressions::GradientOf<DOFVector<T>, Name > gradientOf(DOFVector<T>* vector) 
  { return expressions::GradientOf<DOFVector<T>, Name >(vector); }

  // without Name
  template<typename T>
  expressions::GradientOf<DOFVector<T>, _unknown > gradientOf(DOFVector<T>& vector) 
  { return expressions::GradientOf<DOFVector<T>, _unknown >(vector); }

  template<typename T>
  expressions::GradientOf<DOFVector<T>, _unknown > gradientOf(DOFVector<T>* vector) 
  { return expressions::GradientOf<DOFVector<T>, _unknown >(vector); }


  // Partial derivative of a DOFVector
  // _____________________________________________________________________________

  // with Name
  template<typename Name, int I, typename T>
  expressions::DerivativeOf<I, DOFVector<T>, Name > derivativeOf(DOFVector<T>& vector) 
  { return expressions::DerivativeOf<I, DOFVector<T>, Name >(vector); }

  template<typename Name, int I, typename T>
  expressions::DerivativeOf<I, DOFVector<T>, Name > derivativeOf(DOFVector<T>* vector) 
  { return expressions::DerivativeOf<I, DOFVector<T>, Name >(vector); }

  template<typename Name, typename T>
  expressions::DerivativeOf<-1, DOFVector<T>, Name > derivativeOf(DOFVector<T>& vector, int I0) 
  { return expressions::DerivativeOf<-1, DOFVector<T>, Name >(vector, I0); }

  template<typename Name, typename T>
  expressions::DerivativeOf<-1, DOFVector<T>, Name > derivativeOf(DOFVector<T>* vector, int I0) 
  { return expressions::DerivativeOf<-1, DOFVector<T>, Name >(vector, I0); }


  // without Name
  template<int I, typename T>
  expressions::DerivativeOf<I, DOFVector<T>, _unknown > derivativeOf(DOFVector<T>& vector) 
  { return expressions::DerivativeOf<I, DOFVector<T>, _unknown >(vector); }

  template<int I, typename T>
  expressions::DerivativeOf<I, DOFVector<T>, _unknown > derivativeOf(DOFVector<T>* vector) 
  { return expressions::DerivativeOf<I, DOFVector<T>, _unknown >(vector); }

  template<typename T>
  expressions::DerivativeOf<-1, DOFVector<T>, _unknown > derivativeOf(DOFVector<T>& vector, int I0) 
  { return expressions::DerivativeOf<-1, DOFVector<T>, _unknown >(vector, I0); }

  template<typename T>
  expressions::DerivativeOf<-1, DOFVector<T>, _unknown > derivativeOf(DOFVector<T>* vector, int I0) 
  { return expressions::DerivativeOf<-1, DOFVector<T>, _unknown >(vector, I0); }

} // end namespace AMDiS


// ------- something special needed for gradientOf(DOFVector<WorldVector>)
#include <boost/numeric/mtl/operation/mult_result.hpp>
namespace mtl {
  namespace traits {
  
    typedef AMDiS::WorldVector<AMDiS::WorldVector<double> > WWMatrix;
  
    template <typename Op1, typename IsMatrix>
    struct mult_result_WWMatrix {
	typedef mat_cvec_times_expr<Op1, mtl::dense_vector<WWMatrix> > type;
    };
  
    /// Multiply matrix with column vector
    template <typename Op1>
    struct mult_result<Op1, mtl::dense_vector<WWMatrix> > 
      : public mult_result_WWMatrix<Op1, typename boost::enable_if<is_matrix<Op1> >::type >
    {};
  }
}

#endif // AMDIS_GRADIENT_OF_HPP
