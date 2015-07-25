/** \file gradientOf.hpp */

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

      void initElement(const ElInfo* elInfo,
		      SubAssembler* subAssembler, Quadrature *quad, 
		      const BasisFunction *basisFct = NULL)
      {      
      	if (subAssembler)
      	  subAssembler->getGradientsAtQPs(vecDV, elInfo, quad, vec);
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

      value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("grad(") + vecDV->getName() + ")"; }
    };
      
    
    /// Expressions that extracts the partial derivative of a DOFVector at QPs
    template <int I, class Vector, class Name>
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
      {      
      	if (subAssembler) // TODO: use specialization for derivative instead of gradient!!!
      	  subAssembler->getGradientsAtQPs(vecDV, elInfo, quad, vec); //subAssembler->getDerivativeAtQPs(vecDV, elInfo, quad, comp, vec);
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

      value_type operator()(const int& iq) const { return vec[iq][comp]; }
  //     value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("deriv<") + std::to_string(I) + ">(" + vecDV->getName() + ")"; }
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

      void initElement(const ElInfo* elInfo,
		      SubAssembler* subAssembler, Quadrature *quad, 
		      const BasisFunction *basisFct = NULL)
      {      
      	// if (ot && subAssembler)
      	  // ot->getGradientsAtQPs(vecDV, elInfo, subAssembler, quad, vec);
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

      value_type operator()(const int& iq) const { return vec[iq]; }
    };
      
  #endif
    
  } // end namespace expressions

  
  // gradient of a DOFVector
  // _____________________________________________________________________________

  // with Name
  template <class Name = _unknown, class T>
  expressions::GradientOf<DOFVector<T>, Name > 
  gradientOf(DOFVector<T>& vector) { return {vector}; }

  template <class Name = _unknown, class T>
  expressions::GradientOf<DOFVector<T>, Name > 
  gradientOf(DOFVector<T>* vector) { return {vector}; }

  // Partial derivative of a DOFVector
  // _____________________________________________________________________________

  // with Name
  template <class Name = _unknown, int I, class T>
  expressions::DerivativeOf<I, DOFVector<T>, Name > 
  derivativeOf(DOFVector<T>& vector) { return {vector}; }

  template <class Name = _unknown, int I, class T>
  expressions::DerivativeOf<I, DOFVector<T>, Name > 
  derivativeOf(DOFVector<T>* vector) { return {vector}; }

  template <class Name = _unknown, class T>
  expressions::DerivativeOf<-1, DOFVector<T>, Name > 
  derivativeOf(DOFVector<T>& vector, int I0) { return {vector, I0}; }

  template <class Name = _unknown, class T>
  expressions::DerivativeOf<-1, DOFVector<T>, Name > 
  derivativeOf(DOFVector<T>* vector, int I0) { return {vector, I0}; }

} // end namespace AMDiS


// ------- something special needed for gradientOf(DOFVector<WorldVector>)
#include <boost/numeric/mtl/operation/mult_result.hpp>
namespace mtl 
{
  namespace traits 
  {
  
    typedef AMDiS::WorldVector<AMDiS::WorldVector<double> > WWMatrix;
  
    template <class Op1, class IsMatrix>
    struct mult_result_WWMatrix {
      typedef mat_cvec_times_expr<Op1, mtl::dense_vector<WWMatrix> > type;
    };
  
    /// Multiply matrix with column vector
    template <class Op1>
    struct mult_result<Op1, mtl::dense_vector<WWMatrix> > 
      : public mult_result_WWMatrix<Op1, typename boost::enable_if<is_matrix<Op1> >::type >
    {};
    
  } // end namespace traits
  
} // end namespace mtl
