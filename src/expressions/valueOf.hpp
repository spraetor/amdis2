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



/** \file valueOf.hpp */

#ifndef AMDIS_VALUE_OF_HPP
#define AMDIS_VALUE_OF_HPP

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
#include "DOFVector.h"

#include "traits/category.hpp"
#include "traits/at.hpp"

#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>


namespace AMDiS 
{
  struct _unknown {};

  namespace expressions 
  {
    /// Expressions that extracts the values of a DOFVector at QPs
    template<typename Vector, typename Name, typename Enable = void>
    struct ValueOf : public LazyOperatorTermBase {};
    
    template<typename T, typename Name>
    struct ValueOf<DOFVector<T>, Name> : public LazyOperatorTermBase
    {
      typedef T value_type;
      typedef Name  id;

      DOFVector<T>* vecDV;
      mutable mtl::dense_vector<T> vec;
      mutable mtl::dense_vector<T> coeff;

      ValueOf(DOFVector<T>& vector) : vecDV(&vector) {}
      ValueOf(DOFVector<T>* vector) : vecDV(vector) {}

      template<typename List>
      inline void insertFeSpaces(List& feSpaces) const
      {
	feSpaces.insert(vecDV->getFeSpace());
      }
      
      inline int getDegree() const
      {
	return vecDV->getFeSpace()->getBasisFcts()->getDegree();
      }

      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* elInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {
	if (ot && subAssembler)
	  ot->getVectorAtQPs(vecDV, elInfo, subAssembler, quad, vec);
	else if (quad)
	  vecDV->getVecAtQPs(elInfo, quad, NULL, vec);
	else if (basisFct) {
	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
	  
	  // get coefficients of DOFVector
	  coeff.change_dim(localBasisFct->getNumber());
	  vecDV->getLocalVector(elInfo->getElement(), coeff);
	  
	  // eval basisfunctions of DOFVector at coords of given basisFct
	  size_t nBasisFct = basisFct->getNumber();
	  vec.change_dim(nBasisFct);	
	  for (size_t i = 0; i < nBasisFct; i++)
	    vec[i] = localBasisFct->evalUh(*basisFct->getCoords(i), coeff);
	}
      }


      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {
	if (ot && subAssembler)
	  ot->getVectorAtQPs(vecDV, smallElInfo, largeElInfo, subAssembler, quad, vec);
	else if (quad)
	  vecDV->getVecAtQPs(smallElInfo, largeElInfo, quad, NULL, vec);
	else if (basisFct) {
	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
	  
	  // get coefficients of DOFVector
	  coeff.change_dim(localBasisFct->getNumber());
	  vecDV->getLocalVector(smallElInfo->getElement(), coeff);
	  
	  // eval basisfunctions of DOFVector at coords of given basisFct
	  size_t nBasisFct = basisFct->getNumber();
	  vec.change_dim(nBasisFct);	
	  for (size_t i = 0; i < nBasisFct; i++)
	    vec[i] = localBasisFct->evalUh(*basisFct->getCoords(i), coeff);
	}
      }

      inline value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("value(") + vecDV->getName() + ")"; }
    };
    
    
    /// Expressions that extracts the matrix-value of a Matrix<DOFVector> at QPs
    template<template<class> class Matrix, typename T, typename Name>
    struct ValueOf<Matrix<DOFVector<T>*>, Name, 
		  typename enable_if< traits::is_matrix<Matrix<T> > >::type > 
      : public LazyOperatorTermBase
    {
      typedef Matrix<T> value_type;
      typedef Name  id;

      Matrix<DOFVector<T>*> vecDV;
      mutable mtl::dense_vector<value_type> vec;
      mutable Matrix<mtl::dense_vector<T> > coeff;

      ValueOf(Matrix<DOFVector<T>*> const& vector) : vecDV(vector) 
      {
	resize(coeff, num_rows(vecDV), num_cols(vecDV));
      }

      template<typename List>
      inline void insertFeSpaces(List& feSpaces) const
      {
	for (size_t i = 0; i < num_rows(vecDV); i++)
	for (size_t j = 0; j < num_cols(vecDV); j++)
	  feSpaces.insert(at(vecDV, i, j)->getFeSpace());
      }
      
      inline int getDegree() const
      {
	return at(vecDV, 0, 0)->getFeSpace()->getBasisFcts()->getDegree();
      }

      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* elInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {
	Matrix<mtl::dense_vector<T> > helper; resize(helper, num_rows(vecDV), num_cols(vecDV));
	for (size_t i = 0; i < num_rows(vecDV); i++) {
	for (size_t j = 0; j < num_cols(vecDV); j++) {
	  if (ot && subAssembler)
	    ot->getVectorAtQPs(at(vecDV, i, j), elInfo, subAssembler, quad, at(helper, i, j));
	  else if (quad)
	    at(vecDV, i, j)->getVecAtQPs(elInfo, quad, NULL, at(helper, i, j));	
	  else if (basisFct) {
	    const BasisFunction *localBasisFct = at(vecDV, i, j)->getFeSpace()->getBasisFcts();
	    
	    // get coefficients of DOFVector
	    at(coeff, i, j).change_dim(localBasisFct->getNumber());
	    at(vecDV, i, j)->getLocalVector(elInfo->getElement(), at(coeff, i, j));
	    
	    // eval basisfunctions of DOFVector at coords of given basisFct
	    size_t nBasisFct = basisFct->getNumber();
	    at(helper, i, j).change_dim(nBasisFct);	
	    for (size_t k = 0; k < nBasisFct; k++)
	      at(helper, i, j)[k] = localBasisFct->evalUh(*basisFct->getCoords(k), at(coeff, i, j));
	  }
	}
	}
	vec.change_dim(num_rows(at(helper, 0, 0)));
	for (size_t iq = 0; iq < num_rows(at(helper, 0, 0)); iq++) {
	  value_type tmp; resize(tmp, num_rows(vecDV), num_cols(vecDV));
	  for (size_t i = 0; i < num_rows(vecDV); i++) {
	  for (size_t j = 0; j < num_cols(vecDV); j++) {
	    at(tmp, i, j) = at(helper, i, j)[iq];
	  }
	  }
	  vec[iq] = tmp;
	}
      }


      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {      
	initElement(ot, smallElInfo, subAssembler, quad, basisFct);
      }

      inline value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("value_(") + at(vecDV, 0, 0)->getName() + ")"; }
    };
    
    
    /// Expressions that extracts the vector-value of a Vector<DOFVector> at QPs
    template<template<class> class Vector, typename T, typename Name>
    struct ValueOf<Vector<DOFVector<T>*>, Name, 
		  typename boost::enable_if<typename traits::is_vector<Vector<T> >::type>::type > 
      : public LazyOperatorTermBase
    {
      typedef Vector<T> value_type;
      typedef Name  id;

      Vector<DOFVector<T>*> vecDV;
      mutable mtl::dense_vector<value_type> vec;
      mutable Vector<mtl::dense_vector<T> > coeff;

      ValueOf(Vector<DOFVector<T>*>& vector) : vecDV(vector) 
      {
	resize(coeff, num_rows(vecDV));
      }

      template<typename List>
      inline void insertFeSpaces(List& feSpaces) const
      {
	for (size_t i = 0; i < num_rows(vecDV); i++)
	  feSpaces.insert(at(vecDV, i)->getFeSpace());
      }
      
      inline int getDegree() const
      {
	return at(vecDV, 0)->getFeSpace()->getBasisFcts()->getDegree();
      }

      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* elInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {
	Vector<mtl::dense_vector<T> > helper; resize(helper, num_rows(vecDV));
	for (size_t i = 0; i < num_rows(vecDV); i++) {
	  if (ot && subAssembler)
	    ot->getVectorAtQPs(at(vecDV, i), elInfo, subAssembler, quad, at(helper, i));
	  else if (quad)
	    at(vecDV, i)->getVecAtQPs(elInfo, quad, NULL, at(helper, i));	
	  else if (basisFct) {
	    const BasisFunction *localBasisFct = at(vecDV, i)->getFeSpace()->getBasisFcts();
	    
	    // get coefficients of DOFVector
	    at(coeff, i).change_dim(localBasisFct->getNumber());
	    at(vecDV, i)->getLocalVector(elInfo->getElement(), at(coeff, i));
	    
	    // eval basisfunctions of DOFVector at coords of given basisFct
	    size_t nBasisFct = basisFct->getNumber();
	    at(helper, i).change_dim(nBasisFct);	
	    for (size_t j = 0; j < nBasisFct; j++)
	      at(helper, i)[j] = localBasisFct->evalUh(*basisFct->getCoords(j), at(coeff, i));
	  }
	}
	vec.change_dim(num_rows(at(helper, 0)));
	for (size_t iq = 0; iq < num_rows(at(helper, 0)); iq++) {
	  value_type tmp; resize(tmp, num_rows(vecDV));
	  for (size_t i = 0; i < num_rows(vecDV); i++)
	    at(tmp, i) = at(helper, i)[iq];
	  vec[iq] = tmp;
	}
      }


      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {      
	initElement(ot, smallElInfo, subAssembler, quad, basisFct);
      }

      inline value_type operator()(const int& iq) const { return vec[iq]; }
      
      std::string str() const { return std::string("value_(") + at(vecDV, 0)->getName() + ")"; }
    };
    
    
    /// Expression that extracts the component of a vector-values DOFVector at QPs
    template<typename Vector>
    struct ComponentOf : public LazyOperatorTermBase {};
    
    template<template<class> class Vector, typename T>
    struct ComponentOf<DOFVector<Vector<T> > > : public LazyOperatorTermBase
    {
      typedef T value_type;

      DOFVector<Vector<T> >* vecDV;
      mutable mtl::dense_vector<Vector<T> > vec;
      mutable mtl::dense_vector<Vector<T> > coeff;
      int I;

      ComponentOf(DOFVector<Vector<T> >& vector, int I_) : vecDV(&vector), I(I_) {}
      ComponentOf(DOFVector<Vector<T> >* vector, int I_) : vecDV(vector), I(I_) {}

      template<typename List>
      inline void insertFeSpaces(List& feSpaces) const
      {
	feSpaces.insert(vecDV->getFeSpace());
      }
      
      inline int getDegree() const
      {
	return vecDV->getFeSpace()->getBasisFcts()->getDegree();
      }

      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* elInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {
	if (ot && subAssembler)
	  ot->getVectorAtQPs(vecDV, elInfo, subAssembler, quad, vec);
	else if (quad)
	  vecDV->getVecAtQPs(elInfo, quad, NULL, vec);
	else if (basisFct) {
	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
	  
	  // get coefficients of DOFVector
	  coeff.change_dim(localBasisFct->getNumber());
	  vecDV->getLocalVector(elInfo->getElement(), coeff);
	  
	  // eval basisfunctions of DOFVector at coords of given basisFct
	  size_t nBasisFct = basisFct->getNumber();
	  vec.change_dim(nBasisFct);	
	  for (size_t i = 0; i < nBasisFct; i++)
	    vec[i] = localBasisFct->evalUh(*basisFct->getCoords(i), coeff);
	}
      }


      template<typename OT>
      inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
			      SubAssembler* subAssembler, Quadrature *quad, 
			      const BasisFunction *basisFct = NULL)
      {
	if (ot && subAssembler)
	  ot->getVectorAtQPs(vecDV, smallElInfo, largeElInfo, subAssembler, quad, vec);
	else if (quad)
	  vecDV->getVecAtQPs(smallElInfo, largeElInfo, quad, NULL, vec);
	else if (basisFct) {
	  const BasisFunction *localBasisFct = vecDV->getFeSpace()->getBasisFcts();
	  
	  // get coefficients of DOFVector
	  coeff.change_dim(localBasisFct->getNumber());
	  vecDV->getLocalVector(smallElInfo->getElement(), coeff);
	  
	  // eval basisfunctions of DOFVector at coords of given basisFct
	  size_t nBasisFct = basisFct->getNumber();
	  vec.change_dim(nBasisFct);	
	  for (size_t i = 0; i < nBasisFct; i++)
	    vec[i] = localBasisFct->evalUh(*basisFct->getCoords(i), coeff);
	}
      }

      inline value_type operator()(const int& iq) const { return vec[iq][I]; }
      
      std::string str() const { return std::string("comp<") + boost::lexical_cast<std::string>(I) + ">(" + vecDV->getName() + ")"; }
    };


  } // end namespace expressions

  // value of a DOFVector<T>
  // ___________________________________________________________________________

  // with Name
  template<typename Name, typename T>
  expressions::ValueOf<DOFVector<T>, Name > valueOf(DOFVector<T>& vector) 
  { return expressions::ValueOf<DOFVector<T>, Name >(vector); }

  template<typename Name, typename T>
  expressions::ValueOf<DOFVector<T>, Name > valueOf(DOFVector<T>* vector) 
  { return expressions::ValueOf<DOFVector<T>, Name >(vector); }

  template<typename Name, template<class> class Matrix, typename T>
  typename boost::enable_if<typename traits::is_matrix<Matrix<T> >::type, 
    expressions::ValueOf<Matrix<DOFVector<T>*>, Name > >::type
  valueOf(Matrix<DOFVector<T>*> &mat) 
  { return expressions::ValueOf<Matrix<DOFVector<T>*>, Name >(mat); }

  template<typename Name, template<class> class Vector, typename T>
  typename boost::enable_if<typename traits::is_vector<Vector<T> >::type, 
    expressions::ValueOf<Vector<DOFVector<T>*>, Name > >::type
  valueOf(Vector<DOFVector<T>*> &vector) 
  { return expressions::ValueOf<Vector<DOFVector<T>*>, Name >(vector); }


  // without Name
  template<typename T>
  expressions::ValueOf<DOFVector<T>, _unknown > valueOf(DOFVector<T>& vector) 
  { return expressions::ValueOf<DOFVector<T>, _unknown >(vector); }

  template<typename T>
  expressions::ValueOf<DOFVector<T>, _unknown > valueOf(DOFVector<T>* vector) 
  { return expressions::ValueOf<DOFVector<T>, _unknown >(vector); }

  template<template<class> class Matrix, typename T>
  typename boost::enable_if<typename traits::is_matrix<Matrix<T> >::type, 
    expressions::ValueOf<Matrix<DOFVector<T>*>, _unknown > >::type
  valueOf(Matrix<DOFVector<T>*> &mat) 
  { return expressions::ValueOf<Matrix<DOFVector<T>*>, _unknown >(mat); }

  template<template<class> class Vector, typename T>
  typename boost::enable_if<typename traits::is_vector<Vector<T> >::type, 
    expressions::ValueOf<Vector<DOFVector<T>*>, _unknown > >::type
  valueOf(Vector<DOFVector<T>*> &vector) 
  { return expressions::ValueOf<Vector<DOFVector<T>*>, _unknown >(vector); }

  
  // component of a DOFVector<Vector>
  // ___________________________________________________________________________

  template<template<class> class Vector, typename T>
  expressions::ComponentOf<DOFVector<Vector<T> > > componentOf(DOFVector<Vector<T> >& vector, int I) 
  { return expressions::ComponentOf<DOFVector<Vector<T> > >(vector, I); }

  template<template<class> class Vector, typename T>
  expressions::ComponentOf<DOFVector<Vector<T> > > componentOf(DOFVector<Vector<T> >* vector, int I) 
  { return expressions::ComponentOf<DOFVector<Vector<T> > >(vector, I); }

} // end namespace AMDiS

#endif // AMDIS_VALUE_OF_HPP
