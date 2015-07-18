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



/** \file ElementFunction.h */

#ifndef AMDIS_ELEMENTFUNCTION_H
#define AMDIS_ELEMENTFUNCTION_H

#include "DOFVector.h"
#include "BasisFunction.h"

namespace AMDiS {

  /// Abstract access to functions living on elements.
  template<typename T>
  class ElementFunction // : public AbstractFunction<T, DimVec<double> >
  {
  public:
    /// constructor.
    ElementFunction() 
      : elInfo(NULL) 
    {}

    /// destructor.
    virtual ~ElementFunction() {}

    /// sets \ref elInfo_;
    inline void setElInfo(const ElInfo *elInfo_) 
    { 
      elInfo = elInfo_; 
    }

  protected:
    /// ElInfo the function currently lives on.
    const ElInfo *elInfo;
  };

  /// ElementFunction wich encapsulates the evaluation of an analytical function.
  template<typename T>
  class ElementFunctionAnalytic : public ElementFunction<T>
  {
  public:
    /// constructor
    ElementFunctionAnalytic(const std::function<T(WorldVector<double>)> &fct_)
      : ElementFunction<T>(),
	fct(fct_)
    {}

    /// evaluation at given coordinates.
    T operator()(const DimVec<double>& bary) const 
    {
      WorldVector<double> worldCoords;
      this->elInfo->coordToWorld(bary, worldCoords);
      return fct(worldCoords);
    }

  protected:
    /// function to be avaluated at world coordinates.
    const std::function<T(WorldVector<double>)> &fct;
  };


  /// ElementFunction wich encapsulates the interpolation of an DOFVector.
  template<typename T>
  class ElementFunctionDOFVec : public ElementFunction<T>
  {
  public:
    /// constructor.
    ElementFunctionDOFVec(const DOFVector<T> *vec) 
      : ElementFunction<T>(),
	dofVector(vec)
    {}

    /// evaluation at given coordinates.
    T operator()(const DimVec<double>& bary) const 
    {
      mtl::dense_vector<T> localVec(dofVector->getFeSpace()->getBasisFcts()->getNumber());
      dofVector->getLocalVector(this->elInfo->getElement(), localVec);
      
      return dofVector->getFeSpace()->getBasisFcts()->evalUh(bary, localVec);
    }

  protected:
    /// DOFVector to be interpolated.
    const DOFVector<T> *dofVector;
  };
  
  /// Abstract access to functions living on elements, but using world coordinates
  template<typename T>
  class ElementFunctionWorld // : public AbstractFunction<T, WorldVector<double> >
  {
  public:
    /// constructor.
    ElementFunctionWorld(const DOFVector<T>* vec_, double factor_ = 1.0) 
      : elInfo(NULL),
	vec(vec_),
	localCoeff(vec_->getFeSpace()->getBasisFcts()->getNumber()),
	coords(vec_->getFeSpace()->getBasisFcts()->getDim()),
	factor(factor_)
    {}

    /// destructor.
    virtual ~ElementFunctionWorld() {}

    /// evaluation at given coordinates.
    T operator()(const WorldVector<double>& x) const 
    {
      const BasisFunction *basisFcts = vec->getFeSpace()->getBasisFcts();
      
      elInfo->worldToCoord(x, &coords);
      return basisFcts->evalUh(coords, localCoeff) * factor;
    }
    
    /// sets \ref elInfo_;
    inline void setElInfo(const ElInfo *elInfo_) 
    {
      elInfo = elInfo_; 
      vec->getLocalVector(elInfo->getElement(), localCoeff);
    }

  protected:
    /// ElInfo the function currently lives on.
    const ElInfo *elInfo;
    
    /// DOFVector to be interpolated.
    const DOFVector<T>* vec;
    
    /// local dof-values of DOFVector on current element
    mtl::dense_vector<T> localCoeff;
    
    /// barycentric coordinates
    mutable DimVec<double> coords;
    
    /// global factor multiplied with result of operator()
    double factor;
  };

}

#endif
