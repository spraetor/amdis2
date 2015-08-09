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

#ifndef AMDIS_TRANSFORM_DOF_H
#define AMDIS_TRANSFORM_DOF_H

namespace AMDiS {

/**
 * analogon to std::for_each/std::transform for DOFVectors.
 * The operators are applied to each degree of freedom:
 * 
 * forEachDOF(vec, op)
 * transformDOF(vec1, result, op)
 * transformDOF(vec1, vec2, result, binary_op)
 * transformDOF(vec1, vec2, vec3, result, tertiary_op)
 *
 * transformDOFInterpolation(vec, f(x), result, binary_op)
 *
 * vec1,vec2,vec3 .. {DOFVector<Ti>*, DOFVector<Ti>&, Ti}
 * vec, result .. {DOFVector<TOut>*, DOFVector<TOut>&}
 * op .. AbstractFunction<TOut, T1>*
 * binary_op ..BinaryAbstractFunction<TOut, T1, T2>*
 * tertiary_op .. TertiaryAbstractFunction<TOut, T1, T2, T3>*
 *
 * analogon to std::accumulate for DOFVectors:
 * 
 * result = accumulateDOF_simple(vec, value0, binary_op)
 * result = accumulateDOF_simple(vec1, vec2, value0, tertiary_op)
 **/

// result = op(vec)
template<typename T1, typename T2>
inline void transformDOF_simple(DOFVector<T1> *vec,
  DOFVector<T2> *result,
  AbstractFunction<T2, T1> *op)
{
  TEST_EXIT(vec->getFeSpace() == result->getFeSpace())("FeSpaces must be equal!\n");
  DOFIterator<T1> vecIter(vec,USED_DOFS);
  DOFIterator<T2> resultIter(result,USED_DOFS);
  for(vecIter.reset(),resultIter.reset(); !vecIter.end(); ++vecIter,++resultIter)
  {
    *resultIter = (*op)(*vecIter);
  }
}

template<typename T1,typename T2>
inline void transformDOF_extended(DOFVector<T1> *vec, DOFVector<T2> *result, AbstractFunction<T2, T1> *op)
{
    DOFVector<T2>* res;
    bool useResult = true;
    if(static_cast<void*>(vec)==static_cast<void*>(result)) {
      res= new DOFVector<T2>(result->getFeSpace(), result->getName());
      useResult = false;
    } else
      res= result;

    const FiniteElemSpace *vecFeSpace = vec->getFeSpace();
    const FiniteElemSpace *resFeSpace = result->getFeSpace();

    const BasisFunction *vecBasisFcts = vecFeSpace->getBasisFcts();
    const BasisFunction *resBasisFcts = resFeSpace->getBasisFcts();

    int nVecBasisFcts = vecBasisFcts->getNumber();
    int nResBasisFcts = resBasisFcts->getNumber();

    std::vector<DegreeOfFreedom> resLocalIndices(nResBasisFcts);
    mtl::dense_vector<T1> vecLocalCoeffs(nVecBasisFcts);

    DimVec<double> *coords = NULL;
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(resFeSpace->getMesh(), -1,
            Mesh::CALL_LEAF_EL | 
            Mesh::FILL_COORDS);

    while (elInfo) {
      Element *el = elInfo->getElement();
      
      resBasisFcts->getLocalIndices(el, resFeSpace->getAdmin(), resLocalIndices);
      vec->getLocalVector(el, vecLocalCoeffs);
      
      for (int i = 0; i < nResBasisFcts; i++) {
        coords = resBasisFcts->getCoords(i);
        (*res)[resLocalIndices[i]] = (*op)(vecBasisFcts->evalUh(*coords, vecLocalCoeffs));
      }
      elInfo = stack.traverseNext(elInfo);
    }
    if (!useResult) {
      result->copy(*res);
      delete res;
    }
}

/// result = op(vec)
template<typename T1, typename T2>
inline void transformDOF(DOFVector<T1> *vec,
  DOFVector<T2> *result,
  AbstractFunction<T2, T1> *op)
{
  if (vec->getFeSpace() == result->getFeSpace())
    transformDOF_simple(vec,result,op);
  else
    transformDOF_extended(vec,result,op);
}

template<typename T1,typename T2> inline void transformDOF(DOFVector<T1> &vec, DOFVector<T2> &result, AbstractFunction<T2, T1> &op)
{ transformDOF(&vec, &result, &op); }

template<typename T1>
inline void transformDOF(DOFVector<T1> *vec, AbstractFunction<T1, T1> *op)
{ transformDOF(vec, vec, op); }

template<typename T1> inline void transformDOF(DOFVector<T1> &vec, AbstractFunction<T1, T1> &op)
{ transformDOF(&vec, &op); }

template<typename T1>
inline void forEachDOF(DOFVector<T1> *vec, AbstractFunction<T1, T1> *op)
{ transformDOF(vec, op); }

template<typename T1> inline void forEachDOF(DOFVector<T1> &vec, AbstractFunction<T1, T1> &op)
{ transformDOF(&vec, &op); }


// result = binary_op(vec1, vec2)
template<typename T1, typename T2, typename T3>
inline void transformDOF_simple(DOFVector<T1> *vec1,
  DOFVector<T2> *vec2, 
  DOFVector<T3> *result, 
  BinaryAbstractFunction<T3, T1, T2> *binary_op)
{
  TEST_EXIT(vec1->getFeSpace() == vec2->getFeSpace() && 
      vec2->getFeSpace() == result->getFeSpace())("FeSpaces must be equal!\n");
  DOFIterator<T1> vec1Iter(vec1,USED_DOFS);
  DOFIterator<T2> vec2Iter(vec2,USED_DOFS);
  DOFIterator<T3> resultIter(result,USED_DOFS);
  for(vec1Iter.reset(),vec2Iter.reset(),resultIter.reset(); !vec1Iter.end(); ++vec1Iter,++vec2Iter,++resultIter)
  {
    *resultIter = (*binary_op)(*vec1Iter, *vec2Iter);
  }
}
	
// result = binary_op(vec1, vec2)
template<typename T1, typename T2, typename T3>
inline void transformDOF_extended(DOFVector<T1> *vec1, DOFVector<T2> *vec2, DOFVector<T3> *result, BinaryAbstractFunction<T3, T1, T2> *binary_op)
{
    DOFVector<T3>* res;
    bool useResult = true;
    if(static_cast<void*>(vec1)==static_cast<void*>(result)
      || static_cast<void*>(vec2) == static_cast<void*>(result)) {
      res= new DOFVector<T3>(result->getFeSpace(), result->getName());
      useResult = false;
    } else
      res= result;

    const FiniteElemSpace *vec1FeSpace = vec1->getFeSpace();
    const FiniteElemSpace *vec2FeSpace = vec2->getFeSpace();
    const FiniteElemSpace *resFeSpace = result->getFeSpace();

    const BasisFunction *vec1BasisFcts = vec1FeSpace->getBasisFcts();
    const BasisFunction *vec2BasisFcts = vec2FeSpace->getBasisFcts();
    const BasisFunction *resBasisFcts = resFeSpace->getBasisFcts();

    int nVec1BasisFcts = vec1BasisFcts->getNumber();
    int nVec2BasisFcts = vec2BasisFcts->getNumber();
    int nResBasisFcts = resBasisFcts->getNumber();

    std::vector<DegreeOfFreedom> resLocalIndices(nResBasisFcts);
    mtl::dense_vector<T1> vec1LocalCoeffs(nVec1BasisFcts);
    mtl::dense_vector<T2> vec2LocalCoeffs(nVec2BasisFcts);

    DimVec<double> *coords = NULL;
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(resFeSpace->getMesh(), -1,
            Mesh::CALL_LEAF_EL | 
            Mesh::FILL_COORDS);

    while (elInfo) {
      Element *el = elInfo->getElement();
      
      resBasisFcts->getLocalIndices(el, resFeSpace->getAdmin(), resLocalIndices);
      vec1->getLocalVector(el, vec1LocalCoeffs);
      vec2->getLocalVector(el, vec2LocalCoeffs);
      
      for (int i = 0; i < nResBasisFcts; i++) {
        coords = resBasisFcts->getCoords(i);
        (*res)[resLocalIndices[i]] = (*binary_op)(
          vec1BasisFcts->evalUh(*coords, vec1LocalCoeffs),
          vec2BasisFcts->evalUh(*coords, vec2LocalCoeffs)
        );
      }
      elInfo = stack.traverseNext(elInfo);
    }
    if (!useResult) {
      result->copy(*res);
      delete res;
    }
}

/// result = binary_op(vec1, vec2)
template<typename T1, typename T2, typename T3>
inline void transformDOF(DOFVector<T1> *vec1,
  DOFVector<T2> *vec2,
  DOFVector<T3> *result,
  BinaryAbstractFunction<T3, T1, T2> *binary_op)
{
  if (vec1->getFeSpace() == vec2->getFeSpace() &&
      vec2->getFeSpace() == result->getFeSpace())
    transformDOF_simple(vec1,vec2,result,binary_op);
  else
    transformDOF_extended(vec1,vec2,result,binary_op);
}

template<typename T1, typename T2, typename T3> inline void transformDOF(DOFVector<T1> &vec1, DOFVector<T2> &vec2, DOFVector<T3> &result, BinaryAbstractFunction<T3, T1, T2> &binary_op)
{ transformDOF(&vec1, &vec2, &result, &binary_op); }


// result = binary_op(vec1, value)
template<typename T1, typename T2, typename T3>
inline void transformDOF_simple(DOFVector<T1> *vec1,
  const T2 val,
  DOFVector<T3> *result,
  BinaryAbstractFunction<T3, T1, T2> *binary_op)
{
  TEST_EXIT(vec1->getFeSpace() == result->getFeSpace())("FeSpaces must be equal!\n");
  DOFIterator<T1> vec1Iter(vec1,USED_DOFS);
  DOFIterator<T3> resultIter(result,USED_DOFS);
  for(vec1Iter.reset(),resultIter.reset(); !vec1Iter.end(); ++vec1Iter,++resultIter)
  {
    *resultIter = (*binary_op)(*vec1Iter, val);
  }
}

// result = binary_op(vec1, value)
template<typename T1, typename T2, typename T3>
inline void transformDOF_extended(DOFVector<T1> *vec1, const T2 val, DOFVector<T3> *result, BinaryAbstractFunction<T3, T1, T2> *binary_op)
{
    DOFVector<T3>* res;
    bool useResult = true;
    if(static_cast<void*>(vec1)==static_cast<void*>(result)) {
      res= new DOFVector<T3>(result->getFeSpace(), result->getName());
      useResult = false;
    } else
      res= result;

    const FiniteElemSpace *vec1FeSpace = vec1->getFeSpace();
    const FiniteElemSpace *resFeSpace = result->getFeSpace();

    const BasisFunction *vec1BasisFcts = vec1FeSpace->getBasisFcts();
    const BasisFunction *resBasisFcts = resFeSpace->getBasisFcts();

    int nVec1BasisFcts = vec1BasisFcts->getNumber();
    int nResBasisFcts = resBasisFcts->getNumber();

    std::vector<DegreeOfFreedom> resLocalIndices(nResBasisFcts);
    DenseVector<double> vec1LocalCoeffs(nVec1BasisFcts);

    DimVec<double> *coords = NULL;
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(resFeSpace->getMesh(), -1,
            Mesh::CALL_LEAF_EL | 
            Mesh::FILL_COORDS);

    while (elInfo) {
      Element *el = elInfo->getElement();
      
      resBasisFcts->getLocalIndices(el, resFeSpace->getAdmin(), resLocalIndices);
      vec1->getLocalVector(el, vec1LocalCoeffs);
      
      for (int i = 0; i < nResBasisFcts; i++) {
        coords = resBasisFcts->getCoords(i);
        (*res)[resLocalIndices[i]] = (*binary_op)(
          vec1BasisFcts->evalUh(*coords, vec1LocalCoeffs),
          val
        );
      }
      elInfo = stack.traverseNext(elInfo);
    }
    if (!useResult) {
      result->copy(*res);
      delete res;
    }
}

/// result = binary_op(vec1, value)
template<typename T1, typename T2, typename T3>
inline void transformDOF(DOFVector<T1> *vec1,
  const T2 val,
  DOFVector<T3> *result,
  BinaryAbstractFunction<T3, T1, T2> *binary_op)
{
  if (vec1->getFeSpace() == result->getFeSpace())
    transformDOF_simple(vec1,val,result,binary_op);
  else
    transformDOF_extended(vec1,val,result,binary_op);
}

template<typename T1, typename T2, typename T3> inline void transformDOF(DOFVector<T1> &vec1, const T2 val, DOFVector<T3> &result, BinaryAbstractFunction<T3, T1, T2> &binary_op)
{ transformDOF(&vec1, val, &result, &binary_op); }


// result = binary_op(value, vec2)
template<typename T1, typename T2, typename T3>
inline void transformDOF_simple(const T1 val,
  DOFVector<T2> *vec1,
  DOFVector<T3> *result,
  BinaryAbstractFunction<T3, T1, T2> *binary_op)
{
  TEST_EXIT(vec1->getFeSpace() == result->getFeSpace())("FeSpaces must be equal!\n");
  DOFIterator<T1> vec1Iter(vec1,USED_DOFS);
  DOFIterator<T3> resultIter(result,USED_DOFS);
  for(vec1Iter.reset(),resultIter.reset(); !vec1Iter.end(); ++vec1Iter,++resultIter)
  {
    *resultIter = (*binary_op)(val, *vec1Iter);
  }
}
	
// result = binary_op(value, vec2)
template<typename T1, typename T2, typename T3>
inline void transformDOF_extended(const T1 val, DOFVector<T2> *vec1, DOFVector<T3> *result, BinaryAbstractFunction<T3, T1, T2> *binary_op)
{
    DOFVector<T3>* res;
    bool useResult = true;
    if(static_cast<void*>(vec1)==static_cast<void*>(result)) {
      res= new DOFVector<T3>(result->getFeSpace(), result->getName());
      useResult = false;
    } else
      res= result;

    const FiniteElemSpace *vec1FeSpace = vec1->getFeSpace();
    const FiniteElemSpace *resFeSpace = result->getFeSpace();

    const BasisFunction *vec1BasisFcts = vec1FeSpace->getBasisFcts();
    const BasisFunction *resBasisFcts = resFeSpace->getBasisFcts();

    int nVec1BasisFcts = vec1BasisFcts->getNumber();
    int nResBasisFcts = resBasisFcts->getNumber();

    std::vector<DegreeOfFreedom> resLocalIndices(nResBasisFcts);
    DenseVector<double> vec1LocalCoeffs(nVec1BasisFcts);

    DimVec<double> *coords = NULL;
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(resFeSpace->getMesh(), -1,
            Mesh::CALL_LEAF_EL | 
            Mesh::FILL_COORDS);

    while (elInfo) {
      Element *el = elInfo->getElement();
      
      resBasisFcts->getLocalIndices(el, resFeSpace->getAdmin(), resLocalIndices);
      vec1->getLocalVector(el, vec1LocalCoeffs);
      
      for (int i = 0; i < nResBasisFcts; i++) {
        coords = resBasisFcts->getCoords(i);
        (*res)[resLocalIndices[i]] = (*binary_op)(
          val,
          vec1BasisFcts->evalUh(*coords, vec1LocalCoeffs)
        );
      }
      elInfo = stack.traverseNext(elInfo);
    }
    if (!useResult) {
      result->copy(*res);
      delete res;
    }
}

/// result = binary_op(value, vec2)
template<typename T1, typename T2, typename T3>
inline void transformDOF(const T1 val,
  DOFVector<T2> *vec1,
  DOFVector<T3> *result,
  BinaryAbstractFunction<T3, T1, T2> *binary_op)
{
  if (vec1->getFeSpace() == result->getFeSpace())
    transformDOF_simple(val,vec1,result,binary_op);
  else
    transformDOF_extended(val,vec1,result,binary_op);
}

template<typename T1, typename T2, typename T3> inline void transformDOF(const T1 val, DOFVector<T2> &vec1, DOFVector<T3> &result, BinaryAbstractFunction<T3, T1, T2> &binary_op)
{ transformDOF(val, &vec1, &result, &binary_op); }


// result = tertiary_op(vec1, vec2, vec3)
template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF_simple(DOFVector<T1> *vec1,
  DOFVector<T2> *vec2,
  DOFVector<T3> *vec3,
  DOFVector<T4> *result,
  TertiaryAbstractFunction<T4, T1, T2, T3> *tertiary_op)
{
  TEST_EXIT(vec1->getFeSpace() == vec2->getFeSpace() &&
      vec2->getFeSpace() == vec3->getFeSpace() &&
      vec3->getFeSpace() == result->getFeSpace())("FeSpaces must be equal!\n");
  DOFIterator<T1> vec1Iter(vec1,USED_DOFS);
  DOFIterator<T2> vec2Iter(vec2,USED_DOFS);
  DOFIterator<T3> vec3Iter(vec3,USED_DOFS);
  DOFIterator<T4> resultIter(result,USED_DOFS);
  for(vec1Iter.reset(),vec2Iter.reset(),vec3Iter.reset(),resultIter.reset(); !vec1Iter.end(); ++vec1Iter,++vec2Iter,++vec3Iter,++resultIter)
  {
    *resultIter = (*tertiary_op)(*vec1Iter, *vec2Iter, *vec3Iter);
  }
}

// result = tertiary_op(vec1, vec2, vec3)
template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF_extended(DOFVector<T1> *vec1, DOFVector<T2> *vec2, DOFVector<T3> *vec3, DOFVector<T4> *result, TertiaryAbstractFunction<T4, T1, T2, T3> *tertiary_op)
{
    DOFVector<T4>* res;
    bool useResult = true;
    if(static_cast<void*>(vec1)==static_cast<void*>(result)
      || static_cast<void*>(vec2) == static_cast<void*>(result)
      || static_cast<void*>(vec3) == static_cast<void*>(result)) {
      res= new DOFVector<T4>(result->getFeSpace(), result->getName());
      useResult = false;
    } else
      res= result;

    const FiniteElemSpace *vec1FeSpace = vec1->getFeSpace();
    const FiniteElemSpace *vec2FeSpace = vec2->getFeSpace();
    const FiniteElemSpace *vec3FeSpace = vec3->getFeSpace();
    const FiniteElemSpace *resFeSpace = res->getFeSpace();

    const BasisFunction *vec1BasisFcts = vec1FeSpace->getBasisFcts();
    const BasisFunction *vec2BasisFcts = vec2FeSpace->getBasisFcts();
    const BasisFunction *vec3BasisFcts = vec3FeSpace->getBasisFcts();
    const BasisFunction *resBasisFcts = resFeSpace->getBasisFcts();

    int nVec1BasisFcts = vec1BasisFcts->getNumber();
    int nVec2BasisFcts = vec2BasisFcts->getNumber();
    int nVec3BasisFcts = vec3BasisFcts->getNumber();
    int nResBasisFcts = resBasisFcts->getNumber();

    std::vector<DegreeOfFreedom> resLocalIndices(nResBasisFcts);
    mtl::dense_vector<T1> vec1LocalCoeffs(nVec1BasisFcts);
    mtl::dense_vector<T2> vec2LocalCoeffs(nVec2BasisFcts);
    mtl::dense_vector<T3> vec3LocalCoeffs(nVec3BasisFcts);

    DimVec<double> *coords = NULL;
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(resFeSpace->getMesh(), -1,
            Mesh::CALL_LEAF_EL | 
            Mesh::FILL_COORDS);

    while (elInfo) {
      Element *el = elInfo->getElement();
      
      resBasisFcts->getLocalIndices(el, resFeSpace->getAdmin(), resLocalIndices);
      vec1->getLocalVector(el, vec1LocalCoeffs);
      vec2->getLocalVector(el, vec2LocalCoeffs);
      vec3->getLocalVector(el, vec3LocalCoeffs);
      
      for (int i = 0; i < nResBasisFcts; i++) {
        coords = resBasisFcts->getCoords(i);
        (*res)[resLocalIndices[i]] = (*tertiary_op)(
          vec1BasisFcts->evalUh(*coords, vec1LocalCoeffs),
          vec2BasisFcts->evalUh(*coords, vec2LocalCoeffs),
          vec3BasisFcts->evalUh(*coords, vec3LocalCoeffs)
        );
      }
      elInfo = stack.traverseNext(elInfo);
    }
    
    if (!useResult) {
      result->copy(*res);
      delete res;
    }
}

/// result = tertiary_op(vec1, vec2, vec3)
template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF(DOFVector<T1> *vec1,
  DOFVector<T2> *vec2,
  DOFVector<T3> *vec3,
  DOFVector<T4> *result,
  TertiaryAbstractFunction<T4, T1, T2, T3> *tertiary_op)
{
  if (vec1->getFeSpace() == vec2->getFeSpace() &&
      vec2->getFeSpace() == vec3->getFeSpace() &&
      vec3->getFeSpace() == result->getFeSpace())
    transformDOF_simple(vec1,vec2,vec3,result,tertiary_op);
  else
    transformDOF_extended(vec1,vec2,vec3,result,tertiary_op);
}

template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF(DOFVector<T1> &vec1, DOFVector<T2> &vec2, DOFVector<T3> &vec3, DOFVector<T4> &result, TertiaryAbstractFunction<T4, T1, T2, T3> &tertiary_op)
{ transformDOF(&vec1, &vec2, &vec3, &result, &tertiary_op); }
	
/// result = tertiary_op(vec1, vec2, value)
template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF(DOFVector<T1> *vec1, DOFVector<T2> *vec2, T3 val, DOFVector<T4> *result, TertiaryAbstractFunction<T4, T1, T2, T3> *tertiary_op)
{

    DOFVector<T4>* res;
    bool useResult = true;
    if(static_cast<void*>(vec1)==static_cast<void*>(result)
      || static_cast<void*>(vec2) == static_cast<void*>(result)) {
      res= new DOFVector<T4>(result->getFeSpace(), result->getName());
      useResult = false;
    } else
      res= result;
    
    const FiniteElemSpace *vec1FeSpace = vec1->getFeSpace();
    const FiniteElemSpace *vec2FeSpace = vec2->getFeSpace();
    const FiniteElemSpace *resFeSpace = res->getFeSpace();

    const BasisFunction *vec1BasisFcts = vec1FeSpace->getBasisFcts();
    const BasisFunction *vec2BasisFcts = vec2FeSpace->getBasisFcts();
    const BasisFunction *resBasisFcts = resFeSpace->getBasisFcts();

    int nVec1BasisFcts = vec1BasisFcts->getNumber();
    int nVec2BasisFcts = vec2BasisFcts->getNumber();
    int nResBasisFcts = resBasisFcts->getNumber();

    std::vector<DegreeOfFreedom> resLocalIndices(nResBasisFcts);
    DenseVector<double> vec1LocalCoeffs(nVec1BasisFcts);
    DenseVector<double> vec2LocalCoeffs(nVec2BasisFcts);

    DimVec<double> *coords = NULL;
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(resFeSpace->getMesh(), -1,
            Mesh::CALL_LEAF_EL | 
            Mesh::FILL_COORDS);

    while (elInfo) {
      Element *el = elInfo->getElement();
      
      resBasisFcts->getLocalIndices(el, resFeSpace->getAdmin(), resLocalIndices);
      vec1->getLocalVector(el, vec1LocalCoeffs);
      vec2->getLocalVector(el, vec2LocalCoeffs);
      
      for (int i = 0; i < nResBasisFcts; i++) {
        coords = resBasisFcts->getCoords(i);
        (*res)[resLocalIndices[i]] = (*tertiary_op)(
          vec1BasisFcts->evalUh(*coords, vec1LocalCoeffs),
          vec2BasisFcts->evalUh(*coords, vec2LocalCoeffs),
          val
        );
      }
      elInfo = stack.traverseNext(elInfo);
    }

    if (!useResult) {
      result->copy(*res);
      delete res;
    }
}

template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF(DOFVector<T1> &vec1, DOFVector<T2> &vec2, T3 val, DOFVector<T4> &result, TertiaryAbstractFunction<T4, T1, T2, T3> &tertiary_op)
{ transformDOF(&vec1, &vec2, val, &result, &tertiary_op); }

/// result = tertiary_op(vec1, value, vec3)
template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF(DOFVector<T1> *vec1, T2 val, DOFVector<T2> *vec3, DOFVector<T4> *result, TertiaryAbstractFunction<T4, T1, T2, T3> *tertiary_op)
{

    DOFVector<T4>* res;
    bool useResult = true;
    if(static_cast<void*>(vec1)==static_cast<void*>(result)
      || static_cast<void*>(vec3) == static_cast<void*>(result)) {
      res= new DOFVector<T4>(result->getFeSpace(), result->getName());
      useResult = false;
    } else
      res= result;

    const FiniteElemSpace *vec1FeSpace = vec1->getFeSpace();
    const FiniteElemSpace *vec3FeSpace = vec3->getFeSpace();
    const FiniteElemSpace *resFeSpace = res->getFeSpace();

    const BasisFunction *vec1BasisFcts = vec1FeSpace->getBasisFcts();
    const BasisFunction *vec3BasisFcts = vec3FeSpace->getBasisFcts();
    const BasisFunction *resBasisFcts = resFeSpace->getBasisFcts();

    int nVec1BasisFcts = vec1BasisFcts->getNumber();
    int nVec3BasisFcts = vec3BasisFcts->getNumber();
    int nResBasisFcts = resBasisFcts->getNumber();

    std::vector<DegreeOfFreedom> resLocalIndices(nResBasisFcts);
    DenseVector<double> vec1LocalCoeffs(nVec1BasisFcts);
    DenseVector<double> vec3LocalCoeffs(nVec3BasisFcts);

    DimVec<double> *coords = NULL;
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(resFeSpace->getMesh(), -1,
            Mesh::CALL_LEAF_EL |
            Mesh::FILL_COORDS);

    while (elInfo) {
      Element *el = elInfo->getElement();

      resBasisFcts->getLocalIndices(el, resFeSpace->getAdmin(), resLocalIndices);
      vec1->getLocalVector(el, vec1LocalCoeffs);
      vec3->getLocalVector(el, vec3LocalCoeffs);

      for (int i = 0; i < nResBasisFcts; i++) {
        coords = resBasisFcts->getCoords(i);
        (*res)[resLocalIndices[i]] = (*tertiary_op)(
          vec1BasisFcts->evalUh(*coords, vec1LocalCoeffs),
          val,
          vec3BasisFcts->evalUh(*coords, vec3LocalCoeffs)
	  );
      }
      elInfo = stack.traverseNext(elInfo);
    }

    if (!useResult) {
      result->copy(*res);
      delete res;
    }
}

template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF(DOFVector<T1> &vec1, T2 val, DOFVector<T3> &vec3, DOFVector<T4> &result, TertiaryAbstractFunction<T4, T1, T2, T3> &tertiary_op)
{ transformDOF(&vec1, val, &vec3, &result, &tertiary_op); }

/// result = tertiary_op(value, vec2, vec3)
template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF(T1 val, DOFVector<T2> *vec2, DOFVector<T2> *vec3, DOFVector<T4> *result, TertiaryAbstractFunction<T4, T1, T2, T3> *tertiary_op)
{

    DOFVector<T4>* res;
    bool useResult = true;
    if(static_cast<void*>(vec2)==static_cast<void*>(result)
      || static_cast<void*>(vec3) == static_cast<void*>(result)) {
      res= new DOFVector<T4>(result->getFeSpace(), result->getName());
      useResult = false;
    } else
      res= result;

    const FiniteElemSpace *vec2FeSpace = vec2->getFeSpace();
    const FiniteElemSpace *vec3FeSpace = vec3->getFeSpace();
    const FiniteElemSpace *resFeSpace = res->getFeSpace();

    const BasisFunction *vec2BasisFcts = vec2FeSpace->getBasisFcts();
    const BasisFunction *vec3BasisFcts = vec3FeSpace->getBasisFcts();
    const BasisFunction *resBasisFcts = resFeSpace->getBasisFcts();

    int nVec2BasisFcts = vec2BasisFcts->getNumber();
    int nVec3BasisFcts = vec3BasisFcts->getNumber();
    int nResBasisFcts = resBasisFcts->getNumber();

    std::vector<DegreeOfFreedom> resLocalIndices(nResBasisFcts);
    DenseVector<double> vec2LocalCoeffs(nVec2BasisFcts);
    DenseVector<double> vec3LocalCoeffs(nVec3BasisFcts);

    DimVec<double> *coords = NULL;
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(resFeSpace->getMesh(), -1,
            Mesh::CALL_LEAF_EL |
            Mesh::FILL_COORDS);

    while (elInfo) {
      Element *el = elInfo->getElement();

      resBasisFcts->getLocalIndices(el, resFeSpace->getAdmin(), resLocalIndices);
      vec2->getLocalVector(el, vec2LocalCoeffs);
      vec3->getLocalVector(el, vec3LocalCoeffs);

      for (int i = 0; i < nResBasisFcts; i++) {
        coords = resBasisFcts->getCoords(i);
        (*res)[resLocalIndices[i]] = (*tertiary_op)(
          val,
          vec2BasisFcts->evalUh(*coords, vec2LocalCoeffs),
          vec3BasisFcts->evalUh(*coords, vec3LocalCoeffs)
	  );
      }
      elInfo = stack.traverseNext(elInfo);
    }

    if (!useResult) {
      result->copy(*res);
      delete res;
    }
}

template<typename T1, typename T2, typename T3, typename T4>
inline void transformDOF(T1 val, DOFVector<T2> &vec2, DOFVector<T3> &vec3, DOFVector<T4> &result, TertiaryAbstractFunction<T4, T1, T2, T3> &tertiary_op)
{ transformDOF(val, &vec2, &vec3, &result, &tertiary_op); }


// ===========================================================================================

/// return binary_op(vec, interpol(fct))
template<typename T>
inline void transformDOFInterpolation(
	DOFVector<T> *vec, 
	AbstractFunction<T, WorldVector<double> > *fct,
	BinaryAbstractFunction<T, T, T> *binary_op)
{
	DOFVector<T> helpDOF(vec->getFeSpace(), "temp");
	helpDOF.interpol(fct);
	transformDOF(vec,&helpDOF,vec,binary_op);
}

template<typename T> inline void transformDOFInterpolation(
	DOFVector<T> &vec, AbstractFunction<T, WorldVector<double> > &fct, BinaryAbstractFunction<T, T, T> &binary_op)
{ transformDOFInterpolation(&vec, &fct, &binary_op); }

/// return binary_op(vec, interpol(fct))
template<typename T>
inline void transformDOFInterpolation(
	DOFVector<T> *vec, 
	AbstractFunction<T, WorldVector<double> > *fct,
	DOFVector<T> *result, 
	BinaryAbstractFunction<T, T, T> *binary_op)
{
	DOFVector<T> helpDOF(vec->getFeSpace(), "temp");
	helpDOF.interpol(fct);
	transformDOF(vec,&helpDOF,result,binary_op);
}

template<typename T> inline void transformDOFInterpolation(
	DOFVector<T> &vec, AbstractFunction<T, WorldVector<double> > &fct, DOFVector<T> &result, BinaryAbstractFunction<T, T, T> &binary_op)
{ transformDOFInterpolation(&vec, &fct, &result, &binary_op); }


// ====================================================================================

template<typename S, typename T>
S accumulateDOF_simple(DOFVector<T> *vec,
		       S value0,
		       BinaryAbstractFunction<S, S, T> *binary_op)
{
  DOFIterator<T> vecIter(vec, USED_DOFS);
  S value = value0;
  for(vecIter.reset(); !vecIter.end(); ++vecIter)
  {
    value = (*binary_op)(value, *vecIter);
  }
  return value;
}

template<typename TOut, typename T1, typename T2>
TOut accumulateDOF_simple(DOFVector<T1> *vec1,
				 DOFVector<T2> *vec2,
				 TOut value0,
				 TertiaryAbstractFunction<TOut, TOut, T1, T2> *tertiary_op)
{
  TEST_EXIT(vec1->getFeSpace() == vec2->getFeSpace())("FeSpaces must be equal!\n");
  DOFIterator<T1> vec1Iter(vec1, USED_DOFS);
  DOFIterator<T2> vec2Iter(vec2, USED_DOFS);
  TOut value = value0;
  for(vec1Iter.reset(),vec2Iter.reset(); !vec1Iter.end(); ++vec1Iter,++vec2Iter)
  {
    value = (*tertiary_op)(value, *vec1Iter, *vec2Iter);
  }
  return value;
}

} // end namespace

#endif // AMDIS_TRANSFORM_DOF_H
