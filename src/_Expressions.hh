/** \file Expressions.hh */

#include "Traverse.h"
#include "DualTraverse.h"

#include <traits/basic.hpp>

namespace AMDiS 
{
  template <class T, class Term>
  inline requires::Term< void, Term >
  transformDOF(Term term, DOFVector<T>* result)
  {
    typedef Value_t<Term> TOut;
    STATIC_ASSERT( traits::is_convertible<TOut, T>::value );
   
    GenericOperatorTerm<Term> ot(term);
    std::set<const FiniteElemSpace*> feSpaces = ot.getAuxFeSpaces();
    
    Mesh* mesh = result->getFeSpace()->getMesh();
//     if (feSpaces.size() > 0 && mesh != (*feSpaces.begin())->getMesh())
//       return transformDOF_mm(term, result);
    
    
    DOFVector<TOut> temp(result->getFeSpace(), "temp");
    DOFVector<int> assigned(result->getFeSpace(), "assigned");
    
    const FiniteElemSpace* resultFeSpace = temp.getFeSpace();
    const BasisFunction *basisFcts = resultFeSpace->getBasisFcts();  
    int nBasisFcts = basisFcts->getNumber();
    
    std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1,
					  Mesh::CALL_LEAF_EL | 
					  Mesh::FILL_COORDS | Mesh::FILL_GRD_LAMBDA);
    term.initElement(elInfo, NULL, NULL, basisFcts);
    
    
    TOut tmp(term(0)); nullify(tmp);
    assigned.set(0);
    temp.set(tmp);
    
    while (elInfo) {
      term.initElement(elInfo, NULL, NULL, basisFcts);
      basisFcts->getLocalIndices(elInfo->getElement(), resultFeSpace->getAdmin(), localIndices);	
      
      for (int i = 0; i < nBasisFcts; i++) {
      	temp[localIndices[i]] += term(i);
      	assigned[localIndices[i]]++;
      }
      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(temp);
    Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(assigned);
#endif
    
    DOFIterator<TOut> tempIter(&temp, USED_DOFS);
    DOFIterator<T> resultIter(result, USED_DOFS);
    DOFIterator<int> assignedIter(&assigned, USED_DOFS);
    for (tempIter.reset(), resultIter.reset(), assignedIter.reset(); 
        !resultIter.end(); 
        ++tempIter, ++resultIter, ++assignedIter) {
      *resultIter = (*tempIter);
      *resultIter/= (*assignedIter);
    }
  }

} // end namespace AMDiS
