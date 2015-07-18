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



/** \file Expressions.hh */

namespace AMDiS {

  template <class Expr>
  inline typename boost::enable_if<traits::is_expr<Expr>, typename Expr::value_type>::type
  integrate(Expr expr, Mesh* mesh_)
  {
    typedef typename Expr::value_type TOut;
    
    GenericOperatorTerm<Expr> ot(expr);
    std::set<const FiniteElemSpace*> feSpaces = ot.getAuxFeSpaces();
    
    TEST_EXIT(mesh_ || !feSpaces.empty())
      ("The Expression must contain a DOFVector or FeSpace depended value!\n");
    Mesh* mesh = mesh_ ? mesh_ : (*feSpaces.begin())->getMesh();
    
    int deg = expr.getDegree();
    int dim = mesh->getDim();
    Quadrature* quad = Quadrature::provideQuadrature(dim, deg);
    
    TOut value; nullify(value);

    Flag traverseFlag = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_GRD_LAMBDA | Mesh::FILL_DET;
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, traverseFlag);
    while (elInfo) {
      expr.initElement(&ot, elInfo, NULL, quad);
      TOut tmp; nullify(tmp);
      for (int iq = 0; iq < quad->getNumPoints(); iq++) {
	tmp += quad->getWeight(iq) * expr(iq);
      }
      value += tmp * elInfo->getDet();

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(value);
#endif

    return value;
  }


  template <class Expr, class Functor>
  inline typename boost::enable_if<traits::is_expr<Expr>, typename Expr::value_type>::type
  accumulate(Expr expr, Functor f, typename Expr::value_type value0)
  {
    GenericOperatorTerm<Expr> ot(expr);
    std::set<const FiniteElemSpace*> feSpaces = ot.getAuxFeSpaces();
    
    TEST_EXIT(!feSpaces.empty())
      ("The Expression must contain a DOFVector or FeSpace depended value!\n");
    const FiniteElemSpace* feSpace0 = *feSpaces.begin();  
    Mesh* mesh = feSpace0->getMesh();
    
    const BasisFunction *basisFcts = feSpace0->getBasisFcts();  
    int nBasisFcts = basisFcts->getNumber();
    
    DOFVector<bool> assigned(feSpace0, "assigned");
    assigned.set(false);
    
    std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1,
					  Mesh::CALL_LEAF_EL | 
					  Mesh::FILL_COORDS | Mesh::FILL_GRD_LAMBDA);
    
    while (elInfo) {
      expr.initElement(&ot, elInfo, NULL, NULL, basisFcts);
      basisFcts->getLocalIndices(elInfo->getElement(), feSpace0->getAdmin(), localIndices);	
      
      for (int i = 0; i < nBasisFcts; i++) {
	if (!assigned[localIndices[i]]) {
	  value0 = f(value0, expr(i));
	  assigned[localIndices[i]] = true;
	}
      }
      elInfo = stack.traverseNext(elInfo);
    }
    
    return value0;
  }


  // works only for nodal basis functions!
  template <class T, class Expr>
  inline typename boost::enable_if<
    typename and_<traits::is_expr<Expr>, 
		  traits::is_convertible<typename Expr::value_type, T>
		  >::type
    >::type
  transformDOF(Expr expr, DOFVector<T>* result)
  {
    GenericOperatorTerm<Expr> ot(expr);
    std::set<const FiniteElemSpace*> feSpaces = ot.getAuxFeSpaces();
    
    Mesh* mesh = result->getFeSpace()->getMesh();
    if (feSpaces.size() > 0 && mesh != (*feSpaces.begin())->getMesh())
      return transformDOF_mm(expr, result);
    
    
    typedef typename Expr::value_type TOut;
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
    expr.initElement(&ot, elInfo, NULL, NULL, basisFcts);
    
    
    TOut tmp(expr(0)); nullify(tmp);
    assigned.set(0);
    temp.set(tmp);
    
    while (elInfo) {
      expr.initElement(&ot, elInfo, NULL, NULL, basisFcts);
      basisFcts->getLocalIndices(elInfo->getElement(), resultFeSpace->getAdmin(), localIndices);	
      
      for (int i = 0; i < nBasisFcts; i++) {
	temp[localIndices[i]] += expr(i);
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



  // works only for nodal basis functions!
  template <class T, class Expr>
  inline typename boost::enable_if<
    typename and_<traits::is_expr<Expr>, 
		  traits::is_convertible<typename Expr::value_type, T>
		  >::type
    >::type
  transformDOF_mm(Expr expr, DOFVector<T>* result)
  {
    typedef typename Expr::value_type TOut;
    
    GenericOperatorTerm<Expr> ot(expr);
    std::set<const FiniteElemSpace*> feSpaces = ot.getAuxFeSpaces();
    
    Mesh* mesh1 = result->getFeSpace()->getMesh();
    Mesh* mesh2 = (*feSpaces.begin())->getMesh();
    
    DOFVector<TOut> temp(result->getFeSpace(), "temp");
    DOFVector<int> assigned(result->getFeSpace(), "assigned");
    
    const FiniteElemSpace* resultFeSpace = temp.getFeSpace();
    const BasisFunction *basisFcts = resultFeSpace->getBasisFcts();  
    int nBasisFcts = basisFcts->getNumber();
      
    std::vector<DegreeOfFreedom> localIndices(nBasisFcts);
    mtl::dense_vector<TOut> vecLocalCoeffs(nBasisFcts);

    DimVec<double> *lambda = NULL;
    DimVec<double> *lambda_1 = new DimVec<double>;
    WorldVector<double> coords;
	    
    DualTraverse dualTraverse;
    DualElInfo dualElInfo;
    
    Flag assembleFlag = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_GRD_LAMBDA;
    bool cont = dualTraverse.traverseFirst(mesh1, mesh2, -1, -1, 
					    assembleFlag, assembleFlag, dualElInfo);
    expr.initElement(&ot, dualElInfo.colElInfo, NULL, NULL, basisFcts);
    
    TOut tmp(expr(0)); nullify(tmp);
    assigned.set(0);
    temp.set(tmp);
    
    while (cont) {
      expr.initElement(&ot, dualElInfo.colElInfo, NULL, NULL, basisFcts);
      basisFcts->getLocalIndices(dualElInfo.rowElInfo->getElement(), resultFeSpace->getAdmin(), localIndices);
      
      for (int i = 0; i < nBasisFcts; i++)
	vecLocalCoeffs[i] = expr(i);
      
      for (int i = 0; i < nBasisFcts; i++) {
	lambda = basisFcts->getCoords(i);
	dualElInfo.rowElInfo->coordToWorld(*lambda, coords);
	int inside = dualElInfo.colElInfo->worldToCoord(coords, lambda_1);
	if (inside < 0) {
	  temp[localIndices[i]] += basisFcts->evalUh(*lambda_1, vecLocalCoeffs);
	  assigned[localIndices[i]]++;
	}
      }
      cont = dualTraverse.traverseNext(dualElInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(temp);
    Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(assigned);
#endif
    
    
    DOFIterator<TOut> tempIter(&temp, USED_DOFS);
    DOFIterator<T> resultIter(result, USED_DOFS);
    DOFIterator<int> assignedIter(&assigned, USED_DOFS);
    for (tempIter.reset(), resultIter.reset(), assignedIter.reset(); !resultIter.end(); ++tempIter, ++resultIter, ++assignedIter) {
      *resultIter = (*tempIter);
      *resultIter/= (*assignedIter);
    }
    
    delete lambda_1;
  }

} // end namespace AMDiS
