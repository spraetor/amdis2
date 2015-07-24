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


#include "est/SimpleResidualEstimator.h"
#include "Operator.h"
#include "DOFMatrix.h"
#include "DOFVector.h"
#include "Assembler.h"
#include "Traverse.h"
#include "Initfile.h"

namespace AMDiS {

  SimpleResidualEstimator::SimpleResidualEstimator(std::string name) 
    : Estimator(name, 0),
      C0(0.0), 
      C1(0.0)
  {
    // === Read parameters C0 and C1 from init file. ===

    Parameters::get(name + "->C0", C0);
    Parameters::get(name + "->C1", C1);

    C0 = C0 > 1.e-25 ? sqr(C0) : 0.0;
    C1 = C1 > 1.e-25 ? sqr(C1) : 0.0;
  }


  void SimpleResidualEstimator::init(double)
  {
    double kappa = 0.0;
    Parameters::get("kappa", kappa);
    kappa_inv = 1.0 / kappa;

    // === Create data structures. ===

    basFcts = uh[0]->getFeSpace()->getBasisFcts();
    dim = mesh->getDim();
    degree = basFcts->getDegree() * 2;
    quad = Quadrature::provideQuadrature(dim, degree);
    nPoints = quad->getNumPoints();

    Flag flag = INIT_PHI;
    if (degree > 2)
      flag |= INIT_D2_PHI;    
    quadFast = FastQuadrature::provideFastQuadrature(basFcts, *quad, flag);
  
    uhEl.change_dim(basFcts->getNumber()); 
    uhNeigh.change_dim(basFcts->getNumber());
    riq.change_dim(nPoints);

    // === Clear error indicators and mark elements for jumpRes. ===

    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
    while (elInfo) {
      elInfo->getElement()->setEstimation(0.0, 0);
      elInfo->getElement()->setMark(1);
      elInfo = stack.traverseNext(elInfo);
    }

    est_sum = 0.0;
    est_max = 0.0;

    traverseFlag = 
      Mesh::FILL_NEIGH      |
      Mesh::FILL_COORDS     |
      Mesh::FILL_OPP_COORDS |
      Mesh::FILL_BOUND      |
      Mesh::FILL_GRD_LAMBDA |
      Mesh::FILL_DET        |
      Mesh::CALL_LEAF_EL;
    neighInfo = mesh->createNewElInfo();


    // === Prepare date for computing jump residual. ===

    if (C1 > 0.0 && dim > 1) {
      surfaceQuad = Quadrature::provideQuadrature(dim - 1, degree);
      nPointsSurface = surfaceQuad->getNumPoints();
      grdUhEl.resize(nPointsSurface);
      grdUhNeigh.resize(nPointsSurface);
      jump.resize(nPointsSurface);
      localJump.resize(nPointsSurface);
      nNeighbours = Global::getGeo(NEIGH, dim);
      lambdaNeigh = new DimVec<WorldVector<double> >(dim, NO_INIT);
      lambda = new DimVec<double>(dim, NO_INIT);
    }
  }


  void SimpleResidualEstimator::exit(bool output)
  {
    FUNCNAME("SimpleResidualEstimator::exit()");

    // === Calculate the root of the estimations and make output. ===

    est_sum = std::sqrt(est_sum);

    if (output)
      MSG("estimate = %.8e\n", est_sum);

    /// === Delete data structures. ===

    if (C1 && (dim > 1)) {
      delete lambdaNeigh;
      delete lambda;
    }

    delete neighInfo;
  }


  void SimpleResidualEstimator::estimateElement(ElInfo *elInfo, DualElInfo *)
  {    
    FUNCNAME("SimpleResidualEstimator::estimateElement()");

    TEST_EXIT(matrix[0])("Should not happen!\n");

    // Get pointer to element object.
    Element *el = elInfo->getElement();
    // Get error estimation of this element.
    double est_el = el->getEstimation(0);
    // Shortcut for the system matrix.
    DOFMatrix *dofMat = const_cast<DOFMatrix*>(matrix[0]);
    // Shortcut for the right hand side DOF vector
    DOFVector<double> *dofVec = const_cast<DOFVector<double>*>(fh[0]);
    
    // === Init assembler ===

    std::vector<Operator*>::iterator it;
    std::vector<double*>::iterator itfac;
    // Matrix assembler are only initialized with the corresponding multiply
    // factors are != 0
    for (it = dofMat->getOperatorsBegin(), itfac = dofMat->getOperatorEstFactorBegin();
	 it != dofMat->getOperatorsEnd(); ++it, ++itfac)
      if (*itfac == NULL || **itfac != 0.0)
	(*it)->getAssembler()->initElement(elInfo, NULL, quad);	  


    // Vector assembler are only initialized if C0 is set. Note that the jump 
    // residual (thus C1) does not contain the right hand side.
    if (C0 > 0.0)
      for (it = dofVec->getOperatorsBegin(); it != dofVec->getOperatorsEnd(); ++it)
	(*it)->getAssembler()->initElement(elInfo, NULL, quad);	  


    // === Compute element residuals and time error estimation. ===
    if (C0 > 0.0)
      est_el += computeElementResidual(elInfo);

    // === Compute jump residuals. ===
    if (C1 > 0.0 && dim > 1)
      est_el += computeJumpResidual(elInfo);

    // === Update global residual variables. ===
    el->setEstimation(est_el, 0);
    el->setMark(0);
    est_sum += est_el;
    est_max = std::max(est_max, est_el);
  }


  double SimpleResidualEstimator::computeElementResidual(ElInfo *elInfo)
  {
    double det = elInfo->getDet();
    double h2 = h2_from_det(det, dim);
    riq = 0.0;
    DOFMatrix *dofMat = const_cast<DOFMatrix*>(matrix[0]);
    DOFVector<double> *dofVec = const_cast<DOFVector<double>*>(fh[0]);

    // === If there is a valid left hand side operator get the solution  ===
    // === vector or its derivations on the quadrature points of the     ===
    // === element.                                                      ===
    std::vector<Operator*>::iterator it;
    std::vector<double*>::iterator itfac;      
    for (it = dofMat->getOperatorsBegin(), itfac = dofMat->getOperatorEstFactorBegin();
	 it != dofMat->getOperatorsEnd();  ++it, ++itfac) {
      if (*itfac == NULL || **itfac != 0.0) {
	if ((*it)->zeroOrderTerms()) {
	  uhQP.change_dim(nPoints);
	  uh[0]->getVecAtQPs(elInfo, NULL, quadFast, uhQP);
	}

	if (degree > 2 && (*it)->secondOrderTerms()) { 
	  D2UhQp.change_dim(nPoints);
	  uh[0]->getD2AtQPs(elInfo, NULL, quadFast, D2UhQp);
	}
      }
    }
    
    // === Compute the element residual and store it in irq. ===    
    r(elInfo, nPoints, uhQP, D2UhQp, dofMat, dofVec, quad, riq);
    
    
    // === Add integral over r square. ===
    double result = 0.0;
    for (int iq = 0; iq < nPoints; iq++)
      result += quad->getWeight(iq) * riq[iq] * riq[iq];

    double alpha0 = std::min(h2, kappa_inv);
   
    if (norm == NO_NORM || norm == L2_NORM)
      result = C0 * alpha0 * alpha0 * det * result;
    else
      result = C0 * alpha0 * det * result;

    return result;
  }


  double SimpleResidualEstimator::computeJumpResidual(ElInfo *elInfo)
  {
    // === Init temporary variables. ===
    double result = 0.0;
    Element *el = elInfo->getElement();
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
    double det = elInfo->getDet();

    /// === Compute jump on all faces of the current element. ===
    for (int face = 0; face < nNeighbours; face++) {
      // Pointer to neighbouring element.
      Element *neigh = const_cast<Element*>(elInfo->getNeighbour(face));

      // If there is no neighbour, or we have already visited this pair, 
      // then continue with next one.
      if (!(neigh && neigh->getMark()))
	continue;

      // === The next lines are only to set all information about the ===
      // === neighbouring element.                                    ===
      int oppV = elInfo->getOppVertex(face);	
      el->sortFaceIndices(face, faceIndEl);
      neigh->sortFaceIndices(oppV, faceIndNeigh);	
      neighInfo->setElement(const_cast<Element*>(neigh));
      neighInfo->setFillFlag(Mesh::FILL_COORDS);     
      neighInfo->getCoord(oppV) = elInfo->getOppCoord(face);
      
      for (int i = 0; i < dim; i++) {
	int i1 = faceIndEl[i];
	int i2 = faceIndNeigh[i];
	neighInfo->getCoord(i2) = elInfo->getCoord(i1);
      }      
	
      // Compute determinant of the neighbouring element.
      double detNeigh = abs(neighInfo->calcGrdLambda(*lambdaNeigh));
           
      // Reset data.
      for (int iq = 0; iq < nPointsSurface; iq++)
	jump[iq].set(0.0);     
      
      // Get solution vector on the nodes of the current element and 
      // its neighbour.
      uh[0]->getLocalVector(el, uhEl);
      uh[0]->getLocalVector(neigh, uhNeigh);

      // Compute the jump of the gradients of the solution over the face.
      for (int iq = 0; iq < nPointsSurface; iq++) {
	(*lambda)[face] = 0.0;
	for (int i = 0; i < dim; i++)
	  (*lambda)[faceIndEl[i]] = surfaceQuad->getLambda(iq, i);
	
	basFcts->evalGrdUh(*lambda, grdLambda, uhEl, grdUhEl[iq]);
	
	(*lambda)[oppV] = 0.0;
	for (int i = 0; i < dim; i++)
	  (*lambda)[faceIndNeigh[i]] = surfaceQuad->getLambda(iq, i);
	
	basFcts->evalGrdUh(*lambda, *lambdaNeigh, uhNeigh, grdUhNeigh[iq]);
	
	grdUhEl[iq] -= grdUhNeigh[iq];
      } // for iq				
      

      // === Compute the action of the second order operators on the jump. ===
      std::vector<double*>::iterator fac;
      std::vector<Operator*>::iterator it;
      DOFMatrix *mat = const_cast<DOFMatrix*>(matrix[0]);
      for (it = mat->getOperatorsBegin(), fac = mat->getOperatorEstFactorBegin(); 
	   it != mat->getOperatorsEnd(); ++it, ++fac) {
	
	if (*fac == NULL || **fac != 0.0) {
	  for (int iq = 0; iq < nPointsSurface; iq++)
	    localJump[iq].set(0.0);
	  
	  (*it)->weakEvalSecondOrder(grdUhEl, localJump);
	  
	  double factor = *fac ? **fac : 1.0;
	  if (factor != 1.0)
	    for (int i = 0; i < nPointsSurface; i++)
	      localJump[i] *= factor;
	  
	  for (int i = 0; i < nPointsSurface; i++)
	    jump[i] += localJump[i];
	}		
      } // for (it = ...
      

      // === Compute the squared integral of the jump and multiply with ===
      // === the appropriate constants.                                 ===
      
      double val = 0.0;
      for (int iq = 0; iq < nPointsSurface; iq++)
	val += surfaceQuad->getWeight(iq) * (jump[iq] * jump[iq]);

      double d = 0.5 * (det + detNeigh);
      double h2 = h2_from_det(d, dim);

      double alpha0 = std::min(h2, kappa_inv);
   
      if (norm == NO_NORM || norm == L2_NORM)
	val *= C1 * alpha0 * d;
      else
	val *= C1 * d;

      neigh->setEstimation(neigh->getEstimation(0) + val, 0);
      result += val;
    } // for face
    
    /*
    double val = 
      fh[0]->getBoundaryManager()->boundResidual(elInfo, matrix[0], uh[0]);
    if (norm == NO_NORM || norm == L2_NORM)
      val *= C1 * h2;
    else
      val *= C1;   
    result += val;
    */


    return result;
  }


  void SimpleResidualEstimator::r(const ElInfo *elInfo,
				  int nPoints,
				  const mtl::dense_vector<double>& uhIq,
				  const mtl::dense_vector<WorldMatrix<double> > &D2UhIq,
				  DOFMatrix *A, 
				  DOFVector<double> *fh,
				  Quadrature *quad,
				  mtl::dense_vector<double>& result)
  {
    std::vector<Operator*>::iterator it;
    std::vector<double*>::iterator fac;

    // === Get left hand side on quadrature points. ===
    for (it = A->getOperatorsBegin(), fac = A->getOperatorEstFactorBegin(); 
	 it != A->getOperatorsEnd(); ++it, ++fac) {
     
      double factor = *fac ? **fac : 1.0;

      if (factor) {
	if (num_rows(D2UhIq) > 0)
	  (*it)->evalSecondOrder(nPoints, uhIq, grdUhQp, D2UhIq, result, -factor);

	if (num_rows(uhIq) > 0)
	  (*it)->evalZeroOrder(nPoints, uhIq, grdUhQp, D2UhIq, result, factor);
      }
    }
    

    // === Get right hand side on quadrature points. ===
    for (it = fh->getOperatorsBegin(), fac = fh->getOperatorEstFactorBegin(); 
	 it != fh->getOperatorsEnd(); ++it, ++fac) {

      double factor = *fac ? **fac : 1.0;

      if (factor) {
	ElementVector fx(nPoints, 0.0);
	(*it)->getC(elInfo, nPoints, fx);
	
	for (int iq = 0; iq < nPoints; iq++)
	  result[iq] -= factor * fx[iq];
      }
    }    
  }


}
