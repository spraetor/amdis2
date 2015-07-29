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


#include "ResidualEstimator.h"
#include "Operator.h"
#include "DOFMatrix.h"
#include "DOFVector.h"
#include "Assembler.h"
#include "Traverse.h"
#include "Initfile.h"
#include "Parametric.h"
#include "Math.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#include "parallel/MeshDistributor.h"
#include "parallel/ParallelDebug.h"
#endif

using namespace std;

namespace AMDiS {

  ResidualEstimator::ResidualEstimator(std::string name, int r) 
    : Estimator(name, r),
      C0(0.0), 
      C1(0.0), 
      C2(0.0), 
      C3(0.0),
      jumpResidualOnly(false)
  {
    FUNCNAME("ResidualEstimator::ResidualEstimator()");
    
    using math::sqr;

    Parameters::get(name + "->C0", C0);
    Parameters::get(name + "->C1", C1);
    Parameters::get(name + "->C2", C2);
    Parameters::get(name + "->C3", C3);

    C0 = C0 > 1.e-25 ? sqr(C0) : 0.0;
    C1 = C1 > 1.e-25 ? sqr(C1) : 0.0;
    C2 = C2 > 1.e-25 ? sqr(C2) : 0.0;
    C3 = C3 > 1.e-25 ? sqr(C3) : 0.0;

    if (C1 != 0.0 && C0 == 0.0 && C3 == 0.0)
      jumpResidualOnly = true;
      
    TEST_EXIT(C2 == 0.0)("C2 is not used! Please remove it or set it to 0.0!\n");
  }


  void ResidualEstimator::init(double ts)
  {
    FUNCNAME_DBG("ResidualEstimator::init()");

    timestep = ts;
    nSystems = static_cast<int>(uh.size());

    TEST_EXIT_DBG(nSystems > 0)("no system set\n");

    dim = mesh->getDim();
    basFcts = new const BasisFunction*[nSystems];
    quadFast = new FastQuadrature*[nSystems];

    degree = 0;
    for (int system = 0; system < nSystems; system++) {
      basFcts[system] = uh[system]->getFeSpace()->getBasisFcts();
      degree = std::max(degree, basFcts[system]->getDegree());
    }
    degree *= 2;

    quad = Quadrature::provideQuadrature(dim, degree);
    nPoints = quad->getNumPoints();

    Flag flag = INIT_PHI | INIT_GRD_PHI;
    if (degree > 2)
      flag |= INIT_D2_PHI;    

    for (int system = 0; system < nSystems; system++)
      quadFast[system] = FastQuadrature::provideFastQuadrature(basFcts[system], 
							       *quad, 
							       flag);    
  
    uhEl.resize(nSystems);
    uhNeigh.resize(nSystems);
    if (timestep)
      uhOldEl.resize(nSystems);

    for (int system = 0; system < nSystems; system++) {
      uhEl[system].change_dim(basFcts[system]->getNumber()); 
      uhNeigh[system].change_dim(basFcts[system]->getNumber());
      if (timestep)
        uhOldEl[system].change_dim(basFcts[system]->getNumber());
    }

    if (timestep) {
      uhQP.change_dim(nPoints);
      uhOldQP.change_dim(nPoints);
    }

    riq.change_dim(nPoints);

    // clear error indicators and mark elements for jumpRes
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
    while (elInfo) {
      // SIMON: DEL LINE BELOW
      elInfo->getElement()->setEstimation(0.0, row);

      elInfo->getElement()->setMark(1);
      elInfo = stack.traverseNext(elInfo);
    }

    est_sum = 0.0;
    est_max = 0.0;
    est_t_sum = 0.0;
    est_t_max = 0.0;

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

      secondOrderTerms.resize(nSystems);
      for (int system = 0; system < nSystems; system++) {
	secondOrderTerms[system] = false;

	if (matrix[system] == NULL)
	  continue;

	for (std::vector<Operator*>::iterator it = matrix[system]->getOperators().begin();
	     it != matrix[system]->getOperators().end(); ++it)
	  secondOrderTerms[system] = secondOrderTerms[system] || (*it)->secondOrderTerms();
      }
    }
    
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    initParallel();
#endif
  }


#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
  void ResidualEstimator::initParallel()
  {
    FUNCNAME_DBG("ResidualEstimator::initParallel()");

    if (C1 == 0.0)
      return;

    size_t comp = 0;
    for(; comp < uh.size(); comp++)
      if(uh[comp]->getFeSpace()->getMesh() == mesh)
	break;
      
    TEST_EXIT_DBG(comp < uh.size())("No mesh in the solution is found?.\n");
    TEST_EXIT_DBG(*(uh[comp]->getFeSpace()->getAdmin()) == mesh->getDofAdmin(0))("The admin is not the first one in the mesh.\n");
    
    DOFVector<WorldVector<double> > coords(uh[comp]->getFeSpace(), "tmp");
    mesh->getDofIndexCoords(coords);

    Parallel::InteriorBoundary &intBoundary = 
      Parallel::MeshDistributor::globalMeshDistributor->getIntBoundary(0);
      
#if (DEBUG != 0)
    // Make sure interior boundary is correct
    Parallel::ParallelDebug::testInteriorBoundary(*(Parallel::MeshDistributor::globalMeshDistributor));
#endif

    ElInfo *elInfo = mesh->createNewElInfo();
    elInfo->setFillFlag(Mesh::FILL_COORDS);

    Parallel::StdMpi<vector<double> > stdMpiDet(Parallel::MeshDistributor::globalMeshDistributor->getMpiComm(0));
    Parallel::StdMpi<vector<vector<WorldVector<double> > > > stdMpiGrdUh(Parallel::MeshDistributor::globalMeshDistributor->getMpiComm(0));

    Parallel::RankToBoundMap allBounds = intBoundary.getOther();
    allBounds.insert(intBoundary.getOwn().begin(), intBoundary.getOwn().end());

    for (Parallel::RankToBoundMap::iterator it = allBounds.begin();
	 it != allBounds.end(); ++it) {

      vector<BoundaryObject> subBound;

      for (unsigned int i = 0; i < it->second.size(); i++) {
	BoundaryObject &bObj = it->second[i].rankObj;
	if (bObj.subObj == VERTEX)
	  continue;

	bObj.el = intBoundary.getElementPtr(bObj.elIndex, mesh);
	TEST_EXIT_DBG(bObj.el)("No element found in the map.\n");
	bObj.el->getSubBoundary(bObj, subBound);
      }
      
      if (subBound.size() == 0)
	continue;

      WorldVector<int> faceIndices;
      
      for (unsigned int i = 0; i < subBound.size(); i++) {
	Element *el = subBound[i].el;	  
	int oppV = subBound[i].ithObj;
	
	elInfo->setElement(el);
	el->sortFaceIndices(oppV, faceIndices);	
	for (int k = 0; k <= dim; k++)
	  elInfo->getCoord(k) = coords[el->getDof(k, 0)];
	
	double detNeigh = abs(elInfo->calcGrdLambda(*lambdaNeigh));
	stdMpiDet.getSendData(it->first).push_back(detNeigh);
	
	
	for (int system = 0; system < nSystems; system++) {
	  if (matrix[system] == NULL || secondOrderTerms[system] == false)
	    continue;
	  
	  uh[system]->getLocalVector(el, uhNeigh[system]);
	  
	  for (int iq = 0; iq < nPointsSurface; iq++) {
	    
	    (*lambda)[oppV] = 0.0;
	    for (int k = 0; k < dim; k++)
	      (*lambda)[faceIndices[k]] = surfaceQuad->getLambda(iq, k);
	    
	    basFcts[system]->evalGrdUh(*lambda, *lambdaNeigh, uhNeigh[system], grdUhNeigh[iq]);
	  }
	  
	  stdMpiGrdUh.getSendData(it->first).push_back(grdUhNeigh);	  
	}
      }

      stdMpiDet.recv(it->first);
      stdMpiGrdUh.recv(it->first);
    }

    stdMpiDet.updateSendDataSize();
    stdMpiGrdUh.updateSendDataSize();

    stdMpiDet.startCommunication();
    stdMpiGrdUh.startCommunication();
    
    for (Parallel::RankToBoundMap::iterator it = allBounds.begin();
	 it != allBounds.end(); ++it) {
      vector<BoundaryObject> subBound;
	
      for (unsigned int i = 0; i < it->second.size(); i++) {
	BoundaryObject &bObj = it->second[i].rankObj;
	if (bObj.subObj == VERTEX)
	  continue;

	bObj.el = intBoundary.getElementPtr(bObj.elIndex, mesh);
	TEST_EXIT(bObj.el)("No element found in the map.\n");
	bObj.el->getSubBoundary(bObj, subBound);
      }
      
      if (subBound.size() == 0)
	continue;
    
      // Highly possible mesh is not correct, for example: hanging node, if the following error occurs
      // Highly possible mesh is not correct, for example: hanging node
      TEST_EXIT_DBG(subBound.size() == stdMpiDet.getRecvData(it->first).size())
	("Should not happen: rank %d from rank %d should recv %d (recv %d)\n",
	   MPI::COMM_WORLD.Get_rank(), it->first, subBound.size(), stdMpiDet.getRecvData(it->first).size());
      
      TEST_EXIT_DBG(subBound.size() == stdMpiGrdUh.getRecvData(it->first).size())
	("Should not happen!\n");

      for (unsigned int i = 0; i < subBound.size(); i++) {
	elBoundDet[subBound[i]] = stdMpiDet.getRecvData(it->first)[i];
	elBoundGrdUhNeigh[subBound[i]] = stdMpiGrdUh.getRecvData(it->first)[i];
      }
    }    

    delete elInfo;
  }
#endif


  void ResidualEstimator::exit(bool output)
  {
    FUNCNAME("ResidualEstimator::exit()");

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    double send_est_sum = est_sum;
    double send_est_max = est_max;
    double send_est_t_sum = est_t_sum;
    double send_est_t_max = est_t_max;

    MPI::COMM_WORLD.Allreduce(&send_est_sum, &est_sum, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&send_est_max, &est_max, 1, MPI_DOUBLE, MPI_MAX);
    MPI::COMM_WORLD.Allreduce(&send_est_t_sum, &est_t_sum, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&send_est_t_max, &est_t_max, 1, MPI_DOUBLE, MPI_MAX);
#endif

    est_sum = std::sqrt(est_sum);
    est_t_sum = std::sqrt(est_t_sum);

    if (output) {
      MSG("estimate for component %d = %.8e\n", row, est_sum);
      if (C3) {
        MSG("time estimate for component %d = %.8e\n", row, est_t_sum);
      }
    }

    delete [] basFcts;
    delete [] quadFast;

    if (C1 && (dim > 1)) {
      delete lambdaNeigh;
      delete lambda;
    }

    delete neighInfo;
  }


  void ResidualEstimator::estimateElement(ElInfo *elInfo, DualElInfo *dualElInfo)
  {    
    FUNCNAME_DBG("ResidualEstimator::estimateElement()");

    TEST_EXIT_DBG(nSystems > 0)("no system set\n");

    Element *el = elInfo->getElement();
    double est_el = el->getEstimation(row);
    // SIMON    double est_el = 0.0;
    std::vector<Operator*>::iterator it;
    std::vector<double*>::iterator itfac;

    // === Init assemblers. ===
    for (int system = 0; system < nSystems; system++) {
      if (matrix[system] == NULL) 
	continue;

      DOFMatrix *dofMat = const_cast<DOFMatrix*>(matrix[system]);
      DOFVector<double> *dofVec = const_cast<DOFVector<double>*>(fh[system]);

      for (it = dofMat->getOperatorsBegin(), itfac = dofMat->getOperatorEstFactorBegin();
	   it != dofMat->getOperatorsEnd(); ++it, ++itfac)
	if (*itfac == NULL || **itfac != 0.0) {	  
	  // If the estimator must only compute the jump residual but there are no
	  // second order terms in the operator, it can be skipped.
	  if (jumpResidualOnly && (*it)->secondOrderTerms() == false)
	    continue;
	  
	  if (dualElInfo)
	    (*it)->getAssembler()->initElement(dualElInfo->smallElInfo, 
					       dualElInfo->largeElInfo,
					       quad);
	  else
	    (*it)->getAssembler()->initElement(elInfo, NULL, quad);	  
	}

      if (C0 > 0.0)
	for (it = dofVec->getOperatorsBegin(); it != dofVec->getOperatorsEnd(); ++it) {
	  if (dualElInfo)
	    (*it)->getAssembler()->initElement(dualElInfo->smallElInfo, 
					       dualElInfo->largeElInfo,
					       quad);
	  else
	    (*it)->getAssembler()->initElement(elInfo, NULL, quad);	  
	}
    }


    // === Compute element residuals and time error estimation. ===
    if (C0 || C3)
      est_el += computeElementResidual(elInfo, dualElInfo);

    // === Compute jump residuals. ===
    if (C1 && dim > 1)
      est_el += computeJumpResidual(elInfo, dualElInfo);

    // === Update global residual variables. ===
    el->setEstimation(est_el, row);
    el->setMark(0);
    est_sum += est_el;
    est_max = std::max(est_max, est_el);
  }


  double ResidualEstimator::computeElementResidual(ElInfo *elInfo, 
						   DualElInfo *dualElInfo)
  {
    FUNCNAME("ResidualEstimator::computeElementResidual()");

    TEST_EXIT(!dualElInfo)("Not yet implemented!\n");

    std::vector<Operator*>::iterator it;
    std::vector<double*>::iterator itfac;
    double det = elInfo->getDet();
    double h2 = h2_from_det(det, dim);
    riq = 0.0;

    for (int system = 0; system < nSystems; system++) {
      if (matrix[system] == NULL) 
	continue;

      if (timestep && uhOld[system]) {
	TEST_EXIT_DBG(uhOld[system])("no uhOld\n");
	uhOld[system]->getLocalVector(elInfo->getElement(), uhOldEl[system]);
  
	// === Compute time error. ===

	if (C0 > 0.0 || C3 > 0.0) {   
	  uh[system]->getVecAtQPs(elInfo, NULL, quadFast[system], uhQP);
	  uhOld[system]->getVecAtQPs(elInfo, NULL, quadFast[system], uhOldQP);
	  
	  if (C3 > 0.0 && system == std::max(row, 0)) {
	    double result = 0.0;
	    for (int iq = 0; iq < nPoints; iq++) {
	      double tiq = (uhQP[iq] - uhOldQP[iq]);
	      result += quad->getWeight(iq) * tiq * tiq;
	    }
	    double v = C3 * det * result;
	    est_t_sum += v;
	    est_t_max = std::max(est_t_max, v);
	  }
	}
      }
           
      // === Compute element residual. ===
      if (C0 > 0.0) {
	DOFMatrix *dofMat = const_cast<DOFMatrix*>(matrix[system]);
	DOFVector<double> *dofVec = const_cast<DOFVector<double>*>(fh[system]);
  
	for (it = dofMat->getOperatorsBegin(), itfac = dofMat->getOperatorEstFactorBegin();
	     it != dofMat->getOperatorsEnd();  ++it, ++itfac) {
	  if (*itfac == NULL || **itfac != 0.0) {
	    if ((*it)->zeroOrderTerms()) {
	      uhQP.change_dim(nPoints);
	      uh[system]->getVecAtQPs(elInfo, NULL, quadFast[system], uhQP);
	    }
	    if ((*it)->firstOrderTermsGrdPsi() || (*it)->firstOrderTermsGrdPhi()) {
	      grdUhQp.change_dim(nPoints);
	      uh[system]->getGrdAtQPs(elInfo, NULL, quadFast[system], grdUhQp);
	    }
	    if (degree > 2 && (*it)->secondOrderTerms()) {
	      D2UhQp.change_dim(nPoints);
	      uh[system]->getD2AtQPs(elInfo, NULL, quadFast[system], D2UhQp);
	    }
	  }
	}
	
	// === Compute the element residual and store it in irq. ===

	r(elInfo,
	  nPoints, 
	  uhQP,
	  grdUhQp,
	  D2UhQp,
	  uhOldQP,
	  grdUhOldQp,  // not used
	  D2UhOldQp,  // not used
	  dofMat, 
	  dofVec,
	  quad,
	  riq);
      }     
    }

    // add integral over r square
    double result = 0.0;
    for (int iq = 0; iq < nPoints; iq++)
      result += quad->getWeight(iq) * riq[iq] * riq[iq];
   
    if (timestep != 0.0 || norm == NO_NORM || norm == L2_NORM)
      result = C0 * h2 * h2 * det * result;
    else
      result = C0 * h2 * det * result;
    
    return result;
  }


  double ResidualEstimator::computeJumpResidual(ElInfo *elInfo, 
						DualElInfo *dualElInfo)
  {
    FUNCNAME_DBG("ResidualEstimator::computeJumpResidual()");

    double result = 0.0;
    Element *el = elInfo->getElement();
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
    double det = elInfo->getDet();
    double h2 = h2_from_det(det, dim);
    BoundaryObject blub;

    for (int face = 0; face < nNeighbours; face++) {  
      Element *neigh = const_cast<Element*>(elInfo->getNeighbour(face));

      bool parallelMode = false;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
      if (neigh == NULL) {
	BoundaryObject testObj(el, elInfo->getType(), EDGE, face);
	
	if (elBoundDet.count(testObj)) {
	  parallelMode = true;
	  blub = testObj;
	} else {
	  continue;
	}
      }

      if (!parallelMode && !(neigh && neigh->getMark()))
	continue;
#else
      if (!(neigh && neigh->getMark()))
	continue;
#endif

      int oppV = elInfo->getOppVertex(face);	
      el->sortFaceIndices(face, faceIndEl);

      if (!parallelMode) {
	neigh->sortFaceIndices(oppV, faceIndNeigh);	
	neighInfo->setElement(const_cast<Element*>(neigh));
	neighInfo->setFillFlag(Mesh::FILL_COORDS);     
	neighInfo->getCoord(oppV) = elInfo->getOppCoord(face);
      }
      
      // periodic leaf data ?
      ElementData *ldp = el->getElementData()->getElementData(PERIODIC);	
      bool periodicCoords = false;
      
      if (ldp) {
	typedef std::list<LeafDataPeriodic::PeriodicInfo> PerInfList;
	PerInfList& infoList = dynamic_cast<LeafDataPeriodic*>(ldp)->getInfoList();
	
	for (PerInfList::iterator it = infoList.begin(); it != infoList.end(); ++it) {
	  if (it->elementSide == face) {
	    for (int i = 0; i < dim; i++) {
	      int i1 = faceIndEl[i];
	      int i2 = faceIndNeigh[i];
	      
	      int j = 0;
	      for (; j < dim; j++)
		if (i1 == el->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim), face, j))
		  break;
	      
	      TEST_EXIT_DBG(j != dim)("vertex i1 not on face ???\n");
	      
	      neighInfo->getCoord(i2) = (*(it->periodicCoords))[j];
	    }
	    periodicCoords = true;
	    break;
	  }
	}
      }  // if (ldp)
      
      if (!periodicCoords && !parallelMode) {
	for (int i = 0; i < dim; i++) {
	  int i1 = faceIndEl[i];
	  int i2 = faceIndNeigh[i];
	  neighInfo->getCoord(i2) = elInfo->getCoord(i1);
	}
      }
	
      double detNeigh = 0.0;
      Parametric *parametric = mesh->getParametric();
      if (!parallelMode) {
	if (parametric)
	  neighInfo = parametric->addParametricInfo(neighInfo);	  

	detNeigh = abs(neighInfo->calcGrdLambda(*lambdaNeigh));
      } else {
	TEST_EXIT_DBG(!parametric)("No yet implemented!\n");

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	detNeigh = elBoundDet[blub];
#endif
      }
      
           
      for (int iq = 0; iq < nPointsSurface; iq++)
	jump[iq].set(0.0);     
      
      for (int system = 0; system < nSystems; system++) {
	if (matrix[system] == NULL || secondOrderTerms[system] == false) 
	  continue;
	      
	uh[system]->getLocalVector(el, uhEl[system]);	
	if (!parallelMode)
	  uh[system]->getLocalVector(neigh, uhNeigh[system]);
	  
	for (int iq = 0; iq < nPointsSurface; iq++) {
	  (*lambda)[face] = 0.0;
	  for (int i = 0; i < dim; i++)
	    (*lambda)[faceIndEl[i]] = surfaceQuad->getLambda(iq, i);
	  
	  basFcts[system]->evalGrdUh(*lambda, grdLambda, uhEl[system], grdUhEl[iq]);

	  if (!parallelMode) {
	    (*lambda)[oppV] = 0.0;
	    for (int i = 0; i < dim; i++)
	      (*lambda)[faceIndNeigh[i]] = surfaceQuad->getLambda(iq, i);

	    basFcts[system]->evalGrdUh(*lambda, *lambdaNeigh, uhNeigh[system], grdUhNeigh[iq]);
	  } else {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	    grdUhNeigh[iq] = elBoundGrdUhNeigh[blub][iq];
#endif
	  }
	  

	  
	  grdUhEl[iq] -= grdUhNeigh[iq];
	} // for iq				
	
	std::vector<double*>::iterator fac;
	std::vector<Operator*>::iterator it;
	DOFMatrix *mat = const_cast<DOFMatrix*>(matrix[system]);
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
      } // for system
    
      double val = 0.0;
      for (int iq = 0; iq < nPointsSurface; iq++)
	val += surfaceQuad->getWeight(iq) * (jump[iq] * jump[iq]);

      double d = 0.5 * (det + detNeigh);
   
      if (norm == NO_NORM || norm == L2_NORM)
	val *= C1 * h2_from_det(d, dim) * d;
      else
	val *= C1 * d;

      if (!parallelMode) {
	if (parametric)
	  neighInfo = parametric->removeParametricInfo(neighInfo);      
	neigh->setEstimation(neigh->getEstimation(row) + val, row);
      }

      result += val;
    } // for face
    
    double val = fh[std::max(row, 0)]->getBoundaryManager()->
      boundResidual(elInfo, matrix[std::max(row, 0)], uh[std::max(row, 0)]);    
    if (norm == NO_NORM || norm == L2_NORM)
      val *= C1 * h2;
    else
      val *= C1;   
    result += val;

    return result;
  }


  void r(const ElInfo *elInfo,
	 int nPoints,
	 const DenseVector<double>& uhIq,
	 const DenseVector<WorldVector<double> > &grdUhIq,
	 const DenseVector<WorldMatrix<double> > &D2UhIq,
	 const DenseVector<double>& uhOldIq,
	 const DenseVector<WorldVector<double> > &grdUhOldIq,
	 const DenseVector<WorldMatrix<double> > &D2UhOldIq,
	 DOFMatrix *A, 
	 DOFVector<double> *fh,
	 Quadrature *quad,
	 ElementVector& result)
  {
    std::vector<Operator*>::iterator it;
    std::vector<double*>::iterator fac;

    // lhs
    for (it = A->getOperatorsBegin(), fac = A->getOperatorEstFactorBegin(); 
	 it != A->getOperatorsEnd(); ++it, ++fac) {
     
      double factor = *fac ? **fac : 1.0;

      if (factor) {
	if (num_rows(D2UhIq) > 0)
	  (*it)->evalSecondOrder(nPoints, uhIq, grdUhIq, D2UhIq, result, -factor);	

	if (num_rows(grdUhIq) > 0) {
	  (*it)->evalFirstOrderGrdPsi(nPoints, uhIq, grdUhIq, D2UhIq, result, factor);
	  (*it)->evalFirstOrderGrdPhi(nPoints, uhIq, grdUhIq, D2UhIq, result, factor);
	}
	
	if (num_rows(uhIq) > 0)
	  (*it)->evalZeroOrder(nPoints, uhIq, grdUhIq, D2UhIq, result, factor);	
      }
    }
    
    // rhs
    for (it = fh->getOperatorsBegin(), fac = fh->getOperatorEstFactorBegin(); 
	 it != fh->getOperatorsEnd(); ++it, ++fac) {

      double factor = *fac ? **fac : 1.0;

      if (factor) {
	if ((*it)->getUhOld()) {
	  if (num_rows(D2UhOldIq) > 0)
	    (*it)->evalSecondOrder(nPoints, uhOldIq, grdUhOldIq, D2UhOldIq, result, factor);
	  
	  if (num_rows(grdUhOldIq) > 0) {
	    (*it)->evalFirstOrderGrdPsi(nPoints, uhOldIq, grdUhOldIq, D2UhOldIq, result, -factor);
	    (*it)->evalFirstOrderGrdPhi(nPoints, uhOldIq, grdUhOldIq, D2UhOldIq, result, -factor);
	  }

	  if (num_rows(uhOldIq) > 0)
	    (*it)->evalZeroOrder(nPoints, uhOldIq, grdUhOldIq, D2UhOldIq, result, -factor);	  
	} else {
	  ElementVector fx(nPoints, 0.0);
	  (*it)->getC(elInfo, nPoints, fx);

	  for (int iq = 0; iq < nPoints; iq++)
	    result[iq] -= factor * fx[iq];
	}
      }
    }    
  }


}
