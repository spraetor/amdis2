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


#include <vector>
#include "CFE_NormAndErrorFcts.h"
#include "Mesh.h"
#include "Traverse.h"
#include "SubElInfo.h"

namespace compositeFEM {

  double CFE_NormAndErrorFcts::L2_err_abs = 0.0;
  double CFE_NormAndErrorFcts::L2_u_norm = 0.0;
  double CFE_NormAndErrorFcts::H1_err_abs = 0.0;
  double CFE_NormAndErrorFcts::H1_u_norm = 0.0;

  double
  ElementL1Norm_Analyt::calcElNorm(ElInfo *elInfo, 
				   const double &det,
				   const double &fac)
  {
    double val = 0.0;
    WorldVector<double> worldCoordsAtQP;

    for (int iq = 0; iq < nQPts; ++iq) {
      elInfo->coordToWorld(q->getLambda(iq), worldCoordsAtQP);
      val += q->getWeight(iq) * fabs((*f)(worldCoordsAtQP));
    }
    double nrm = det * val;

    return nrm;
  }

  double
  ElementL2Norm_Analyt::calcElNorm(ElInfo *elInfo, 
				   const double &det,
				   const double &fac)
  {
    double val = 0.0;
    WorldVector<double> worldCoordsAtQP;

    for (int iq = 0; iq < nQPts; ++iq) {
      elInfo->coordToWorld(q->getLambda(iq), worldCoordsAtQP);
      val += q->getWeight(iq) * sqr((*f)(worldCoordsAtQP));
    }
    double nrm = det * val;

    return nrm;
  }

  double
  ElementH1Norm_Analyt::calcElNorm(ElInfo *elInfo, 
				   const double &det,
				   const double &fac)
  {
    double val = 0.0;
    double norm_grd2 = 0.0;
    WorldVector<double> worldCoordsAtQP;

    for (int iq = 0; iq < nQPts; ++iq) {
      elInfo->coordToWorld(q->getLambda(iq), worldCoordsAtQP);
    
      norm_grd2 = 0.0;
      for (int j = 0; j < dim; j++)
	norm_grd2 += sqr(((*grd)(worldCoordsAtQP))[j]);
    
      val += q->getWeight(iq) * norm_grd2;
    }
    double nrm = det * val;

    return nrm;
  }


  double ElementL1Norm_DOF::calcElNorm(ElInfo *elInfo, 
				       const double &det,
				       const double &fac)
  {
    double val = 0.0;
    mtl::dense_vector<double> dofAtQPs(q->getNumPoints());
    dofVec->getVecAtQPs(elInfo, q, NULL, dofAtQPs); 

    for (int iq = 0; iq < nQPts; ++iq)
      val += q->getWeight(iq) * fabs(dofAtQPs[iq]);
    
    return det * val;
  }


  double ElementL2Norm_DOF::calcElNorm(ElInfo *elInfo, 
				       const double &det,
				       const double &fac)
  {
    double val = 0.0;
    mtl::dense_vector<double> dofAtQPs(q->getNumPoints());
    dofVec->getVecAtQPs(elInfo, q, NULL, dofAtQPs); 

    for (int iq = 0; iq < nQPts; ++iq)
      val += q->getWeight(iq) * sqr(dofAtQPs[iq]);
    
    return det * val;
  }

  double
  ElementH1Norm_DOF::calcElNorm(ElInfo *elInfo, 
				const double &det,
				const double &fac)
  {
    double val = 0.0;
    double norm_grd2;
    mtl::dense_vector<WorldVector<double> > grdDofAtQPs;
    dofVec->getGrdAtQPs(elInfo, q, NULL, grdDofAtQPs);

    for (int iq = 0; iq < nQPts; ++iq) {
    
      norm_grd2 = 0.0;
      for (int j = 0; j < dim; ++j)
	norm_grd2 += sqr(grdDofAtQPs[iq][j]);
    
      val += q->getWeight(iq) * norm_grd2;
    }
    double nrm = det * val;

    return nrm;
  }


  double ElementL2Err::calcElNorm(ElInfo *elInfo, 
				  const double &det, 
				  const double &fac)
  {
    double val = 0.0;
    double val_nrm = 0.0;
    mtl::dense_vector<double> uhAtQPs(q->getNumPoints());
    uh->getVecAtQPs(elInfo, q, NULL, uhAtQPs); 
    WorldVector<double> worldCoordsAtQP;

    for (int iq = 0; iq < nQPts; ++iq) {
      elInfo->coordToWorld(q->getLambda(iq), worldCoordsAtQP);
      val += q->getWeight(iq) * sqr((*u)(worldCoordsAtQP) - uhAtQPs[iq]);
    
      if (relErr) 
	val_nrm += q->getWeight(iq) * sqr((*u)(worldCoordsAtQP));
    }

    double err = det * val;

    if (relErr)
      nrmU += fac * det * val_nrm;

    return err;
  }


  double ElementH1Err::calcElNorm(ElInfo *elInfo, 
				  const double &det,
				  const double &fac)
  {
    double val = 0.0;
    double val_nrm = 0.0;
    double norm_err_grd2;
    double norm_grd2;
    mtl::dense_vector<WorldVector<double> > grdUhAtQPs;
    uh->getGrdAtQPs(elInfo, q, NULL, grdUhAtQPs);
    WorldVector<double> worldCoordsAtQP;

    for (int iq = 0; iq < nQPts; ++iq) {
      elInfo->coordToWorld(q->getLambda(iq), worldCoordsAtQP);

      norm_err_grd2 = 0.0;
      for (int j = 0; j < dim; ++j)
	norm_err_grd2 += 
	  sqr(((*grdu)(worldCoordsAtQP))[j] - grdUhAtQPs[iq][j]);
    
      val += q->getWeight(iq) * norm_err_grd2;
    
      if (relErr) {
	norm_grd2 = 0.0;
	for (int j = 0; j < dim; ++j)
	  norm_grd2 += sqr(((*grdu)(worldCoordsAtQP))[j]);
      
	val_nrm += q->getWeight(iq) * norm_grd2;
      }
    }
    double err = det * val;
  
    if (relErr) 
      nrmGrdU += fac * det * val_nrm;
  
    return err;
  }

  double
  CFE_NormAndErrorFcts::Norm_IntNoBound(ElementNorm *elNorm,
					ElementLevelSet *elLS,
					Flag fillFlag,
					int deg, 
					Quadrature *q)
  {
    int dim = elLS->getDim();
    Mesh *mesh = elLS->getMesh();
    double nrm = 0.0;
    int elStatus;

    // ===== Get quadratures. =====
    if (!q) {
      q = Quadrature::provideQuadrature(dim, deg);
    }
    elNorm->setQuadrature(q);


    // ===== Traverse mesh and calculate integral on each element. =====
    TraverseStack stack;

    ElInfo *elInfo = stack.traverseFirst(mesh, -1, fillFlag);

    while(elInfo) {

      // Check whether current element is cut by the zero level set.
      elStatus = elLS->createElementLevelSet(elInfo);

      if (elStatus == ElementLevelSet::LEVEL_SET_INTERIOR) {

	// -------------------------------------------------------------------
	//  Element is in the domain with negative level set function values.
	// -------------------------------------------------------------------

	elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);

	nrm += elNorm->calcElNorm(elInfo, fabs(elInfo->getDet()));
      }

      elInfo = stack.traverseNext(elInfo);

    }  // end of: mesh traverse

    return nrm;
  }

  double
  CFE_NormAndErrorFcts::Norm_IntBound(ElementNorm *elNorm,
				      ElementLevelSet *elLS,
				      Flag fillFlag,
				      int deg, 
				      Quadrature *q)
  {
    int dim = elLS->getDim();
    Mesh *mesh = elLS->getMesh();
    double nrm = 0.0;
    double el_norm;
    VectorOfFixVecs<DimVec<double> > *intersecPts;
    int numIntersecPts;
    SubPolytope *subPolytope;
    ScalableQuadrature *scalQuad;
//     int nScalQPts;
    int elStatus;

    // ===== Get quadratures. =====
    if (!q) {
      q = Quadrature::provideQuadrature(dim, deg);
    }
    scalQuad = new ScalableQuadrature(q);
//     nScalQPts = scalQuad->getNumPoints();


    // ===== Traverse mesh and calculate integral on each element. =====
    TraverseStack stack;

    ElInfo *elInfo = stack.traverseFirst(mesh, -1, fillFlag);

    while(elInfo) {
      el_norm = 0.0;

      // Check whether current element is cut by the zero level set.
      elStatus = elLS->createElementLevelSet(elInfo);

      if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {

	// -------------------------------------------------------------------
	//  Element is cut by the zero level set.
	// -------------------------------------------------------------------

	// Calculate norm on subpolyope.
	intersecPts = elLS->getElIntersecPoints();
	numIntersecPts = elLS->getNumElIntersecPoints();

	// -----------------------------------------------------------------
	//  Subelement may be inside the domain with negative level set
	//  function value as well as inside the domain with positive
	//  function value.
	//
	//  Whether a subelement is in the domain with negative or positive
	//  level set function values is checked by the level set function
	//  value of the first vertex of the subelement. (The subelements 
	//  are created in such a way that this vertex always is a vertex 
	//  of the element and not an intersection point. Thus the level set 
	//  function value of this vertex really is unequal to zero.)

	subPolytope = new SubPolytope(elInfo, 
				      intersecPts, 
				      numIntersecPts, 
				      0);

	elLS->setLevelSetDomain(
				elLS->getVertexPos(subPolytope->getSubElement(0)->getLambda(0)));

	el_norm = calcSubPolNorm(elInfo,
				 subPolytope,
				 elNorm,
				 scalQuad);

	// Calculate integral on the other element part.
	if (elLS->getLevelSetDomain() == ElementLevelSet::LEVEL_SET_EXTERIOR) 
	  elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);
	else 
	  elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_EXTERIOR);

	elNorm->setQuadrature(q);
	el_norm += elNorm->calcElNorm(elInfo, fabs(elInfo->getDet()));

	el_norm -= calcSubPolNorm(elInfo,
				  subPolytope,
				  elNorm,
				  scalQuad,
				  -1.0);

	nrm += el_norm;
    
	// Free data.
	delete subPolytope;
      }
      else if (elStatus == ElementLevelSet::LEVEL_SET_INTERIOR) {

	// -------------------------------------------------------------------
	//  Element is in the domain with negative level set function values.
	// -------------------------------------------------------------------

	elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);

	elNorm->setQuadrature(q);
	nrm += elNorm->calcElNorm(elInfo, fabs(elInfo->getDet()));
      }

      elInfo = stack.traverseNext(elInfo);

    }  // end of: mesh traverse

    delete scalQuad;

    return nrm;
  }

  double
  CFE_NormAndErrorFcts::Norm_Int(ElementNorm *elNorm,
				 ElementLevelSet *elLS,
				 Flag fillFlag,
				 int deg, 
				 Quadrature *q)
  {
    int dim = elLS->getDim();
    Mesh *mesh = elLS->getMesh();
    double nrm = 0.0;
    double el_norm;
    VectorOfFixVecs<DimVec<double> > *intersecPts;
    int numIntersecPts;
    int vertex_interior;
    SubPolytope *subPolytope;
    ScalableQuadrature *scalQuad;
//     int nScalQPts;
    int elStatus;

    // ===== Get quadratures. =====
    if (!q) {
      q = Quadrature::provideQuadrature(dim, deg);
    }
    scalQuad = new ScalableQuadrature(q);
//     nScalQPts = scalQuad->getNumPoints();


    // ===== Traverse mesh and calculate integral on each element. =====
    TraverseStack stack;

    ElInfo *elInfo = stack.traverseFirst(mesh, -1, fillFlag);

    while(elInfo) {
      el_norm = 0.0;

      // Check whether current element is cut by the zero level set.
      elStatus = elLS->createElementLevelSet(elInfo);

      if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {

	// -------------------------------------------------------------------
	//  Element is cut by the zero level set.
	// -------------------------------------------------------------------

	// Create subelements.
	intersecPts = elLS->getElIntersecPoints();
	numIntersecPts = elLS->getNumElIntersecPoints();

	if (dim == 1 || (dim == 3 && numIntersecPts == 4)) {

	  // -----------------------------------------------------------------
	  //  Subelement(s) are inside the domain with negative level set
	  //  function value.

	  // Get vertex with negative level set function value.
	  for (int i=0; i<=dim; ++i) {
	    if (elLS->getElVertLevelSetVec(i) < 0) {
	      vertex_interior = i;
	      break;
	    }
	  }

	  subPolytope = new SubPolytope(elInfo, 
					intersecPts, 
					numIntersecPts, 
					vertex_interior);

	  elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);
	}
	else {

	  // -----------------------------------------------------------------
	  //  Subelement may be inside the domain with negative level set
	  //  function value as well as inside the domain with positive
	  //  function value.
	  //
	  //  Whether a subelement is in the domain with negative or positive
	  //  level set function values is checked by the level set function
	  //  value of the first vertex of the subelement. (The subelements 
	  //  are created in such a way that this vertex always is a vertex 
	  //  of the element and not an intersection point. Thus the level set 
	  //  function value of this vertex really is unequal to zero.)

	  subPolytope = new SubPolytope(elInfo, 
					intersecPts, 
					numIntersecPts, 
					0);

	  elLS->setLevelSetDomain(
				  elLS->getVertexPos(subPolytope->getSubElement(0)->getLambda(0)));
	}

	// Calculate norm on subpolytope.
	if (elLS->getLevelSetDomain() == ElementLevelSet::LEVEL_SET_INTERIOR)
	  el_norm = calcSubPolNorm(elInfo,
				   subPolytope,
				   elNorm,
				   scalQuad);
	else 
	  el_norm = calcSubPolNorm(elInfo,
				   subPolytope,
				   elNorm,
				   scalQuad,
				   -1.0);

	// -------------------------------------------------------------------
	// In case the subelement is in the domain with positive level set
	// function values:
	// Calculate the integral on the element part with negative
	// level set function values by substracting the integral on the
	// subelement from the integral on the complete element.
	if (elLS->getLevelSetDomain() == ElementLevelSet::LEVEL_SET_EXTERIOR) {

	  elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);

	  elNorm->setQuadrature(q);
	  el_norm *= -1.0;
	  el_norm += elNorm->calcElNorm(elInfo, fabs(elInfo->getDet()));
	}

	// Free data.
	delete subPolytope;
      }
      else if (elStatus == ElementLevelSet::LEVEL_SET_INTERIOR) {

	// -------------------------------------------------------------------
	//  Element is in the domain with negative level set function values.
	// -------------------------------------------------------------------

	elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);

	elNorm->setQuadrature(q);
	el_norm = elNorm->calcElNorm(elInfo, fabs(elInfo->getDet()));
      }

      nrm += el_norm;
      elInfo = stack.traverseNext(elInfo);

    }  // end of: mesh traverse

    delete scalQuad;

    return nrm;
  }

  double
  CFE_NormAndErrorFcts::Norm_Bound(ElementNorm *elNorm,
				   ElementLevelSet *elLS,
				   Flag fillFlag,
				   int deg, 
				   Quadrature *q)
  {
    int dim = elLS->getDim();
    Mesh *mesh = elLS->getMesh();
    double nrm = 0.0;
    double el_norm;
    VectorOfFixVecs<DimVec<double> > *intersecPts;
    int numIntersecPts;
    SubPolytope *subPolytope;
    ScalableQuadrature *scalQuad;
//     int nScalQPts;
    int elStatus;

    // ===== Get quadratures. =====
    if (!q) {
      q = Quadrature::provideQuadrature(dim, deg);
    }
    scalQuad = new ScalableQuadrature(q);
//     nScalQPts = scalQuad->getNumPoints();


    // ===== Traverse mesh and calculate integral on each element. =====
    TraverseStack stack;

    ElInfo *elInfo = stack.traverseFirst(mesh, -1, fillFlag);

    while(elInfo) {
      el_norm = 0.0;

      // Check whether current element is cut by the zero level set.
      elStatus = elLS->createElementLevelSet(elInfo);

      if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {

	// -------------------------------------------------------------------
	//  Element is cut by the zero level set.
	// -------------------------------------------------------------------

	// Calculate norm on subpolyope.
	intersecPts = elLS->getElIntersecPoints();
	numIntersecPts = elLS->getNumElIntersecPoints();

	// -----------------------------------------------------------------
	//  Subelement may be inside the domain with negative level set
	//  function value as well as inside the domain with positive
	//  function value.
	//
	//  Whether a subelement is in the domain with negative or positive
	//  level set function values is checked by the level set function
	//  value of the first vertex of the subelement. (The subelements 
	//  are created in such a way that this vertex always is a vertex 
	//  of the element and not an intersection point. Thus the level set 
	//  function value of this vertex really is unequal to zero.)

	subPolytope = new SubPolytope(elInfo, 
				      intersecPts, 
				      numIntersecPts, 
				      0);

	elLS->setLevelSetDomain(
				elLS->getVertexPos(subPolytope->getSubElement(0)->getLambda(0)));

	el_norm = calcSubPolNorm(elInfo,
				 subPolytope,
				 elNorm,
				 scalQuad);

	// Calculate integral on the other element part.
	if (elLS->getLevelSetDomain() == ElementLevelSet::LEVEL_SET_EXTERIOR) 
	  elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);
	else 
	  elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_EXTERIOR);

	elNorm->setQuadrature(q);
	el_norm += elNorm->calcElNorm(elInfo, fabs(elInfo->getDet()));

	el_norm -= calcSubPolNorm(elInfo,
				  subPolytope,
				  elNorm,
				  scalQuad,
				  -1.0);

	nrm += el_norm;
    
	// Free data.
	delete subPolytope;
      }

      elInfo = stack.traverseNext(elInfo);

    }  // end of: mesh traverse

    delete scalQuad;

    return nrm;
  }

  double
  CFE_NormAndErrorFcts::Norm_Complete(ElementNorm *elNorm,
				      ElementLevelSet *elLS,
				      Flag fillFlag,
				      int deg, 
				      Quadrature *q)
  {
    int dim = elLS->getDim();
    Mesh *mesh = elLS->getMesh();
    double nrm = 0.0;
    double el_norm;
    VectorOfFixVecs<DimVec<double> > *intersecPts;
    int numIntersecPts;
    SubPolytope *subPolytope;
    ScalableQuadrature *scalQuad;
//     int nScalQPts;
    int elStatus;

    // ===== Get quadratures. =====
    if (!q) {
      q = Quadrature::provideQuadrature(dim, deg);
    }
    scalQuad = new ScalableQuadrature(q);
//     nScalQPts = scalQuad->getNumPoints();


    // ===== Traverse mesh and calculate integral on each element. =====
    TraverseStack stack;

    ElInfo *elInfo = stack.traverseFirst(mesh, -1, fillFlag);

    while(elInfo) {
      el_norm = 0.0;

      // Check whether current element is cut by the zero level set.
      elStatus = elLS->createElementLevelSet(elInfo);

      if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {

	// -------------------------------------------------------------------
	//  Element is cut by the zero level set.
	// -------------------------------------------------------------------

	// Calculate norm on subpolyope.
	intersecPts = elLS->getElIntersecPoints();
	numIntersecPts = elLS->getNumElIntersecPoints();

	// -----------------------------------------------------------------
	//  Subelement may be inside the domain with negative level set
	//  function value as well as inside the domain with positive
	//  function value.
	//
	//  Whether a subelement is in the domain with negative or positive
	//  level set function values is checked by the level set function
	//  value of the first vertex of the subelement. (The subelements 
	//  are created in such a way that this vertex always is a vertex 
	//  of the element and not an intersection point. Thus the level set 
	//  function value of this vertex really is unequal to zero.)

	subPolytope = new SubPolytope(elInfo, 
				      intersecPts, 
				      numIntersecPts, 
				      0);

	elLS->setLevelSetDomain(
				elLS->getVertexPos(subPolytope->getSubElement(0)->getLambda(0)));

	el_norm = calcSubPolNorm(elInfo,
				 subPolytope,
				 elNorm,
				 scalQuad);

	// Calculate integral on the other element part.
	if (elLS->getLevelSetDomain() == ElementLevelSet::LEVEL_SET_EXTERIOR) 
	  elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);
	else 
	  elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_EXTERIOR);

	elNorm->setQuadrature(q);
	el_norm += elNorm->calcElNorm(elInfo, fabs(elInfo->getDet()));

	el_norm -= calcSubPolNorm(elInfo,
				  subPolytope,
				  elNorm,
				  scalQuad,
				  -1.0);

	nrm += el_norm;
    
	// Free data.
	delete subPolytope;
      }
      else {

	// -------------------------------------------------------------------
	//  Element is either completely in the domain with negative 
	//  or positive level set function values.
	// -------------------------------------------------------------------

	elLS->setLevelSetDomain(elStatus);

	elNorm->setQuadrature(q);
	nrm += elNorm->calcElNorm(elInfo, fabs(elInfo->getDet()));
      }

      elInfo = stack.traverseNext(elInfo);

    }  // end of: mesh traverse

    delete scalQuad;

    return nrm;
  }

  double 
  CFE_NormAndErrorFcts::L1Norm_Analyt(
				      AbstractFunction<double, WorldVector<double> > *f, 
				      ElementLevelSet *elLS,
				      int domainFlag, 
				      int deg, 
				      Quadrature *q)
  {
    FUNCNAME("CFE_NormAndErrorFcts::L1Norm_Analyt");

    ElementL1Norm_Analyt *elNorm = new ElementL1Norm_Analyt(q, f);
    int dim = elLS->getDim();

    TEST_EXIT(dim == Global::getGeo(WORLD))
      ("doesn't work for dimension of problem != dimension of world!\n");

    Flag fillFlag = Mesh::CALL_LEAF_EL | 
      Mesh::FILL_COORDS |
      Mesh::FILL_DET;

    double nrm = 0.0;
    switch(domainFlag) {
    case -3: nrm = Norm_IntNoBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -2: nrm = Norm_IntBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -1: nrm = Norm_Int(elNorm, elLS, fillFlag, deg, q);
      break;
    case 0: nrm = Norm_Bound(elNorm, elLS, fillFlag, deg, q);
      break;
    case 1: nrm = Norm_Complete(elNorm, elLS, fillFlag, deg, q);
      break;
    default: ERROR_EXIT("illegal flag !\n");
      break;
    }

    delete elNorm;

    return nrm;  
  }

  double 
  CFE_NormAndErrorFcts::L2NormSquare_Analyt(
					    AbstractFunction<double, WorldVector<double> > *f, 
					    ElementLevelSet *elLS,
					    int domainFlag, 
					    int deg, 
					    Quadrature *q)
  {
    FUNCNAME("CFE_NormAndErrorFcts::L2NormSquare_Analyt");

    ElementL2Norm_Analyt *elNorm = new ElementL2Norm_Analyt(q, f);
    int dim = elLS->getDim();

    TEST_EXIT(dim == Global::getGeo(WORLD))
      ("doesn't work for dimension of problem != dimension of world!\n");

    Flag fillFlag = Mesh::CALL_LEAF_EL | 
      Mesh::FILL_COORDS |
      Mesh::FILL_DET;

    double nrm = 0.0;
    switch(domainFlag) {
    case -3: nrm = Norm_IntNoBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -2: nrm = Norm_IntBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -1: nrm = Norm_Int(elNorm, elLS, fillFlag, deg, q);
      break;
    case 0: nrm = Norm_Bound(elNorm, elLS, fillFlag, deg, q);
      break;
    case 1: nrm = Norm_Complete(elNorm, elLS, fillFlag, deg, q);
      break;
    default: ERROR_EXIT("illegal flag !\n");
      break;
    }

    delete elNorm;

    return nrm;  
  }

  double 
  CFE_NormAndErrorFcts::L2Norm_Analyt(
				      AbstractFunction<double, WorldVector<double> > *f, 
				      ElementLevelSet *elLS,
				      int domainFlag, 
				      int deg, 
				      Quadrature *q)
  {
    return sqrt(L2NormSquare_Analyt(f, elLS, domainFlag, deg, q));
  }

  double 
  CFE_NormAndErrorFcts::H1NormSquare_Analyt(
					    AbstractFunction<WorldVector<double>, WorldVector<double> > *grd, 
					    ElementLevelSet *elLS,
					    int domainFlag, 
					    int deg, 
					    Quadrature *q)
  {
    FUNCNAME("CFE_NormAndErrorFcts::H1NormSquare_Analyt");

    int dim = elLS->getDim();
    ElementH1Norm_Analyt *elNorm = new ElementH1Norm_Analyt(q, grd, dim);

    TEST_EXIT(dim == Global::getGeo(WORLD))
      ("doesn't work for dimension of problem != dimension of world!\n");

    Flag fillFlag = Mesh::CALL_LEAF_EL | 
      Mesh::FILL_COORDS |
      Mesh::FILL_DET;

    double nrm = 0.0;
    switch(domainFlag) {
    case -3: nrm = Norm_IntNoBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -2: nrm = Norm_IntBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -1: nrm = Norm_Int(elNorm, elLS, fillFlag, deg, q);
      break;
    case 0: nrm = Norm_Bound(elNorm, elLS, fillFlag, deg, q);
      break;
    case 1: nrm = Norm_Complete(elNorm, elLS, fillFlag, deg, q);
      break;
    default: ERROR_EXIT("illegal flag !\n");
      break;
    }

    delete elNorm;

    return nrm;  
  }

  double 
  CFE_NormAndErrorFcts::H1Norm_Analyt(
				      AbstractFunction<WorldVector<double>, WorldVector<double> > *grd, 
				      ElementLevelSet *elLS,
				      int domainFlag, 
				      int deg, 
				      Quadrature *q)
  {
    return sqrt(H1NormSquare_Analyt(grd, elLS, domainFlag, deg, q));
  }

  double 
  CFE_NormAndErrorFcts::L1Norm_DOF(DOFVector<double> *dof, 
				   ElementLevelSet *elLS,
				   int domainFlag, 
				   int deg, 
				   Quadrature *q)
  {
    FUNCNAME("CFE_NormAndErrorFcts::L1Norm_DOF");

    ElementL1Norm_DOF *elNorm = new ElementL1Norm_DOF(q, dof);
    int dim = elLS->getDim();

    TEST_EXIT(dim == Global::getGeo(WORLD))
      ("doesn't work for dimension of problem != dimension of world!\n");

    Flag fillFlag = Mesh::CALL_LEAF_EL | 
      Mesh::FILL_COORDS |
      Mesh::FILL_DET;

    double nrm = 0.0;
    switch(domainFlag) {
    case -3: nrm = Norm_IntNoBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -2: nrm = Norm_IntBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -1: nrm = Norm_Int(elNorm, elLS, fillFlag, deg, q);
      break;
    case 0: nrm = Norm_Bound(elNorm, elLS, fillFlag, deg, q);
      break;
    case 1: nrm = Norm_Complete(elNorm, elLS, fillFlag, deg, q);
      break;
    default: ERROR_EXIT("illegal flag !\n");
      break;
    }

    delete elNorm;

    return nrm;  
  }

  double 
  CFE_NormAndErrorFcts::L2NormSquare_DOF(DOFVector<double> *dof, 
					 ElementLevelSet *elLS,
					 int domainFlag, 
					 int deg, 
					 Quadrature *q)
  {
    FUNCNAME("CFE_NormAndErrorFcts::L2NormSquare_DOF");

    ElementL2Norm_DOF *elNorm = new ElementL2Norm_DOF(q, dof);
    int dim = elLS->getDim();

    TEST_EXIT(dim == Global::getGeo(WORLD))
      ("doesn't work for dimension of problem != dimension of world!\n");

    Flag fillFlag = Mesh::CALL_LEAF_EL | 
      Mesh::FILL_COORDS |
      Mesh::FILL_DET;

    double nrm = 0.0;
    switch(domainFlag) {
    case -3: nrm = Norm_IntNoBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -2: nrm = Norm_IntBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -1: nrm = Norm_Int(elNorm, elLS, fillFlag, deg, q);
      break;
    case 0: nrm = Norm_Bound(elNorm, elLS, fillFlag, deg, q);
      break;
    case 1: nrm = Norm_Complete(elNorm, elLS, fillFlag, deg, q);
      break;
    default: ERROR_EXIT("illegal flag !\n");
      break;
    }

    delete elNorm;

    return nrm;  
  }

  double 
  CFE_NormAndErrorFcts::L2Norm_DOF(DOFVector<double> *dof, 
				   ElementLevelSet *elLS,
				   int domainFlag, 
				   int deg, 
				   Quadrature *q)
  {
    return sqrt(L2NormSquare_DOF(dof, elLS, domainFlag, deg, q));
  }

  double 
  CFE_NormAndErrorFcts::H1NormSquare_DOF(DOFVector<double> *dof, 
					 ElementLevelSet *elLS,
					 int domainFlag, 
					 int deg, 
					 Quadrature *q)
  {
    FUNCNAME("CFE_NormAndErrorFcts::H1NormSquare_DOF");

    int dim = elLS->getDim();
    ElementH1Norm_DOF *elNorm = new ElementH1Norm_DOF(q, dof, dim);

    TEST_EXIT(dim == Global::getGeo(WORLD))
      ("doesn't work for dimension of problem != dimension of world!\n");

    Flag fillFlag = Mesh::CALL_LEAF_EL | 
      Mesh::FILL_COORDS |
      Mesh::FILL_DET | 
      Mesh::FILL_GRD_LAMBDA;

    double nrm = 0.0;
    switch(domainFlag) {
    case -3: nrm = Norm_IntNoBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -2: nrm = Norm_IntBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -1: nrm = Norm_Int(elNorm, elLS, fillFlag, deg, q);
      break;
    case 0: nrm = Norm_Bound(elNorm, elLS, fillFlag, deg, q);
      break;
    case 1: nrm = Norm_Complete(elNorm, elLS, fillFlag, deg, q);
      break;
    default: ERROR_EXIT("illegal flag !\n");
      break;
    }

    delete elNorm;

    return nrm;  
  }

  double 
  CFE_NormAndErrorFcts::H1Norm_DOF(DOFVector<double> *dof, 
				   ElementLevelSet *elLS,
				   int domainFlag, 
				   int deg, 
				   Quadrature *q)
  {
    return sqrt(H1NormSquare_DOF(dof, elLS, domainFlag, deg, q));
  }

  double 
  CFE_NormAndErrorFcts::L2Err(
			      AbstractFunction<double, WorldVector<double> > *u,
			      DOFVector<double> *uh,
			      ElementLevelSet *elLS,
			      int domainFlag,
			      int relErr,
			      int deg,
			      Quadrature *q)
  {
    FUNCNAME("CFE_NormAndErrorFcts::L2Err()");

    ElementL2Err *elNorm = new ElementL2Err(q, u, uh, relErr);
    int dim = elLS->getDim();

    TEST_EXIT(dim == Global::getGeo(WORLD))
      ("doesn't work for dimension of problem != dimension of world!\n");

    Flag fillFlag = Mesh::CALL_LEAF_EL | 
      Mesh::FILL_COORDS |
      Mesh::FILL_DET;

    double err = 0.0;
    switch(domainFlag) {
    case -3: err = Norm_IntNoBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -2: err = Norm_IntBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -1: err = Norm_Int(elNorm, elLS, fillFlag, deg, q);
      break;
    case 0: err = Norm_Bound(elNorm, elLS, fillFlag, deg, q);
      break;
    case 1: err = Norm_Complete(elNorm, elLS, fillFlag, deg, q);
      break;
    default: ERROR_EXIT("illegal flag !\n");
      break;
    }

    L2_err_abs = sqrt(err);
    L2_u_norm = sqrt(elNorm->getNormU());

    if (relErr)
      err = L2_err_abs / (L2_u_norm + 1.e-15);
    else 
      err = L2_err_abs;

    delete elNorm;

    return err;  
  }

  double 
  CFE_NormAndErrorFcts::H1Err(
			      AbstractFunction<WorldVector<double>, WorldVector<double> > *u,
			      DOFVector<double> *uh,
			      ElementLevelSet *elLS,
			      int domainFlag,
			      int relErr,
			      int deg,
			      Quadrature *q)
  {
    FUNCNAME("CFE_NormAndErrorFcts::H1Err()");

    int dim = elLS->getDim();
    ElementH1Err *elNorm = new ElementH1Err(q, u, uh, relErr, dim);

    TEST_EXIT(dim == Global::getGeo(WORLD))
      ("doesn't work for dimension of problem != dimension of world!\n");

    Flag fillFlag = Mesh::CALL_LEAF_EL | 
      Mesh::FILL_COORDS |
      Mesh::FILL_DET | 
      Mesh::FILL_GRD_LAMBDA;

    double err = 0.0;
    switch(domainFlag) {
    case -3: err = Norm_IntNoBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -2: err = Norm_IntBound(elNorm, elLS, fillFlag, deg, q);
      break;
    case -1: err = Norm_Int(elNorm, elLS, fillFlag, deg, q);
      break;
    case 0: err = Norm_Bound(elNorm, elLS, fillFlag, deg, q);
      break;
    case 1: err = Norm_Complete(elNorm, elLS, fillFlag, deg, q);
      break;
    default: ERROR_EXIT("illegal flag !\n");
      break;
    }

    H1_err_abs = sqrt(err);
    H1_u_norm = sqrt(elNorm->getNormGrdU());

    if (relErr)
      err = H1_err_abs / (H1_u_norm + 1.e-15);
    else 
      err = H1_err_abs;

    delete elNorm;

    return err;  
  }

  double
  CFE_NormAndErrorFcts::calcSubPolNorm(ElInfo *elInfo,
				       SubPolytope *subPolytope,
				       ElementNorm *elNorm,
				       ScalableQuadrature *scalQuad,
				       const double &subPolFac)
  {
    double nrm = 0.0;

    for (std::vector<SubElInfo *>::iterator it = 
	   subPolytope->getSubElementsBegin(); 
	 it != subPolytope->getSubElementsEnd();
	 it++) {

      scalQuad->scaleQuadrature(**it);
      elNorm->setQuadrature(scalQuad);
    
      nrm += elNorm->calcElNorm(elInfo, 
				fabs((*it)->getDet()),
				subPolFac);
    }

    return nrm;
  }

}
