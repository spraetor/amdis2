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


#include "CFE_Integration.h"
#include "Mesh.h"
#include "SurfaceQuadrature.h"
#include "Traverse.h"
#include "ScalableQuadrature.h"
#include "SubElInfo.h"
#include "SubPolytope.h"

namespace compositeFEM {

  double 
  CFE_Integration::integrate_onNegLs(ElementFunction<double> *f, 
				     ElementLevelSet *elLS,
				     int deg, 
				     Quadrature *q)
  {
    int dim = elLS->getDim();
    double int_val = 0.0;
    double el_int_val;
    double subEl_int_val;
    VectorOfFixVecs<DimVec<double> > *intersecPts;
    SubPolytope *subPolytope;
    int numIntersecPts;
    int iq;
    double val;
    int numQuadPts;
    int vertex_interior;
    ScalableQuadrature *loc_scalQuad;
    int numScalQuadPts;
    bool subPolIsExterior = false;
    int elStatus;
    Mesh *mesh = elLS->getMesh();

    // ===== Get quadratures. =====
    if (!q) {
      q = Quadrature::provideQuadrature(dim, deg);
    }
    numQuadPts = q->getNumPoints();
    loc_scalQuad = new ScalableQuadrature(q);
    numScalQuadPts = loc_scalQuad->getNumPoints();


    // ===== Traverse mesh and calculate integral on each element. =====
    TraverseStack stack;

    ElInfo *loc_elInfo = stack.traverseFirst(mesh,
					     -1, 
					     Mesh::CALL_LEAF_EL | 
					     Mesh::FILL_COORDS |
					     Mesh::FILL_DET);
    while(loc_elInfo) {
      el_int_val = 0.0;
      subEl_int_val = 0.0;
      subPolIsExterior = false;

      // Check whether current element is cut by the zero level set.
      elStatus = elLS->createElementLevelSet(loc_elInfo);

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

	  subPolytope = new SubPolytope(loc_elInfo, 
					intersecPts, 
					numIntersecPts, 
					vertex_interior);
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

	  subPolytope = new SubPolytope(loc_elInfo, 
					intersecPts, 
					numIntersecPts, 
					0);

	  if(elLS->getVertexPos(subPolytope->getSubElement(0)->getLambda(0)) ==
	     ElementLevelSet::LEVEL_SET_EXTERIOR) {
	    subPolIsExterior = true;
	  }
	}

	elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);

	// Calculate integral on subelement(s).
	f->setElInfo(loc_elInfo);
	for (std::vector<SubElInfo *>::iterator it = 
	       subPolytope->getSubElementsBegin(); 
	     it != subPolytope->getSubElementsEnd();
	     it++) {

	  loc_scalQuad->scaleQuadrature(**it);

	  for (val = iq = 0; iq < numScalQuadPts; ++iq) {
	    val += loc_scalQuad->getWeight(iq)*(*f)(loc_scalQuad->getLambda(iq));
	  }
	  el_int_val += fabs((*it)->getDet()) * val;
	}

	// -------------------------------------------------------------------
	// In case the subelement is in the domain with positive level set
	// function values:
	// Calculate the integral on the element part with negative
	// level set function values by substracting the integral on the
	// subelement from the integral on the complete element.
	if (subPolIsExterior) {
	  subEl_int_val = el_int_val;
	  el_int_val = 0.0;

	  f->setElInfo(loc_elInfo);
	  for (val = iq = 0; iq < numQuadPts; ++iq) {
	    val += q->getWeight(iq)*(*f)(q->getLambda(iq));
	  }
	  el_int_val =  loc_elInfo->calcDet()*val - subEl_int_val;
	}

	// Free data.
	delete subPolytope;
      }
      else if (elStatus == ElementLevelSet::LEVEL_SET_INTERIOR) {

	// -------------------------------------------------------------------
	//  Element is in the domain with negative level set function values.
	// -------------------------------------------------------------------

	elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);

	f->setElInfo(loc_elInfo);
	for (val = iq = 0; iq < numQuadPts; ++iq) {
	  val += q->getWeight(iq)*(*f)(q->getLambda(iq));
	}
	el_int_val = loc_elInfo->calcDet() * val;
      }

      int_val += el_int_val;
      loc_elInfo = stack.traverseNext(loc_elInfo);

    }  // end of: mesh traverse

    delete loc_scalQuad;

    return int_val;
  }

  double
  CFE_Integration::integrate_onZeroLs(ElementFunction<double> *f, 
				      ElementLevelSet *elLS,
				      int deg, 
				      Quadrature *q)
  {
    int dim = elLS->getDim();
    double int_val = 0.0;
    VectorOfFixVecs<DimVec<double> > *intersecPts;
    int numIntersecPts;
    int iq;
    double val;
    SurfaceQuadrature *surfQuad;
    int numQuadPts;
    VectorOfFixVecs<DimVec<double> > tmpPts(dim, dim, NO_INIT);
    DimVec<double> tmpPt(dim, 0.0);
    int elStatus;
    Mesh *mesh = elLS->getMesh();

    // ===== Define default points for surface quadrature. =====
    for (int i=0; i<dim; ++i) {
      tmpPt[i] = 1.0;
      tmpPts[i] = tmpPt;
      tmpPt[i] = 0.0;
    }

    // ===== Get quadratures. =====
    if (!q) {
      q = Quadrature::provideQuadrature(dim-1, deg);
    }
    surfQuad = new SurfaceQuadrature(q, tmpPts);
    numQuadPts = surfQuad->getNumPoints();


    // ===== Traverse mesh and calculate integral on each element cut by
    //       the zero level set. =====
    TraverseStack stack;

    ElInfo *loc_elInfo = stack.traverseFirst(mesh,
					     -1, 
					     Mesh::CALL_LEAF_EL | 
					     Mesh::FILL_COORDS |
					     Mesh::FILL_DET);
    while(loc_elInfo) {

      // Check whether current element is cut by the zero level set.
      elStatus = elLS->createElementLevelSet(loc_elInfo);

      if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {
 
	// Get intersection points.
	intersecPts = elLS->getElIntersecPoints();
	numIntersecPts = elLS->getNumElIntersecPoints();

	// Calculate surface integral on intersection plane.
	f->setElInfo(loc_elInfo);

	// Note: The vector *intersecPts has always MAX_INTERSECTION_POINTS
	//       entries.
	for (int i=0; i<dim; ++i) 
	  tmpPts[i] = (*intersecPts)[i];

	surfQuad->scaleSurfaceQuadrature(tmpPts);

	for (val = iq = 0; iq < numQuadPts; ++iq) {
	  val += surfQuad->getWeight(iq)*(*f)(surfQuad->getLambda(iq));
	}
	int_val += calcSurfaceDet(loc_elInfo, tmpPts) * val;

	if (dim == 3 && numIntersecPts == 4) {

	  // -----------------------------------------------------------------
	  // Intersection plane must be divided into two simplices. Calculate
	  // the surface integral for the second simplex.
	  //
	  // Note: The intersection points S0, S1, S2, S3 are supposed to be 
	  //       alligned in such a manner that a line through S1 and S2 
	  //       divides the intersection plane.
	  tmpPts[0] = (*intersecPts)[3];

	  surfQuad->scaleSurfaceQuadrature(tmpPts);

	  for (val = iq = 0; iq < numQuadPts; ++iq) {
	    val += surfQuad->getWeight(iq)*(*f)(surfQuad->getLambda(iq));
	  }
	  int_val += calcSurfaceDet(loc_elInfo, tmpPts) * val;
	}
      }

      loc_elInfo = stack.traverseNext(loc_elInfo);

    }  // end of: mesh traverse

   /// delete surfQuad;  // deleting surfQuad prevents multiple executions of integrate_onZeroLs, so better don't delete it here

    return int_val;
  }

  double
  CFE_Integration::calcSurfaceDet(ElInfo *loc_elInfo,
				  VectorOfFixVecs<DimVec<double> > &surfVert)
  {
    double surfDet;
    int dim = surfVert[0].getSize() - 1;
    FixVec<WorldVector<double>, VERTEX> worldCoords(dim - 1);
  
    // transform barycentric coords to world coords
    for (int i = 0; i < dim; i++) {
      loc_elInfo->coordToWorld(surfVert[i], worldCoords[i]);
    }
  
    // calculate determinant for surface
    surfDet = loc_elInfo->calcDet(worldCoords);

    return surfDet;
  }

}
