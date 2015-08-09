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


#include <boost/numeric/mtl/mtl.hpp>
#include "CompositeFEMOperator.h"
#include "SubElementAssembler.h"
#include "SubElInfo.h"
#include "SubPolytope.h"

namespace compositeFEM {

using namespace std;
using namespace AMDiS;

void CompositeFEMOperator::getElementMatrix(const ElInfo *elInfo, 
					    ElementMatrix& userMat, 
					    double factor)
{
  FUNCNAME("CompositeFEMOperator::getElementMatrix");

  VectorOfFixVecs<DimVec<double> > *intersecPoints = NULL;
  SubPolytope *subPolytope = NULL;
  double levelSetSubPolytope;
  DimVec<double> subElVertexBarCoords(elInfo->getMesh()->getDim());

  /**
   * Get element status. Does element lie completely inside the integration 
   * domain, completely outside of the integration domain or is it 
   * intersected by the boundary ?
   */
  elStatus = elLS->createElementLevelSet(elInfo);

  /**
   * element status == completely inside or outside  
   *                                       --->  take the "normal" 
   *                                             integration routine
   *                                             Operator::getElementMatrix
   * element status == lies on boundary  ---> integration on subpolytopes and 
   *                                          subelements
   */
  if (elStatus == ElementLevelSet::LEVEL_SET_INTERIOR  ||  
      elStatus == ElementLevelSet::LEVEL_SET_EXTERIOR) {

    elLS->setLevelSetDomain(elStatus);
    Operator::getElementMatrix(elInfo, userMat, factor);
    return;
  }

  /***************************************************************************
   * Integration on intersected element.
   *
   * The integral is calculated as the sum of integrals on the two 
   * subpolytopes given by the intersection. 
   * We only calculate the integral on one of the subpolytopes. The 
   * integral on the second subpolytope then is the difference between the 
   * integral on the complete element and the integral on the first 
   * subpolytope.
   */

  if(!subElementAssembler) {
    subElementAssembler = new SubElementAssembler(this, 
						  rowFeSpace, 
						  colFeSpace);
  }

  // Get intersection points.
  intersecPoints = elLS->getElIntersecPoints();
  subPolytope = new SubPolytope(elInfo, 
				intersecPoints, 
				elLS->getNumElIntersecPoints());
  
  /**
   * Calculate integral on element.
   *
   * Whether a subpolytope lies inside or outside the integration domain is 
   * decided using the level set of the first vertex in the first subelement 
   * of the subpolytope. (The subelements of a subpolytope are created in 
   * such a way that this vertex always is a vertex of the element 
   * and not an intersection point. Thus the level set of this vertex really 
   * is unequal to zero.)
   */

  /**
   * Integration on subPolytope.
   */
  subElVertexBarCoords = subPolytope->getSubElement(0)->getLambda(0);
  levelSetSubPolytope = elLS->getVertexPos(
			(const DimVec<double>) subElVertexBarCoords);

  if (levelSetSubPolytope < 0) {
    elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);
  }
  else if (levelSetSubPolytope > 0) {
    elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_EXTERIOR);
  }
  else {
    ERROR_EXIT("cannot get position of subpolytope\n");
  }

  ElementMatrix subPolMat1(subElementAssembler->getNRow(),
			   subElementAssembler->getNCol());
  set_to_zero(subPolMat1);
  subElementAssembler->getSubPolytopeMatrix(subPolytope,
					    subElementAssembler,
					    elInfo,
					    subPolMat1);  

  /**
   * Integration on second subpolytope produced by the intersection.
   */
  ElementMatrix elMat(subElementAssembler->getNRow(),
		      subElementAssembler->getNCol());
  set_to_zero(elMat);
  ElementMatrix subPolMat2(subElementAssembler->getNRow(),
			   subElementAssembler->getNCol());
  set_to_zero(subPolMat2);

  if (!assembler.get()) {
    Assembler *aptr = new StandardAssembler(this, NULL, NULL, NULL, NULL, rowFeSpace, colFeSpace);
    assembler.set(aptr);
  }

  if (elLS->getLevelSetDomain() == 
      ElementLevelSet::LEVEL_SET_INTERIOR) {
    elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_EXTERIOR);
  }
  else {
    elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);
  }

  assembler.get()->calculateElementMatrix(elInfo, elMat, 1.0);
  subElementAssembler->getSubPolytopeMatrix(subPolytope,
					    subElementAssembler,
					    elInfo,
					    subPolMat2);

  elMat -= subPolMat2;

  // Get integral on element as sum of the two integrals on subpolytopes.
  elMat += subPolMat1;

  // Add integral to userMat.
  userMat += factor * elMat;

  // Free data
  delete subPolytope;
}


void CompositeFEMOperator::getElementVector(const ElInfo *elInfo, 
					    DenseVector<double>& userVec, 
					    double factor)
{
  FUNCNAME("CompositeFEMOperator::getElementVector");

  VectorOfFixVecs<DimVec<double> >*intersecPoints = NULL;
  SubPolytope *subPolytope = NULL;
  double levelSetSubPolytope;
  DimVec<double> subElVertexBarCoords(elInfo->getMesh()->getDim());

  /**
   * Get element status. Does element lie completely inside the integration 
   * domain, completely outside of the integration domain or is it 
   * intersected by the boundary ?
   */
  elStatus = elLS->createElementLevelSet(elInfo);

  /**
   * element status == completely inside or outside  
   *                                        --->  take the "normal" 
   *                                              integration routine  
   *                                              Operator::getElementVector
   * element status == lies on boundary  ---> integration on subpolytopes and 
   *                                          subelements
   */
  if (elStatus == ElementLevelSet::LEVEL_SET_INTERIOR  ||  
      elStatus == ElementLevelSet::LEVEL_SET_EXTERIOR) {

    elLS->setLevelSetDomain(elStatus);
    Operator::getElementVector(elInfo, userVec, factor);
    return;
  }

  /*********************************************************************************
   * Integration on intersected element.
   *
   * The integral is calculated as the sum of integrals on the two 
   * subpolytopes given by the intersection. 
   * We only calculate the integral on one of the subpolytopes. The integral 
   * on the second subpolytope then is the difference between the integral on 
   * the complete element and the integral on the first subpolytope.
   */

  if(!subElementAssembler) {
    subElementAssembler = new SubElementAssembler(this, 
						  rowFeSpace, 
						  colFeSpace);
  }

  /**
   * Get intersection points.
   */
  intersecPoints = elLS->getElIntersecPoints();
  subPolytope = new SubPolytope(elInfo, 
				intersecPoints, 
				elLS->getNumElIntersecPoints());

  /**
   * Calculate integral on element.
   *
   * Whether a subpolytope lies inside or outside the integration domain is 
   * decided using the level set of the first vertex in the first subelement 
   * of the subpolytope. (The subelements of a subpolytope are created in 
   * such a way that this vertex is always a vertex of the element and not 
   * an intersection point. Thus the level set of this vertex really is 
   * unequal to zero.)
   */

  /**
   * Integration on ubPolytope.
   */
  subElVertexBarCoords = subPolytope->getSubElement(0)->getLambda(0);
  levelSetSubPolytope = elLS->getVertexPos(
			(const DimVec<double>) subElVertexBarCoords);

  if (levelSetSubPolytope < 0) {
    elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);
  }
  else if (levelSetSubPolytope > 0) {
    elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_EXTERIOR);
  }
  else {
    ERROR_EXIT("cannot get position of subpolytope\n");
  }

  DenseVector<double> subPolVec1(subElementAssembler->getNRow());
  set_to_zero(subPolVec1);
  subElementAssembler->getSubPolytopeVector(subPolytope,
					    subElementAssembler,
					    elInfo,
					    subPolVec1);  

  // Integration on second subpolytope produced by the intersection.
  DenseVector<double> elVec(subElementAssembler->getNRow());
  set_to_zero(elVec);
  DenseVector<double> subPolVec2(subElementAssembler->getNRow());
  set_to_zero(subPolVec2);

  if (!assembler.get()) {
    Assembler *aptr = new StandardAssembler(this, NULL, NULL, NULL, NULL, rowFeSpace, colFeSpace);
    assembler.set(aptr);      
  }

  if (elLS->getLevelSetDomain() == 
      ElementLevelSet::LEVEL_SET_INTERIOR) {
    elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_EXTERIOR);
  } else {
    elLS->setLevelSetDomain(ElementLevelSet::LEVEL_SET_INTERIOR);
  }

  assembler.get()->calculateElementVector(elInfo, elVec, 1.0);
  subElementAssembler->getSubPolytopeVector(subPolytope,
					    subElementAssembler,
					    elInfo,
					    subPolVec2);

  elVec -= subPolVec2;

  // Get integral on element as sum of the two integrals on subpolytopes.
  elVec += subPolVec1;

  // Add integral to userVec.
  userVec += factor * elVec;

  // Free data
  delete subPolytope;
}

}
