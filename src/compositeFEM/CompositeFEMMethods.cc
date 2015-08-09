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


#include "BasisFunction.h"
#include "DOFAdmin.h"
#include "Element.h"
#include "ElInfo.h"
#include "FiniteElemSpace.h"
#include "Global.h"
#include "Mesh.h"
#include "Initfile.h"
#include "Traverse.h"

#include "CompositeFEMMethods.h"

namespace compositeFEM {

using namespace AMDiS;
using namespace std;

void CompositeFEMMethods::setPosLsToVal(DOFVector<double> *dof,
					const double &val,
					const DOFVector<double> *lsFct_dof)
{
  DOFVector<double>::Iterator it_dof(dof, USED_DOFS);
  DOFVector<double>::Iterator it_lsFct_dof(
                     const_cast<DOFVector<double> *>(lsFct_dof), 
		     USED_DOFS);
  for (it_dof.reset(), it_lsFct_dof.reset();
       !it_dof.end();
       ++it_dof, ++it_lsFct_dof) {

    // Is vertex in domain with positive level set function values ?
    if (*it_lsFct_dof > 0)
      *it_dof = val;
  }
}


void CompositeFEMMethods::setPosLsToFct(DOFVector<double> *dof,
					const AbstractFunction<double, WorldVector<double> > *fct,
					const DOFVector<double> *lsFct_dof)
{
  const BasisFunction *basisFcts = dof->getFeSpace()->getBasisFcts();
  const DOFAdmin *admin = dof->getFeSpace()->getAdmin();
  const int dim = dof->getFeSpace()->getMesh()->getDim();

  TraverseStack stack;
  ElInfo *elInfo = stack.traverseFirst(dof->getFeSpace()->getMesh(), -1, 
				       Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
  DenseVector<double> locVec(basisFcts->getNumber());
  std::vector<DegreeOfFreedom> locInd(basisFcts->getNumber());

  while (elInfo) {
    const Element *el = elInfo->getElement();
    
    // Get level set function values for all vertices of element.
    lsFct_dof->getLocalVector(el, locVec);

    // Get dof indices of vertices.
    basisFcts->getLocalIndices(el, admin, locInd);

    for (int i = 0; i <= dim; i++) {
      // Is vertex in domain with positive level set function values ?
      if (locVec[i] > 0)
	(*dof)[locInd[i]] = (*fct)(elInfo->getCoord(i));
    }

    elInfo = stack.traverseNext(elInfo);
  }
}


void CompositeFEMMethods::printBoundaryElements(const std::string fn_str,
						ElementLevelSet *elLS,
						FiniteElemSpace *feSpace)
{
  int dim = feSpace->getMesh()->getDim();
  std::string fn_main;
  std::string fn;
  Parameters::get(fn_str + "->filename", fn_main);
  fn = fn_main + "boundary_elements";
  int elStatus;

  ofstream boundaryOut(fn.c_str());

  // ===== Traverse mesh and print boundary elements. =====
  TraverseStack stack;
  int boundEl_cntr = 0;
  WorldVector<double> coord;

  const int nBasFcts = feSpace->getBasisFcts()->getNumber();
  std::vector<DegreeOfFreedom> locInd(nBasFcts);

  ElInfo *loc_elInfo = stack.traverseFirst(feSpace->getMesh(),
					   -1, 
					   Mesh::CALL_LEAF_EL | 
					   Mesh::FILL_BOUND |
					   Mesh::FILL_COORDS);
  while (loc_elInfo) {

    // Get local indices of vertices.
    feSpace->getBasisFcts()->getLocalIndices(const_cast<Element *>(loc_elInfo->getElement()),
					     const_cast<DOFAdmin *>(feSpace->getAdmin()),
					     locInd);

    // Get element status.
    elStatus = elLS->createElementLevelSet(loc_elInfo);

    // Is element cut by the interface ?
    if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {
      ++boundEl_cntr;

      boundaryOut << loc_elInfo->getElement()->getIndex() << ": \t";

      for (int i=0; i<=dim; ++i) {

	coord = loc_elInfo->getCoord(i);

	if (i > 0) {
	  boundaryOut << "\n\t";
	}

	boundaryOut << "\t" << coord[0];
	for (int j=1; j<dim; ++j) {
	  boundaryOut << "; " << coord[j] ;
	}
	boundaryOut << " (" << locInd[i] << ")";
      }
      boundaryOut << "\n";
    }

    loc_elInfo = stack.traverseNext(loc_elInfo);

  }  // end of: mesh traverse

  boundaryOut << "\nNumber of boundary elements: \t" << boundEl_cntr << "\n";

  boundaryOut.close();
}

}
