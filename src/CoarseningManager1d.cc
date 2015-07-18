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


#include "CoarseningManager1d.h"
#include "Mesh.h"
#include "AdaptStationary.h"
#include "AdaptInstationary.h"
#include "Traverse.h"
#include "MacroElement.h"
#include "RCNeighbourList.h"
#include "FixVec.h"
#include "DOFIndexed.h"

namespace AMDiS {

  Flag CoarseningManager1d::coarsenMesh(Mesh *aMesh)
  {
    mesh = aMesh;

    std::deque<MacroElement*>::const_iterator mel;
    int n_elements = mesh->getNumberOfLeaves();

    for (mel = mesh->firstMacroElement(); 
	 mel != mesh->endOfMacroElements(); 
	 mel++) 
      {
	coarsenRecursive(dynamic_cast<Line*>(const_cast<Element*>((*mel)->getElement())));
      }

    return(mesh->getNumberOfLeaves() < n_elements ? MESH_COARSENED : Flag(0));
  }

  int CoarseningManager1d::coarsenRecursive(Line *parent)
  {
    FUNCNAME("CoarseningManager1d::coarsenRecursive");
    int  mark;

    INFO(0,2)("\n");
    if (parent->isLeaf())  // leaf element, clean element marker
      {
	mark = parent->getMark();
	parent->setMark(0);
	return(mark);
      }
    else {                 // no leaf element, first coarse children
      Line *child[2];

      int  mark0 = coarsenRecursive(dynamic_cast<Line*>(const_cast<Element*>(parent->getChild(0))));
      int  mark1 = coarsenRecursive(dynamic_cast<Line*>(const_cast<Element*>( parent->getChild(1))));

      mark = std::max(mark0, mark1);

      child[0] = dynamic_cast<Line*>(const_cast<Element*>( parent->getChild(0)));
      child[1] = dynamic_cast<Line*>(const_cast<Element*>( parent->getChild(1))); 

      if (mark >= 0) {     // element must not be coarsend
	parent->setMark(0);
	if (child[0]) {
	  child[0]->setMark(0);
	  child[1]->setMark(0);
	}

	return(0);
      }

      /*--------------------------------------------------------------------------*/
      /*--- and now coarsen child[0] and child[1] into parent                  ---*/
      /*--------------------------------------------------------------------------*/
  
      /*--------------------------------------------------------------------------*/
      /*--- hand DOFs from children to parent, and add DOF at center           ---*/
      /*--------------------------------------------------------------------------*/

      if (mesh->getNumberOfDofs(VERTEX) && !mesh->queryCoarseDOFs()) {
	int node = mesh->getNode(VERTEX);
	parent->setDof(node+0, const_cast<DegreeOfFreedom*>( child[0]->getDof(node+0)));
	parent->setDof(node+1, const_cast<DegreeOfFreedom*>( child[1]->getDof(node+1)));
      }

      if (mesh->getNumberOfDofs(CENTER) && !mesh->queryCoarseDOFs()) {
	parent->setDof(mesh->getNode(CENTER), mesh->getDof(CENTER));
      }

      /*--------------------------------------------------------------------------*/
      /*  restrict dof vectors to the parents on the patch                        */
      /*--------------------------------------------------------------------------*/


      RCNeighbourList coarsenList(1); // = {{nil, 0, 0}};
      coarsenList.setElement(0, parent);
      int iadmin;
      int nrAdmin = mesh->getNumberOfDOFAdmin();
      for(iadmin = 0; iadmin < nrAdmin; iadmin++) {
	std::list<DOFIndexedBase*>::iterator it;
	DOFAdmin* admin = const_cast<DOFAdmin*>(&(mesh->getDofAdmin(iadmin)));
	std::list<DOFIndexedBase*>::iterator end = admin->endDOFIndexed();
	for(it = admin->beginDOFIndexed(); it != end; ++it)
	  (*it)->coarseRestrict(coarsenList, 1);
      }
 

      /*--------------------------------------------------------------------------*/
      /*--- remove all DOFs of children that are not used anymore              ---*/
      /*--------------------------------------------------------------------------*/

      if (mesh->getNumberOfDofs(VERTEX))    /*---  midpoint of parent          ---*/
	{
	  mesh->freeDof(const_cast<DegreeOfFreedom*>( child[1]->getDof(mesh->getNode(VERTEX))), VERTEX);
	}

      if (mesh->getNumberOfDofs(CENTER))    /*--- center of the children       ---*/
	{
	  mesh->freeDof(const_cast<DegreeOfFreedom*>( child[0]->getDof(mesh->getNode(CENTER))), CENTER);
	  mesh->freeDof(const_cast<DegreeOfFreedom*>( child[1]->getDof(mesh->getNode(CENTER))), CENTER);
	}

      parent->coarsenElementData(child[0], child[1]);

      parent->setFirstChild(NULL);
      parent->setSecondChild(NULL);

      mesh->freeElement(child[0]);
      mesh->freeElement(child[1]);

      mesh->incrementNumberOfLeaves(-1);
      mesh->incrementNumberOfElements(-2);
      mesh->incrementNumberOfVertices(-1);

      mark++;
      parent->setMark(std::min(mark,0));

      return(parent->getMark());
    }

    return(0);  /*--- statement never reached                              ---*/
  }

}
