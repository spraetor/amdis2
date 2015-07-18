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


#include "ComponentTraverseInfo.h"

namespace AMDiS {

  void SingleComponentInfo::updateStatus()
  {
    if (rowFeSpace == NULL) {
      status = SingleComponentInfo::EMPTY;
      return;
    }

    if (colFeSpace == NULL || 
	(colFeSpace != NULL && rowFeSpace->getMesh() == colFeSpace->getMesh())) {
      if (auxFeSpaces.size() == 0) {
	status = SingleComponentInfo::EQ_SPACES_NO_AUX;
      } else {
	status = SingleComponentInfo::EQ_SPACES_WITH_AUX;

	for (std::set<const FiniteElemSpace*>::iterator it = auxFeSpaces.begin();
	     it != auxFeSpaces.end(); ++it) {
	  if ((*it)->getMesh() != rowFeSpace->getMesh()) {     
	    status = SingleComponentInfo::EQ_SPACES_WITH_DIF_AUX;
	    break;
	  }
	}
      }
    } else {
      if (auxFeSpaces.size() == 0) {
	status = SingleComponentInfo::DIF_SPACES_NO_AUX;
      } else {
	status = SingleComponentInfo::DIF_SPACES_WITH_AUX;
	for (std::set<const FiniteElemSpace*>::iterator it = auxFeSpaces.begin();
	     it != auxFeSpaces.end(); ++it) {
	  if ((*it)->getMesh() != rowFeSpace->getMesh() &&
	      (*it)->getMesh() != colFeSpace->getMesh()) {
	    status = SingleComponentInfo::DIF_SPACES_WITH_DIF_AUX;
	    break;
	  }
	}
      }	
    }    
  }

  
  const FiniteElemSpace* ComponentTraverseInfo::getRowFeSpace(int row) const
  {
    FUNCNAME_DBG("ComponentTraverseInfo::getRowFeSpace()");
    
    TEST_EXIT_DBG(row < nComponents)("No component traverse info for this row!\n");
    TEST_EXIT_DBG(matrixComponents[row][row].getRowFeSpace() ==
		  matrixComponents[row][row].getColFeSpace())
      ("Should not happen!\n");
    
    return matrixComponents[row][row].getRowFeSpace();      
  }


  const FiniteElemSpace* ComponentTraverseInfo::getNonRowFeSpace(int row) const
  {
    FUNCNAME_DBG("ComponentTraverseInfo::getNonRowFeSpace()");
    TEST_EXIT_DBG(row < nComponents)("No component traverse info for this row!\n");
    
    const FiniteElemSpace* rowFeSpace = getRowFeSpace(row);
    TEST_EXIT_DBG(rowFeSpace != NULL)("No row FE space!\n");
    
    for (int i = 0; i < nComponents; i++) {
      if (matrixComponents[row][i].getColFeSpace() && matrixComponents[row][i].getColFeSpace() != rowFeSpace)
	return matrixComponents[row][i].getColFeSpace();
      if (matrixComponents[row][i].getAuxFeSpace() && matrixComponents[row][i].getAuxFeSpace() != rowFeSpace)
	return matrixComponents[row][i].getAuxFeSpace();
    }
    
    if (vectorComponents[row].getAuxFeSpace() != rowFeSpace)
      return vectorComponents[row].getAuxFeSpace();
    
    return NULL;
  }

}
