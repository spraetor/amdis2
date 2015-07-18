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
#ifndef AMDIS_PERIODICCONNECTION_H
#define AMDIS_PERIODICCONNECTION_H

		       Element *el1,
		       DimVec<int> &vertices0,
		       DimVec<int> &vertices1) 
    {
//
// Software License for AMDiS
//
// Copyright (c) 2010 Dresden University of Technology 
// All rights reserved.
// Authors: Simon Vey, Thomas Witkowski et al.
//
// This file is part of AMDiS
//
// See also license.opensource.txt in the distribution.


      FUNCNAME("PeriodicConnection::PeriodicConnection()");
      TEST_EXIT(el0)("element0 not set\n");
      TEST_EXIT(el1)("element1 not set\n");
      TEST_EXIT(el0->getMesh()->getDim() == dim)("invalid dim of el0\n");
      TEST_EXIT(el1->getMesh()->getDim() == dim)("invalid dim of el1\n");
      TEST_EXIT(vertices0.getSize() == dim)("invalid size of vertices0\n");
      TEST_EXIT(vertices1.getSize() == dim)("invalid size of vertices1\n");
      int i;
      for(i = 0; i < dim + 1; i++) {
	TEST_EXIT(vertices0[i] < dim + 1 && vertices0[i] >= 0)
	  ("invalid value in vertices0\n");
	TEST_EXIT(vertices1[i] < dim + 1 && vertices1[i] >= 0)
	  ("invalid value in vertices1\n");
      }

      el_[0] = el0;
      el_[1] = el1;
      vertices_[0] = new DimVec<int>(vertices0);
      vertices_[1] = new DimVec<int>(vertices1);
      child_[0] = child_[1] = NULL;
    };

    ~PeriodicConnection() {
      delete vertices_[0];
      delete vertices_[1];
    };

    Element *getElement(int index) { 
      FUNCNAME("PeriodicConnection::getElement()");
      TEST_EXIT(index >= 0 && index <= 2)("invalid index\n");
      return el_[index];
    };

    const DimVec<int>& getVertices(int index) { 
      FUNCNAME("PeriodicConnection::getVertices()");
      TEST_EXIT(index >= 0 && index <= 2)("invalid index\n");
      return *(vertices_[index]);
    };

    bool removeElement(int index) {
      FUNCNAME("PeriodicConnection::removeElement()");
      TEST_EXIT(index >= 0 && index <= 2)("invalid index\n");
      TEST_EXIT(el_[index])("element already removed\n");
      el_[index] = NULL;
      return (el_[abs(index-1)] == NULL);
    };

    bool refineElement(int index, 
		       PeriodicConnection **child0,
		       PeriodicConnection **child1) 
    {
      FUNCNAME("PeriodicConnection::refineElement()");

      if(child_[0] || child_[1]) {
	TEST_EXIT(child_[0] && child_[1])
	  ("only one child\n");

	TEST_EXIT(el_[abs(index-1)] == NULL)
	  ("connection already refined but other element != NULL\n");

	*child0 = child_[0];
	*child1 = child_[1];

	return true;
      }

    

      return false;
    };

  protected:
    Element *el_[2];
    DimVec<int> *vertices_[2];
    PeriodicConnection *child_[2];
  };

}

#endif
