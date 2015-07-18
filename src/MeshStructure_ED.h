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



/** \file MeshStructure.h */

#ifndef AMDIS_MESHSTRUCTUREED_H
#define AMDIS_MESHSTRUCTUREED_H

#include <vector>
#include "MeshStructure.h"
#include "ElementData.h"
#include "Element.h"

namespace AMDiS {

  const int MESH_STRUCTURE = 10;

  class MeshStructure_ED : public ElementData
  {
  public:
    MeshStructure_ED(ElementData *decorated = NULL) 
      : ElementData(decorated),
	structure(NULL)
    {}

    virtual ~MeshStructure_ED() 
    {
      delete structure;
    }

    virtual bool refineElementData(Element* parent, 
				   Element* child1,
				   Element* child2,
				   int elType)
    {
      FUNCNAME_DBG("MeshStructureED::refineElementData()");

      ElementData::refineElementData(parent, child1, child2, elType);

      structure->reset();

      TEST_EXIT_DBG(structure)("No structure!\n");
      TEST_EXIT_DBG(!structure->isLeafElement())("Don't store leaf structures!\n");

      structure->nextElement();

      if (!structure->isLeafElement()) {
	MeshStructure *structure1 = new MeshStructure();
	structure->skipBranch(structure1);
	structure1->commit();
	MeshStructure_ED *elData1 = new MeshStructure_ED(child1->getElementData());
	elData1->setStructure(structure1);
	child1->setElementData(elData1);
      } else {
	structure->nextElement();
      }

      if (!structure->isLeafElement()) {
	MeshStructure *structure2 = new MeshStructure();
	structure->skipBranch(structure2);
	structure2->commit();
	MeshStructure_ED *elData2 = new MeshStructure_ED(child2->getElementData());
	elData2->setStructure(structure2);
	child2->setElementData(elData2);
      }

      return true;
    }

    virtual int getTypeID() const 
    {
      return MESH_STRUCTURE;
    }

    virtual bool isOfType(int typeID) const 
    {
      if (typeID == MESH_STRUCTURE) 
	return true;

      return false;
    }

    inline void setStructure(MeshStructure *s) 
    {
      structure = s;
    }

    inline MeshStructure *getStructure() 
    {
      return structure;
    }

  protected:
    MeshStructure *structure;
  };

}

#endif
