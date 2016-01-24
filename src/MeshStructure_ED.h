/** \file MeshStructure.h */

#pragma once

#include <vector>

#include "MeshStructure.h"
#include "ElementData.h"
#include "Element.h"

namespace AMDiS
{

#define AMDIS_MESH_STRUCTURE 10

  class MeshStructure_ED : public ElementData
  {
  public:
    MeshStructure_ED(ElementData* decorated = NULL)
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

      if (!structure->isLeafElement())
      {
        MeshStructure* structure1 = new MeshStructure();
        structure->skipBranch(structure1);
        structure1->commit();
        MeshStructure_ED* elData1 = new MeshStructure_ED(child1->getElementData());
        elData1->setStructure(structure1);
        child1->setElementData(elData1);
      }
      else
      {
        structure->nextElement();
      }

      if (!structure->isLeafElement())
      {
        MeshStructure* structure2 = new MeshStructure();
        structure->skipBranch(structure2);
        structure2->commit();
        MeshStructure_ED* elData2 = new MeshStructure_ED(child2->getElementData());
        elData2->setStructure(structure2);
        child2->setElementData(elData2);
      }

      return true;
    }

    virtual int getTypeID() const
    {
      return AMDIS_MESH_STRUCTURE;
    }

    virtual bool isOfType(int typeID) const
    {
      if (typeID == AMDIS_MESH_STRUCTURE)
        return true;

      return false;
    }

    void setStructure(MeshStructure* s)
    {
      structure = s;
    }

    MeshStructure* getStructure()
    {
      return structure;
    }

  protected:
    MeshStructure* structure;
  };

} // end namespace AMDiS
