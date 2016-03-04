#include "Parametric.hpp"

#include "DOFVector.hpp"
#include "ElInfo.hpp"
#include "FixVec.hpp"

namespace AMDiS
{
  ElInfo* ParametricFirstOrder::addParametricInfo(ElInfo* elInfo)
  {
    elInfo->setParametric(true);
    const int dow = Global::getGeo(WORLD);
    Element* element = elInfo->getElement();
    const DegreeOfFreedom** dof = element->getDof();

    for (int i = 0; i < elInfo->getElement()->getGeo(VERTEX); i++)
    {
      if (elInfo->getFillFlag().isSet(Mesh::FILL_COORDS))
        for (int j = 0; j < dow; j++)
          elInfo->getCoord(i)[j] = (*(*dofCoords_)[j])[dof[i][0]];

      if (elInfo->getFillFlag().isSet(Mesh::FILL_OPP_COORDS))
      {
        TEST_EXIT(elInfo->getFillFlag().isSet(Mesh::FILL_NEIGH))
        ("FILL_NEIGH not set\n");

        if (elInfo->getNeighbour(i))
        {
          const DegreeOfFreedom** neighDof = elInfo->getNeighbour(i)->getDof();
          for (int j = 0; j < dow; j++)
            elInfo->getOppCoord(i)[j] =
              (*(*dofCoords_)[j])[neighDof[elInfo->getOppVertex(i)][0]];
        }
      }
    }

    return elInfo;
  }

  ElInfo* ParametricFirstOrder::removeParametricInfo(ElInfo* elInfo)
  {
    elInfo->setParametric(false);
    return elInfo;
  }

} // end namespace AMDiS
