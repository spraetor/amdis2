#include "Triangle.hpp"

#include "CoarseningManager.hpp"
#include "DOFAdmin.hpp"
#include "ElementDofIterator.hpp"
#include "FiniteElemSpace.hpp"
#include "FixVec.hpp"
#include "Mesh.hpp"

using namespace std;

namespace AMDiS
{
  constexpr int Triangle::vertexOfEdge[3][2];
  constexpr int Triangle::sideOfChild[2][3];
  constexpr int Triangle::vertexOfParent[2][3];

  
  bool Triangle::hasSide(Element* sideElem) const
  {
    FUNCNAME("Triangle::hasSide");
    TEST_EXIT_DBG(sideElem->isLine())("called for sideElem-type != Line\n");
    ERROR_EXIT("not yet\n");
    return false;
  }


  int Triangle::getVertexOfPosition(GeoIndex position, int positionIndex,
                                    int vertexIndex) const
  {
    FUNCNAME("Triangle::getVertexOfPosition");
    switch(position)
    {
    case VERTEX:
      return positionIndex;
      break;
    case EDGE:
      return vertexOfEdge[positionIndex][vertexIndex];
      break;
    case CENTER:
      return vertexIndex;
      break;
    default:
      ERROR_EXIT("invalid position\n");
      return 0;
    }
  }


  void Triangle::sortFaceIndices(int face, FixVec<int, WORLD>& vec) const
  {
    // TODO: REMOVE STATIC
    static MatrixOfFixVecs<FixVec<int, WORLD>>* sorted_2d = NULL;

    if (sorted_2d == NULL)
    {
      sorted_2d = new MatrixOfFixVecs<FixVec<int, WORLD>>(2, 3, 2, NO_INIT);

      (*sorted_2d)[1][0][1] = (*sorted_2d)[1][1][0] =
                                (*sorted_2d)[2][0][0] = (*sorted_2d)[2][1][1] = 0;
      (*sorted_2d)[0][0][0] = (*sorted_2d)[0][1][1] =
                                (*sorted_2d)[2][0][1] = (*sorted_2d)[2][1][0] = 1;
      (*sorted_2d)[0][0][1] = (*sorted_2d)[0][1][0] =
                                (*sorted_2d)[1][0][0] = (*sorted_2d)[1][1][1] = 2;
    }

    const int* vof = vertexOfEdge[face];
    int no = ((dof[vof[0]][0] < dof[vof[1]][0]) ? 0 : 1);
    vec = (*sorted_2d)[face][no];
  }


  void Triangle::getNodeDofs(const FiniteElemSpace* feSpace,
                             BoundaryObject bound,
                             DofContainer& dofs,
                             bool baseDofPtr) const
  {
    FUNCNAME("Triangle::getNodeDofs()");

    // Get displacement if more than one FE space is defined on mesh.
    int n0 = (baseDofPtr ? 0 : feSpace->getAdmin()->getNumberOfPreDofs(VERTEX));

    if (bound.subObj == VERTEX)
    {
      dofs.push_back(&(dof[bound.ithObj][n0]));
      return;
    }

    TEST_EXIT_DBG(bound.subObj == EDGE)("This should not happen!\n");

    BoundaryObject nextBound = bound;

    switch (bound.ithObj)
    {
    case 0:
    {
      if (child[1] && child[1]->getFirstChild())
      {
        const DegreeOfFreedom** elDofs = child[1]->getFirstChild()->getDof();

        if (bound.reverseMode)
        {
          nextBound.ithObj = 1;
          child[1]->getSecondChild()->getNodeDofs(feSpace, nextBound, dofs,
                                                  baseDofPtr);
          dofs.push_back(&(elDofs[2][n0]));
          nextBound.ithObj = 0;
          child[1]->getFirstChild()->getNodeDofs(feSpace, nextBound, dofs,
                                                 baseDofPtr);
        }
        else
        {
          nextBound.ithObj = 0;
          child[1]->getFirstChild()->getNodeDofs(feSpace, nextBound, dofs,
                                                 baseDofPtr);
          dofs.push_back(&(elDofs[2][n0]));
          nextBound.ithObj = 1;
          child[1]->getSecondChild()->getNodeDofs(feSpace, nextBound, dofs,
                                                  baseDofPtr);
        }
      }
    }
    break;
    case 1:
    {
      if (child[0] && child[0]->getFirstChild())
      {
        const DegreeOfFreedom** elDofs = child[0]->getFirstChild()->getDof();

        if (bound.reverseMode)
        {
          nextBound.ithObj = 1;
          child[0]->getSecondChild()->getNodeDofs(feSpace, nextBound, dofs,
                                                  baseDofPtr);
          dofs.push_back(&(elDofs[2][n0]));
          nextBound.ithObj = 0;
          child[0]->getFirstChild()->getNodeDofs(feSpace, nextBound, dofs,
                                                 baseDofPtr);
        }
        else
        {
          nextBound.ithObj = 0;
          child[0]->getFirstChild()->getNodeDofs(feSpace, nextBound, dofs,
                                                 baseDofPtr);
          dofs.push_back(&(elDofs[2][n0]));
          nextBound.ithObj = 1;
          child[0]->getSecondChild()->getNodeDofs(feSpace, nextBound, dofs,
                                                  baseDofPtr);
        }
      }
    }
    break;
    case 2:
      if (child[0])
      {
        const DegreeOfFreedom** elDofs = child[0]->getDof();

        if (bound.reverseMode)
        {
          nextBound.ithObj = 1;
          child[1]->getNodeDofs(feSpace, nextBound, dofs, baseDofPtr);
          dofs.push_back(&(elDofs[2][n0]));
          nextBound.ithObj = 0;
          child[0]->getNodeDofs(feSpace, nextBound, dofs, baseDofPtr);
        }
        else
        {
          nextBound.ithObj = 0;
          child[0]->getNodeDofs(feSpace, nextBound, dofs, baseDofPtr);
          dofs.push_back(&(elDofs[2][n0]));
          nextBound.ithObj = 1;
          child[1]->getNodeDofs(feSpace, nextBound, dofs, baseDofPtr);
        }
      }
      break;
    default:
      ERROR_EXIT("Should never happen!\n");
    }
  }


  void Triangle::getHigherOrderDofs(const FiniteElemSpace* feSpace,
                                    BoundaryObject bound,
                                    DofContainer& dofs,
                                    bool baseDofPtr,
                                    vector<GeoIndex>* dofGeoIndex) const
  {
    FUNCNAME("Triange::getHigherOrderDofs()");

    if (bound.subObj == VERTEX)
      return;

    TEST_EXIT_DBG(bound.subObj == EDGE)("This should not happen!\n");

    bool addThisEdge = false;
    BoundaryObject nextBound = bound;

    switch (bound.ithObj)
    {
    case 0:
      if (child[1])
      {
        nextBound.ithObj = 2;
        child[1]->getHigherOrderDofs(feSpace, nextBound, dofs,
                                     baseDofPtr, dofGeoIndex);
      }
      else
      {
        addThisEdge = true;
      }

      break;
    case 1:
      if (child[0])
      {
        nextBound.ithObj = 2;
        child[0]->getHigherOrderDofs(feSpace, nextBound, dofs,
                                     baseDofPtr, dofGeoIndex);
      }
      else
      {
        addThisEdge = true;
      }

      break;
    case 2:
      if (child[0])
      {
        if (bound.reverseMode)
        {
          nextBound.ithObj = 1;
          child[1]->getHigherOrderDofs(feSpace, nextBound, dofs,
                                       baseDofPtr, dofGeoIndex);
          nextBound.ithObj = 0;
          child[0]->getHigherOrderDofs(feSpace, nextBound, dofs,
                                       baseDofPtr, dofGeoIndex);
        }
        else
        {
          nextBound.ithObj = 0;
          child[0]->getHigherOrderDofs(feSpace, nextBound, dofs,
                                       baseDofPtr, dofGeoIndex);
          nextBound.ithObj = 1;
          child[1]->getHigherOrderDofs(feSpace, nextBound, dofs,
                                       baseDofPtr, dofGeoIndex);
        }
      }
      else
      {
        addThisEdge = true;
      }

      break;
    default:
      ERROR_EXIT("Should never happen!\n");
    }

    if (addThisEdge)
    {
      DofContainer addDofs;
      ElementDofIterator elDofIter(feSpace, true);
      elDofIter.reset(this);

      if (baseDofPtr)
      {
        do
        {
          if (elDofIter.getPosIndex() == EDGE &&
              elDofIter.getCurrentElementPos() == bound.ithObj)
            addDofs.push_back(elDofIter.getBaseDof());
        }
        while (elDofIter.nextStrict());
      }
      else
      {
        do
        {
          if (elDofIter.getPosIndex() == EDGE &&
              elDofIter.getCurrentElementPos() == bound.ithObj)
            addDofs.push_back(elDofIter.getDofPtr());
        }
        while (elDofIter.next());
      }

      if (bound.reverseMode)
      {
        for (int i = int( addDofs.size() ) - 1; i >= 0; i--)
        {
          dofs.push_back(addDofs[i]);
          if (dofGeoIndex != NULL)
            dofGeoIndex->push_back(EDGE);
        }
      }
      else
      {
        for (size_t i = 0; i < addDofs.size(); i++)
        {
          dofs.push_back(addDofs[i]);
          if (dofGeoIndex != NULL)
            dofGeoIndex->push_back(EDGE);
        }
      }
    }
  }


  void Triangle::getSubBoundary(BoundaryObject bound,
                                vector<BoundaryObject>& subBound) const
  {
    TEST_EXIT_DBG(bound.subObj == EDGE)("This should not happen!\n");

    if (!child[0])
    {
      subBound.push_back(bound);
      return;
    }

    BoundaryObject nextBound0 = bound;
    BoundaryObject nextBound1 = bound;
    prepareNextBound(nextBound0, 0);
    prepareNextBound(nextBound1, 1);

    if (nextBound0.ithObj >= 0 && nextBound1.ithObj >= 0)
    {
      if (bound.reverseMode)
      {
        child[1]->getSubBoundary(nextBound1, subBound);
        child[0]->getSubBoundary(nextBound0, subBound);
      }
      else
      {
        child[0]->getSubBoundary(nextBound0, subBound); // TODO: check this!
        child[1]->getSubBoundary(nextBound1, subBound);
      }
    }
    else
    {
      if (nextBound0.ithObj >= 0)
        child[0]->getSubBoundary(nextBound0, subBound);

      if (nextBound1.ithObj >= 0)
        child[1]->getSubBoundary(nextBound1, subBound);
    }
  }


  void Triangle::prepareNextBound(BoundaryObject& bound, int ithChild) const
  {
    FUNCNAME_DBG("Triangle::prepareNextBound()");

    TEST_EXIT_DBG(bound.el == this)("Wrong element!\n");
    TEST_EXIT_DBG(child[0])("Has no child!\n");

    bound.ithObj = sideOfChild[ithChild][bound.ithObj];
    bound.el = child[ithChild];
    bound.elIndex =  child[ithChild]->getIndex();
  }

} // end namespace AMDiS
