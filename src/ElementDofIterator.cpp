#include "ElementDofIterator.hpp"
#include "Mesh.hpp"
#include "DOFAdmin.hpp"
#include "Element.hpp"
#include "BasisFunction.hpp"
#include "FiniteElemSpace.hpp"

namespace AMDiS
{

  ElementDofIterator::ElementDofIterator(const FiniteElemSpace* feSpace, bool inOrderPos)
    : admin(feSpace->getAdmin()),
      basisFcts(feSpace->getBasisFcts()),
      mesh(feSpace->getMesh()),
      dim(mesh->getDim()),
      inOrder(inOrderPos)
  { }


  const DegreeOfFreedom*
  ElementDofIterator::getDofPtr() const
  {
    if (inOrder)
      return &dofs[node0 + elementPos][n0 + orderPosition[dofPos]];
    else
      return &dofs[node0 + elementPos][n0 + dofPos];
  }


  DegreeOfFreedom
  ElementDofIterator::getDof() const
  {
    if (inOrder)
      return dofs[node0 + elementPos][n0 + orderPosition[dofPos]];
    else
      return dofs[node0 + elementPos][n0 + dofPos];
  }

  void ElementDofIterator::reset(const Element* el)
  {
    FUNCNAME_DBG("ElementDofIterator::reset()");

    TEST_EXIT_DBG(el->getMesh() == mesh)
    ("Mesh and element does not fit together!\n");
    TEST_EXIT_DBG(el)("No element!\n");

    element = el;
    dofs = element->getDof();

    // Start with vertices.
    pos = 0;
    elementPos = 0;
    dofPos = 0;

    // Get geo index of vertices in the given dimension.
    posIndex = INDEX_OF_DIM(pos, dim);
    // Get number of DOFs per vertex (should be one in all cases).
    nDofs = admin->getNumberOfDofs(posIndex);

    TEST_EXIT_DBG(nDofs != 0)("Mh, I've to think about this situation!\n");

    // Calculate displacement. Is used if there is more than one DOF admin
    // on the mesh.
    n0 = admin->getNumberOfPreDofs(posIndex);
    // Get first DOF index position for vertices.
    node0 = mesh->getNode(posIndex);
    // Get number of vertices in this dimension.
    nElements = Global::getGeo(posIndex, mesh->getDim());

    if (inOrder)
      orderPosition = basisFcts->orderOfPositionIndices(element, posIndex, 0);
  }


  bool ElementDofIterator::next()
  {
    // First iterate over the DOFs of one element (vertex, edge, face).
    dofPos++;

    if (dofPos >= nDofs)
    {
      // We are finished with all DOFs of on element. Go to the next one.
      dofPos = 0;
      elementPos++;

      if (elementPos >= nElements)
      {
        // We are finished with all element.
        elementPos = 0;

        // If we have iterated over all positions, we can finish the iteration.
        if (pos >= dim)
          return false;

        // Increase position, i.e., go from vertices to edges to faces and search
        // for the next position with DOFs.
        do
        {
          pos++;
          // Get geo index posistion.
          posIndex = INDEX_OF_DIM(pos, dim);
          // Get number of DOFs in this position.
          nDofs = admin->getNumberOfDofs(posIndex);
        }
        while (nDofs == 0 && pos < dim);

        if (nDofs > 0 && pos <= dim)
        {
          // We have found on more position with DOFs.

          // Get number of elements in this position, i.e, the number of
          // vertices, edges and faces in the given dimension.
          nElements = Global::getGeo(posIndex, dim);

          // Calculate displacement. Is used if there is more than one DOF
          // admin on the mesh.
          n0 = admin->getNumberOfPreDofs(posIndex);

          // Get first DOF index position for the geo index position.
          node0 = mesh->getNode(posIndex);

          if (inOrder)
            orderPosition =
              basisFcts->orderOfPositionIndices(element, posIndex, 0);

        }
        else
        {
          // That's all, we jave traversed all DOFs of the mesh element.
          return false;
        }

      }
      else
      {
        if (inOrder)
          orderPosition =
            basisFcts->orderOfPositionIndices(element, posIndex, elementPos);
      }
    }

    return true;
  }

} // end namespace AMDiS
