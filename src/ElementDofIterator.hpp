/** \file ElementDofIterator.h */

#pragma once

#include "AMDiS_fwd.hpp"
#include "Global.hpp"
#include "Mesh.hpp"

namespace AMDiS
{
  /** \brief
   * This class implements an iterator to iterate over all DOFs of one element
   * independet of dimension and the degree of basis functions.
   *
   * Should be used in the following way:
   *
   *    ElementDofIterator elDofIter(feSpace);
   *    elDofIter.reset(el);
   *    do {
   *       elDofIter.getDof();
   *    } while (elDofIter.next());
   *
   */
  class ElementDofIterator
  {
  public:
    /// Constructor.
    ElementDofIterator(const FiniteElemSpace* feSpace, bool inOrderPos = false);

    /// Start a new traverse with the given element.
    void reset(const Element* el);

    /// Go to next DOF. Returns false, if there is no DOF anymore.
    // ATTENTION: works only, if #Vertex-DOFs == 1
    bool next();

    bool nextStrict()
    {
      dofPos = nDofs;
      return next();
    }

    /// Returns the DOF index of the current DOF.
    DegreeOfFreedom getDof() const;

    /// Returns a pointer to the current DOF.
    const DegreeOfFreedom* getDofPtr() const;

    /// Returns a pointer to the starting position of the current DOF
    /// array. Makes only sence, if \ref nextStrict() is used for traverse.
    const DegreeOfFreedom* getBaseDof() const
    {
      return dofs[node0 + elementPos];
    }

    /// Returns \ref pos, the current position (vertex, edge, face) of
    /// the traverse.
    int getCurrentPos() const
    {
      return pos;
    }

    /// Returns \ref node0
    int getCurrentNode0() const
    {
      return node0;
    }

    /// Returns \ref elementPos, the number of vertex, edge or face that
    /// is traversed.
    int getCurrentElementPos() const
    {
      return elementPos;
    }

    GeoIndex getPosIndex() const
    {
      return posIndex;
    }


  protected:
    /// The DOF admin for which DOFs should be traversed.
    const DOFAdmin* admin;

    const BasisFunction* basisFcts;

    /// Pointer to the DOFs that should be traversed.
    const DegreeOfFreedom** dofs;

    /// Mesh on which the element is defined.
    Mesh* mesh;

    /// Dimension of the mesh.
    int dim;

    bool inOrder;

    int* orderPosition;

    const Element* element;

    /// Current position (i.e., vertex, edge, face) of the traverse.
    int pos;

    /// Dimension dependent geo index of the current position in traverse.
    GeoIndex posIndex;

    /// Number of DOFs at the current traverse position. Examples: independent of
    /// dimension and  degree of basis functions there is only one DOF per vertex.
    /// But in 2D and with 3rd degree lagrange basis functions there are two
    /// DOFs per edge.
    int nDofs;

    /// Displacement of DOF indices. Used if more than one DOF admin is defined
    /// on the mesh.
    int n0;

    /// DOF index of the first DOF at this geo index position.
    int node0;

    /// Number of elements in the current geo position. Examples: 3 vertices in
    /// 2d, 1 face in 2d, 4 faces in 3d, etc.
    int nElements;

    /// Current element, i.e., ith vertex, edge or face, that is traversed.
    int elementPos;

    /// Currrent DOF that is traversed on the current element;
    int dofPos;
  };

} // end namespace AMDiS
