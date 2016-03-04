#pragma once

#include <deque>
#include <vector>

#include "AMDiS_fwd.hpp"
#include "AMDiS_base.hpp"
#include "Global.hpp"

namespace AMDiS
{
  /** \ingroup Adaption
   * \brief
   * Stores information about coarsening and refinement patches. For refining
   * and coarsening we need information of the elements at the refinement and
   * coarsening edge. Thus, we have to collect all elements at this edge. In 2d
   * we have at most the current element and its neighbour across this edge, if
   * the edge is not part of the boundary. In 3d we have to loop around this
   * edge to collect all the elements. Every element at the edge has at most two
   * neighbours sharing the same edge. Defining an orientation for this edge,
   * we can define the right and left neighbour.
   * For every element at the refinement/coarsening edge we have an entry in a
   * vector. The elements of this vector build the refinement/coarsening patch.
   * In 2d the vector has length 2 and in 3d length mesh->getMaxEdgeNeigh()
   * since this is the maximal number of elements sharing the same edge in the
   * Mesh.
   */
  class RCNeighbourList
  {
  public:
    /// Constructs a RCNeighbourList of size maxEdgeNeigh
    RCNeighbourList(int maxEdgeNeigh);

    RCNeighbourList() {}

    /// Destructor
    virtual ~RCNeighbourList();

    /// Sets flag of \ref rclist[i] = true
    void setCoarsePatch(int i)
    {
      rclist[i]->flag = true;
    }

    /// Sets flag of \ref rclist[i] = f
    void setCoarsePatch(int i, bool f)
    {
      rclist[i]->flag = f;
    }

    /// Returns \ref rclist[i].flag
    bool isPatchCoarse(int i) const
    {
      return rclist[i]->flag;
    }

    /// If \ref rclist[i].neigh[j] is not a NULL pointer
    /// \ref rclist[i].neigh[j]->ith will be returned. Otherwise the return value is -1
    int getNeighbourNr(int i, int j) const
    {
      return rclist[i]->neigh[j] ? rclist[i]->neigh[j]->ith : -1;
    }

    /// If \ref rclist[i].neigh[j] is not a NULL pointer
    /// \ref rclist[i].neigh[j]->el will be returned. Otherwise the return value
    /// is NULL
    Element* getNeighbourElement(int i, int j) const
    {
      return rclist[i]->neigh[j] ? rclist[i]->neigh[j]->el : NULL;
    }

    /// Returns \ref rclist[i].el
    Element* getElement(int i) const
    {
      if (static_cast<int>(rclist.size()) <= i)
        return NULL;

      return rclist[i]->el;
    }

    /// Returns number of elements in list
    int getSize() const
    {
      return int(rclist.size());
    }

    /// Sets \ref rclist[i].el to el and \ref rclist[i].flag to cp.
    const Element* setElement(int i, const Element* el, bool cp = false)
    {
      rclist[i]->el = const_cast<Element*>(el);
      rclist[i]->flag = cp;

      return el;
    }

    /// Returns \ref rclist[i].elType
    int getType(int i) const
    {
      return rclist[i]->elType;
    }

    void setType(int i, int type) const
    {
      rclist[i]->elType = type;
    }

    /// If patch can be coarsend return true, else false and reset the
    /// element marks.
    virtual bool doCoarsePatch(int n_neigh);

    /// Sets \ref rclist[i].oppVertex[j] = k
    void setOppVertex(int i, int j, int k)
    {
      rclist[i]->oppVertex[j] = k;
    }

    /// Returns \ref rclist[i].oppVertex[j]
    int getOppVertex(int i, int j)
    {
      return rclist[i]->oppVertex[j];
    }

    /// Sets \ref rclist[i].elType = t
    void setElType(int i, unsigned char t)
    {
      rclist[i]->elType = t;
    }

    /// Sets \ref coarseningManager = cm
    void setCoarseningManager(CoarseningManager* cm)
    {
      coarseningManager = cm;
    }

    /// Fills \ref rclist[i].neigh and \ref rclist[i].oppVertex  infos (0 <= i
    /// < n_neigh)
    void fillNeighbourRelations(int n_neigh, int bound);

    /// Adds those dof's on the parent that are handed on by the
    /// children and adds the dof in the midpoint of the coarsening edge (3d)
    void addDOFParent(int elIndex, DegreeOfFreedom* dof);

    /// If DOFs for higher order are been removed on parents during refinement
    /// they are now added again (2d)
    void addDOFParents(int n_neigh);

    /// Removes DOFs during refinement (3d)
    void removeDOFParent(int index);

    /// Removes DOFs during refinement (2d)
    void removeDOFParents(int n_neigh);

    ///
    void periodicSplit(DegreeOfFreedom* edge[2],
                       DegreeOfFreedom* nextEdge[2],
                       int* n_neigh,
                       int* n_neigh_periodic,
                       RCNeighbourList& periodicList);

    ///
    void clearList();

  protected:
    /// Information about one Element of the patch
    class RCListElement
    {
    public:
      /// Pointer to the Element
      Element* el;

      /// This is the i-th entry in the RCNeighbourList
      int ith;

      /// Only used in the coarsening module: flag is true if the coarsening
      /// edge of the Element is the coarsening edge of the patch, otherwise
      /// flag is false;
      bool flag;

      /// neigh[0/1] neighbour of element to the right/left in the orientation
      /// of the edge, or a NULL pointer in the case of a boundary face (only 3d)
      RCListElement* neigh[2];

      /// opp vertex[0/1] the opposite vertex of neigh[0/1] (only 3d)
      int oppVertex[2];

      /// The element type; is set during looping around the
      /// refinement/coarsening edge; if neighbour information is produced by the
      /// traversal routines, information about the type of an element can not be
      /// accessed via el->getType() and thus has to be stored in the
      /// RCListElement vector (only 3d)
      unsigned char elType;
    };

    /// Refinement/coarsening patch
    std::vector<RCListElement*> rclist;

    /// Pointer to the CoarseningManager
    CoarseningManager* coarseningManager;
  };

} // end namespace
