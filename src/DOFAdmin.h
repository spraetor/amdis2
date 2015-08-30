/** \file DOFAdmin.h */

/** \defgroup DOFAdministration DOF adaministration module
 * @{ <img src="dof.png"> @}
 * \brief
 * Contains all classes used for the DOF administration.
 */

#pragma once

#include <vector>
#include <memory>
#include <list>

#include "AMDiS_fwd.h"
#include "AMDiS_base.h"
#include "Global.h"
#include "MatrixVector.h"

namespace AMDiS
{
  /** \ingroup DOFAdministration
   * \brief
   * Holds all data about one set of DOFs. It includes information about used and
   * unused DOF indices, as well as lists of DOFIndexed objects and DOFContainer
   * objects, that are automatically resized and resorted during mesh changes.
   */
  class DOFAdmin
  {
  public:
    DOFAdmin();

    /// Constructor
    DOFAdmin(Mesh* m);

    /// Constructor
    DOFAdmin(Mesh* m, std::string aName);

    /// Copy constructor
    DOFAdmin(const DOFAdmin&);

    /// Destructor
    ~DOFAdmin();

    /** \brief
     * Enlarges the number of DOFs that can be managed at least to minsize by
     * a step size of \ref sizeIncrement.
     */
    void enlargeDofLists(int minsize = 0);

    /// assignment operator
    DOFAdmin& operator=(const DOFAdmin&);

    /// Compares two DOFAdmins by their names.
    bool operator==(const DOFAdmin&) const;

    /// Compares two DOFAdmins by their names.
    bool operator!=(const DOFAdmin& ad) const
    {
      return !(ad == *this);
    }

    /** \brief
     * Adds a DOFIndexedBase object to the DOFAdmin. This object will be
     * managed by DOFAdmin from now on.
     */
    void addDOFIndexed(DOFIndexedBase* dofIndexed);

    /// Adds a DOFContainer object to the DOFAdmin.
    void addDOFContainer(DOFContainer* dofContainer);

    /// Removes the given DOFIndexedBase object from DOFAdmin.
    void removeDOFIndexed(DOFIndexedBase* dofIndexed);

    /// Removes the given DOFContainer object from DOFAdmin.
    void removeDOFContainer(DOFContainer* dofContainer);

    /** \brief
     * Removes all holes of unused DOF indices by compressing the used range of
     * indices (it does not resize the vectors). While the global index of a DOF
     * may change, the relative order of DOF indices remains unchanged during
     * compression. This method is usually called after mesh adaption involving
     * higher order elements or coarsening.
     */
    void compress(std::vector<DegreeOfFreedom>& newDofIndex);

    /// Returns an iterator to the begin of \ref dofIndexedList
    std::list<DOFIndexedBase*>::iterator beginDOFIndexed()
    {
      return dofIndexedList.begin();
    }

    /// Returns an iterator to the end of \ref dofIndexedList
    std::list<DOFIndexedBase*>::iterator endDOFIndexed()
    {
      return dofIndexedList.end();
    }

    /** \name getting methods
     * \{
     */

    /// Returns \ref sizeUsed.
    DofIndex::size_type getUsedSize() const
    {
      return sizeUsed;
    }

    /// Returns \ref size
    DofIndex::size_type getSize() const
    {
      return size;
    }

    /// Returns \ref usedCount
    DofIndex::size_type getUsedDofs() const
    {
      return usedCount;
    }

    /// Returns \ref holeCount
    DofIndex::size_type getHoleCount() const
    {
      return holeCount;
    }

    /// Returns \ref name
    std::string getName() const
    {
      return name;
    }

    /// Returns \ref nDof[i], i.e., the number of DOFs for the position i.
    int getNumberOfDofs(int i) const
    {
      return nDof[i];
    }

    /// Returns \ref nDof
    const DimVec<int>& getNumberOfDofs() const
    {
      return nDof;
    }

    /// Returns \ref nPreDof[i]
    int getNumberOfPreDofs(int i) const
    {
      return nPreDof[i];
    }

    /// Returns \ref nPreDof
    const DimVec<int>& getNumberOfPreDofs() const
    {
      return nPreDof;
    }

    /// Returns \ref mesh
    const Mesh* getMesh() const
    {
      return mesh;
    }

    /// Returns \ref dofFree, the array denoting DOFs to be either free or used.
    const std::vector<bool>& getDofFree() const
    {
      return dofFree;
    }

    std::vector<bool>& getDofFree()
    {
      return dofFree;
    }

    /// Returns if the given DOF is free.
    bool isDofFree(DofIndex::size_type i) const
    {
      return dofFree[i];
    }

    /// Sets a DOF to be free or not.
    void setDofFree(DofIndex::size_type i, bool b)
    {
      dofFree[i] = b;
    }

    /// Sets \ref usedSize.
    void setUsedSize(DofIndex::size_type i)
    {
      sizeUsed = i;
    }

    /// Sets \ref usedCount.
    void setUsedCount(DofIndex::size_type i)
    {
      usedCount = i;
    }

    /// Sets \ref firstHole
    void setFirstHole(DofIndex::size_type i)
    {
      TEST_EXIT_DBG(dofFree[i])("There is no hole!\n");

      firstHole = i;
    }

    /** \} */

    /** \name setting methods
     * \{
     */

    /// Sets \ref nDof[i] = v
    void setNumberOfDofs(int i, int v);

    /// Sets all values of \ref nDof
    void setNumberOfDofs(DimVec<int> v)
    {
      nDof = v;
    }

    /// Sets \ref nPreDof[i] = v
    void setNumberOfPreDofs(int i, int v);

    /// Sets \ref name = n
    void setName(std::string n)
    {
      name = n;
    }

    /// Sets \ref mesh = m
    void setMesh(Mesh* m)
    {
      mesh = m;
    }

    int calcMemoryUsage() const;

    void reset()
    {
      init();
    }

    /** \} */

    /** \brief
     * Adds one index to all DOF lists. Used by Mesh::getDof() to provide
     * DOFS for a specific position
     */
    DegreeOfFreedom getDOFIndex(MeshAccessor const&);

    /// Frees index DOF. Used by Mesh::getDof()
    void freeDofIndex(MeshAccessor const&, DegreeOfFreedom dof);

  protected:
    /// Initializes this DOFAdmin
    void init();

  protected:
    /// name of this DOFAdmin
    std::string name;

    /// the mesh this DOFAdmin belongs to
    Mesh* mesh;

    /// The dofFree vector stores for each index whether it is used or not
    std::vector<bool> dofFree;

    /// index of the first free value in \ref dofFree
    DegreeOfFreedom firstHole;

    /// allocated size of managed vectors and matrices
    DofIndex::size_type size;

    /// number of used DOF indices
    DofIndex::size_type usedCount;

    /// number of FREED DOF indices (NOT size - sizeUsed)
    DofIndex::size_type holeCount;

    /// > max. index of a used entry
    DofIndex::size_type sizeUsed;

    /** \brief
     * Number of dofs for each position, i.e., vertex, edge, ..., center,
     * for this DOFAdmin.
     */
    DimVec<int> nDof;

    /// DOFs from previous DOFAdmins
    DimVec<int> nPreDof;

    /// List of all managed DOFIndexed objects.
    std::list<DOFIndexedBase*> dofIndexedList;

    /// List of all managed DOFContainer objects
    std::list<DOFContainer*> dofContainerList;

    /// Size increment used by \ref enlargeDOFLists.
    static const int sizeIncrement;
  };

} // end namespace AMDiS
