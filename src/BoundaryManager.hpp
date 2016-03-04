#pragma once

// std c++ headers
#include <map>
#include <vector>

// AMDiS includes
#include "AMDiS_fwd.hpp"
#include "AMDiS_base.hpp"
#include "Boundary.hpp"

#include "traits/not_null.hpp"

namespace AMDiS
{
  // forward declaration
  class BoundaryCondition;

  using BoundaryIndexMap = std::map<BoundaryType, BoundaryCondition*>;

  /**
   * \ingroup Assembler
   *
   * \brief
   * A BoundaryManager handles a set of boundary conditions and applies
   * this conditions to DOFVectorBase and DOFMatrix objects. Each DOFVectorBase
   * and each DOFMatrix has its own BoundaryManager.
   */
  class BoundaryManager
  {
  public:
    /// Constructor. Reads the nr. of basis functions and allocates an
    /// array \ref localBound of this size.
    BoundaryManager(not_null<const FiniteElemSpace*> feSpace);

    /// Copy constructor.
    BoundaryManager(BoundaryManager const& bm);

    /// Destructor. Frees the \ref localBound array
    virtual ~BoundaryManager();

    /// Adds a local boundary condition to the list of managed conditions.
    void addBoundaryCondition(BoundaryCondition* localBC);

    void initMatrix(DOFMatrix* matrix);

    void exitMatrix(DOFMatrix* matrix);

    void initVector(DOFVectorBase<double>* vector);

    void exitVector(DOFVectorBase<double>* vector);

    /// Calls DOFVectorBase::fillBoundaryCondition() for each local boundary
    /// condition in \ref localBCs.
    void fillBoundaryConditions(ElInfo* elInfo, DOFVectorBase<double>* vec);

    /// Calls DOFMatrix::fillBoundaryCondition() for each local boundary condition
    /// in \ref localBCs.
    void fillBoundaryConditions(ElInfo* elInfo, DOFMatrix* mat);

    /// Calls BoundaryCondition::boundResidual() for each boundary condition in
    /// \ref localBCs.
    double boundResidual(ElInfo* elInfo, DOFMatrix* matrix,
                         const DOFVectorBase<double>* dv);

    BoundaryCondition* getBoundaryCondition(BoundaryType type)
    {
      return localBCs[type];
    }

    const BoundaryIndexMap& getBoundaryConditionMap() const
    {
      return localBCs;
    }

    void setBoundaryConditionMap(const BoundaryIndexMap& bcs)
    {
      localBCs = bcs;
    }

    /// Returns true, if there is at least one boundary object with the given
    /// boundary id, which implements a periodic boundary.
    static bool isBoundaryPeriodic(BoundaryType b);

  protected:
    /// Map of managed local boundary conditions.
    BoundaryIndexMap localBCs;

    /// Temporary variable for functions fillBoundaryconditions.
    BoundaryType* localBound;

    /// Temporary variable for functions fillBoundaryconditions.
    std::vector<DegreeOfFreedom> dofVec;

    /// Stores the number of byte that were allocated in the constructor for
    /// each localBounds value. Is used to free the memory in the destructor.
    int allocatedMemoryLocalBounds;

    /** \brief
     * For every boundary id we store here all possible boundary object (although
     * it's not clear if it is meaningful to have different boundary conditions on
     * the same boundary id).
     *
     * We have to use this global variable, because the mesh traverse interface
     * does not provide more information about traversed boundaries at elements
     * than the boundary id.
     *
     * TODO: Change interface such that mesh traverse returns the boundary objects
     * directly and we can remove this global variable. The biggest problem will be
     * than serialization and deserialization of the mesh.
     */
    static std::map<BoundaryType, std::vector<BoundaryCondition*>> globalBoundaryMap;
  };

} // end namespace AMDiS
