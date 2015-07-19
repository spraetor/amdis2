/** \file BoundaryCondition.h */

#pragma once

#include "Boundary.h"
#include "FiniteElemSpace.h"
#include "AMDiS_fwd.h"

namespace AMDiS {

  /**
   * \ingroup Assembler
   *
   * \brief
   * Sub class of BoundaryCondition. Local boundary conditions are filled
   * while mesh traversal.
   */
  class BoundaryCondition
  {
  public:
    /// Constructor.
    BoundaryCondition(BoundaryType type, 
		      const FiniteElemSpace *rowFeSpace_,
		      const FiniteElemSpace *colFeSpace_ = NULL) 
      : boundaryType(type),
      	rowFeSpace(rowFeSpace_),
      	colFeSpace(colFeSpace_)
    {
      if (!colFeSpace) 
        colFeSpace = rowFeSpace;
    }

    /// Returns \ref boundaryType.
    BoundaryType getBoundaryType() const
    { 
      return boundaryType; 
    }

    /// Returns \ref rowFeSpace.
    const FiniteElemSpace *getRowFeSpace() const
    { 
      return rowFeSpace; 
    }

    /// Returns \ref rowFeSpace.
    const FiniteElemSpace *getColFeSpace() const
    { 
      return colFeSpace; 
    }

    virtual void initMatrix(DOFMatrix*) {}

    virtual void exitMatrix(DOFMatrix*) {}

    virtual void initVector(DOFVectorBase<double>*) {}

    virtual void exitVector(DOFVectorBase<double>*) {}

    /// Destructor.
    virtual ~BoundaryCondition() {}

    /// Adds the local boundary condition for elInfo to object.
    /// The dofIndices and localBound as well as nBasFcts are determined by
    // the calling BoundaryManager.
    virtual void fillBoundaryCondition(DOFMatrix             *matrix,
				       ElInfo                *elInfo,
				       const DegreeOfFreedom *dofIndices,
				       const BoundaryType    *localBound,
				       int                    nBasFcts) {}
  
    /// Adds the local boundary condition for elInfo to vector.
    /// The dofIndices and localBound as well as nBasFcts are determined by
    /// the calling BoundaryManager.
    virtual void fillBoundaryCondition(DOFVectorBase<double>     *vector, 
				       ElInfo                *elInfo,
				       const DegreeOfFreedom *dofIndices,
				       const BoundaryType    *localBound,
				       int                    nBasFcts) {}
  
    /// Returns the boundary residual for the given element. Called by estimator.
    virtual double boundResidual(ElInfo *elInfo, 
                        				 DOFMatrix *matrix,
                        				 const DOFVectorBase<double> *dv) 
    { 
      return 0.0; 
    }

    /// Returns whether the condition must be treated as Dirichlet condition
    /// while assemblage.
    virtual bool isDirichlet() const
    { 
      return false; 
    }

    /// Returns whether the boundary condition is a periodic condition or not.
    virtual bool isPeriodic() const
    {
      return false;
    }

    /// In some situations it may be required to set Dirichlet boundary 
    /// conditions, but not to apply them to the matrix. This is for example the
    /// case, if the boundary condition is set to a couple matrix. Then, the
    /// boundary conditions must be applied to the couple matrix, but they are
    /// set to all matrices in this row (to ensure that there are no other
    /// element entries in the Dirichlet boundary condition rows).
    virtual bool applyBoundaryCondition() const
    {
      return true;
    }

  protected:
    /// Speciefies for which parts of the boundary the condition holds.
    /// This id corresponds to the boundary numbers spcified in the
    /// macro file. 
    BoundaryType boundaryType;

    /// FiniteElemSpace for this BoundaryCondition.
    const FiniteElemSpace *rowFeSpace;

    /// FiniteElemSpace for this BoundaryCondition.
    const FiniteElemSpace *colFeSpace;
  };

} // end namespace AMDiS
