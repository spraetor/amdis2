/** \file PeriodicBC.h */

#pragma once

#include <map>
#include <vector>

#include "AMDiS_fwd.h"
#include "BoundaryCondition.h"
#include "FixVec.h"

namespace AMDiS
{

  template <class T>
  class DimVecLess
  {
  public:
    bool operator()(DimVec<T> const& v1, DimVec<T> const& v2) const
    {
      int  size = v1.getSize();
      for (int i = 0; i < size; i++)
      {
        if (v1[i] < v2[i])
          return true;
        if (v1[i] > v2[i])
          return false;
      }
      return false;
    }
  };

  class PeriodicDOFMapping
  {
  public:
    PeriodicDOFMapping(BasisFunction const* basFcts);

    ~PeriodicDOFMapping();

  public:
    static PeriodicDOFMapping* providePeriodicDOFMapping(BasisFunction const* basFcts);

    const DegreeOfFreedom* getDOFPermutation(DimVec<int> const& vertexPermutation);

  protected:
    /// Basis functions the mapping object is defined on.
    const BasisFunction* basFcts;

    std::map<DimVec<int>, DegreeOfFreedom*, DimVecLess<int>> dofPermutation;

    /// Maps from coordinates (in barycentric coordinates) of a local basis
    /// function to the local index of that basis function.
    std::map<DimVec<double>, int, DimVecLess<double>> indexOfCoords;

    /// Global array that provids for each existing basis functions a unique
    /// mapping object.
    static std::vector<PeriodicDOFMapping*> mappings;
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * Sub class of BoundaryCondition. Implements Periodic boundary conditions.
   */
  class PeriodicBC : public BoundaryCondition
  {
  public:
    /// Constructor.
    PeriodicBC(BoundaryType type, FiniteElemSpace const* rowFeSpace);

    ~PeriodicBC();

    virtual void initMatrix(DOFMatrix* matrix) override;

    virtual void fillBoundaryCondition(DOFMatrix* matrix,
                                       ElInfo* elInfo,
                                       DegreeOfFreedom const* dofIndices,
                                       BoundaryType const* localBound,
                                       int nBasFcts) override;

    virtual void exitMatrix(DOFMatrix* matrix) override;

    virtual void exitVector(DOFVectorBase<double>* vector) override;

    /// We are defining periodic boundary conditions, so return always true here.
    virtual bool isPeriodic() const override
    {
      return true;
    }

  protected:
    VertexVector* associated;

    PeriodicDOFMapping* periodicDOFMapping;

    DOFMatrix* masterMatrix;
  };

} // end namespace AMDiS
