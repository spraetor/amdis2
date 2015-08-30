/** \file FiniteElemSpace.h */

/** \defgroup FEMSpace FEMSpace
 * @{ <img src="femspace.png"> @}
 */

#pragma once

#include <string>
#include <vector>

#include "AMDiS_fwd.h"

namespace AMDiS
{

  /** \ingroup FEMSpace
   * \brief
   * A FiniteElemSpace is a triple of a DOFAdmin and a set of BasisFunction on a mesh.
   */
  class FiniteElemSpace
  {
  public:
    /// Create an empty fe space.
    FiniteElemSpace() {}

    /// Destructor.
    ~FiniteElemSpace() {}

    ///
    static FiniteElemSpace* provideFeSpace(DOFAdmin* admin,
                                           const BasisFunction* basFcts,
                                           Mesh* mesh,
                                           std::string name = "");

    static void destroyFeSpaces();

#if DEBUG
    /// For debugging it may be useful to get some FE space for a given mesh at a
    /// position in code where it is not possible to access the FE space directly. The
    /// function assumes that there is only one FE space defined for the mesh.

    static FiniteElemSpace* provideFeSpace(Mesh* mesh);
#endif

    FiniteElemSpace& operator=(const FiniteElemSpace& feSpace);

    /// Returns \ref name.
    std::string getName() const
    {
      return name;
    }

    /// Returns \ref admin.
    DOFAdmin* getAdmin() const
    {
      return admin;
    }

    /// Set a new DOF admin.
    void setAdmin(DOFAdmin* a)
    {
      admin = a;
    }

    /// Returns \ref basFcts
    const BasisFunction* getBasisFcts() const
    {
      return basFcts;
    }

    /// Returns \ref mesh
    Mesh* getMesh() const
    {
      return mesh;
    }

    int calcMemoryUsage();

    static void clear();

    /// Returns for a set of FE spaces that FE space having basis functions with
    /// the highest degree.
    static const FiniteElemSpace*
    getHighest(std::vector<const FiniteElemSpace*>& feSpaces);

  protected:
    /// Constructs a FiniteElemSpace with name name_ and the given DOFAdmin,
    /// BasisFunction and Mesh.
    FiniteElemSpace(DOFAdmin* admin,
                    const BasisFunction* basisFcts,
                    Mesh* mesh,
                    std::string name = "");

  protected:
    /// Name of this FiniteElemSpace.
    const std::string name;

    /// DOFAdmin corresponding to this FiniteElemSpace.
    DOFAdmin* admin;

    /// Set of BasisFunction of this FiniteElemSpace.
    const BasisFunction* basFcts;

    /// The Mesh this FiniteElemSpace belongs to.
    Mesh* mesh;

    ///
    static std::vector<FiniteElemSpace*> feSpaces;
  };

} // end namespace AMDiS
