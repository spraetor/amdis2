/** \file DOFContainer.h */

#pragma once

#include "Global.h"

namespace AMDiS 
{
  /** \ingroup DOFAdministration
   * \brief
   * DOFContainer is the base class for objects that stores DOF indices.
   * After a DOFContainer object is registered to a DOFAdmin, the DOF
   * indices of the container will be managed during DOF compression. The
   * DOFAdmin then calls the compress method of every registered DOFContainer.
   */
  class DOFContainer
  {
  public:
    virtual ~DOFContainer() {}

    /** \brief
     * Returns the DOF index at position i. Must be overriden by a concrete
     * DOFContainer.
     */
    virtual DegreeOfFreedom& operator[](DegreeOfFreedom i) = 0;

    virtual void freeDofIndex(DegreeOfFreedom dof) {}

    /** \brief
     * Used by DOFAdmin to actualize the DOF indices in this container after
     * DOF compression.
     */
    virtual void compressDofContainer(int size, std::vector<DegreeOfFreedom> &newDOF)
    {
      FUNCNAME_DBG("DOFContainer::compressDofContainer()");

      for (DegreeOfFreedom i = 0; i < size; i++) {
	DegreeOfFreedom j = newDOF[operator[](i)];

	TEST_EXIT_DBG(j >= 0)
	  ("Invalid DOF %d in DOF container! (%d %d)\n", j, i, operator[](i));

	operator[](i) = j;
      }
    }

  };
  
} // end namespace AMDiS
