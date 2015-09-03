/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors:
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 *
 ******************************************************************************/



/** \file CompositeFEMMethods.h */

#ifndef AMDIS_COMPOSITEFEMMETHODS_H
#define AMDIS_COMPOSITEFEMMETHODS_H

#include "AbstractFunction.h"
#include "DOFVector.h"
#include "FiniteElemSpace.h"
#include "ElementLevelSet.h"

namespace compositeFEM
{

  using namespace AMDiS;

  class CompositeFEMMethods
  {
  public:
    /// Set all dof-values on domain with positive level set function values to val.
    static void setPosLsToVal(
      DOFVector<double>* dof,
      const double& val,
      const DOFVector<double>* lsFct_dof);

    /**
     * Set all dof-values on domain with positive level set function values
     * to values of function fct.
     */
    static void setPosLsToFct(
      DOFVector<double>* dof,
      const AbstractFunction<double, WorldVector<double>>* fct,
      const DOFVector<double>* lsFct_dof);

    /**
     * Print coordinates of all boundary elements to file. Name of file is
     * read from init file.
     */
    static void printBoundaryElements(const std::string fn_str,
                                      ElementLevelSet* elLS,
                                      FiniteElemSpace* feSpace);
  };

}

using compositeFEM::CompositeFEMMethods;

#endif  // AMDIS_COMPOSITEFEMMETHODS_H
