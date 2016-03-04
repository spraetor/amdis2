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



/** \file MeshManipulation.h */

#ifndef AMDIS_MESHMANIPULATION_H
#define AMDIS_MESHMANIPULATION_H

#include <set>
#include <map>

#include "AMDiS_fwd.h"
#include "parallel/ElementObjectDatabase.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    class MeshManipulation
    {
    public:
      MeshManipulation(Mesh* m);

      ~MeshManipulation();

      void deleteDoubleDofs(std::vector<const FiniteElemSpace*>& feSpaces,
                            std::set<MacroElement*>& newMacroEl,
                            ElementObjectDatabase& elObj);

      void deleteDoubleDofs(std::vector<const FiniteElemSpace*>& feSpaces,
                            std::vector<MacroElement*>& newMacroEl,
                            ElementObjectDatabase& elObjDb);

      /** \brief
       * Starts the procedure to fit a given edge/face of a macro element with a
       * mesh structure code. This functions prepares some data structures and
       * call then \ref fitElementToMeshCode, that mainly refines the element such
       * that it fits to the mesh structure code.
       *
       * \param[in] code      The mesh structure code to which the edge/face of
       *                      an element must be fitted.
       * \param[in] boundEl   Defines the element and its sub-object, i.e., edge or
       *                      face that must be fitted.
       *
       * \return   Returns true, if the mesh was changed, i.e. refined, to fit the
       *           element to the structure code.
       */
      bool fitElementToMeshCode(MeshStructure& code, BoundaryObject& boundEl);

    private:
      /** \brief
       * Recursively fits a given mesh structure code to an edge/face of an element.
       * This function is always initialy called from \ref startFitElementToMeshCode.
       *
       * \param[in] subObj       Defines whether an edge or a face must be fitted.
       * \param[in] ithObj       Defines which edge/face must be fitted.
       *
       * \return   Returns true, if the mesh was changed, i.e. refined, to fit the
       *           element to the structure code.
       */
      bool fitElementToMeshCode(GeoIndex subObj, int ithObj);



    private:
      Mesh* mesh;

      RefinementManager* refineManager;

      MeshStructure* pcode;

      TraverseStack* pstack;

      bool rMode;
    };

  }
}

#endif

