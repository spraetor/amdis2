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



/** \file MacroWriter.h */

#ifndef AMDIS_MACROWRITER_H
#define AMDIS_MACROWRITER_H

#include <list>
#include <vector>
#include "AMDiS_fwd.h"
#include "VertexInfo.h"
#include "ElementInfo.h"
#include "DataCollector.h"
#include "FixVec.h"
#include "Boundary.h"
#include "Projection.h"
#include "Flag.h"
#include "Mesh.h"

namespace AMDiS
{
  namespace io
  {

    /**
     * \ingroup Output
     *
     * \brief
     * Writes the current leaf elements of a mesh as macro triangulation to
     * a text file. Pure static class.
     */
    struct MacroWriter
    {
      /// Writes the leaf elements of a Mesh as a macro triangulation to a file.
      static int writeMacro(DataCollector<>* dc,
                            std::string name,
                            double time = 0.0,
                            int level = -1,
                            Flag traverseFlag = Mesh::CALL_LEAF_EL,
                            bool (*writeElem)(ElInfo*) = NULL);

      /// Init \ref periodicFile for the next macro to be written.
      static void writePeriodicFile(DataCollector<>* dc, std::string filename);
    };

  }
} // end namespace io, AMDiS

#endif
