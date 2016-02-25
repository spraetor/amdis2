#pragma once

#include <list>
#include <vector>

#include "AMDiS_fwd.hpp"
#include "Boundary.hpp"
#include "DataCollector.hpp"
#include "ElementInfo.hpp"
#include "FixVec.hpp"
#include "Flag.hpp"
#include "Mesh.hpp"
#include "Projection.hpp"
#include "VertexInfo.hpp"

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
