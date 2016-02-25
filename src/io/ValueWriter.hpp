#pragma once

#include "Global.hpp"
#include "Flag.hpp"
#include "Mesh.hpp"
#include "io/DataCollector.hpp"
#include "AMDiS_fwd.hpp"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Output
     * \brief
     * ValueWriter is a namespace which provied writers to write the values of a DOFVector
     * to an ascii file named values->name.'dat'. This output is done
     * via two leaf-traversals of values->feSpace->mesh. In the first traversal
     * the values at the vertices are printed, in the second these at the
     * interpolation points of each element. For a closer disription of the
     * output format see (...link fehlt noch)
     */
    namespace ValueWriter
    {
      /// Writes DOFVector values to a file.
      void writeValues(DataCollector<>* dc, std::string filename,
                       double time = 0.0, int level = -1,
                       Flag traverseFlag = Mesh::CALL_LEAF_EL,
                       bool (*writeElem)(ElInfo*) = NULL);

    } // end namespace ValueWriter
  }
} // end namespace io, AMDiS
