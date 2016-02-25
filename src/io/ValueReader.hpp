#pragma once

#include <cstring>

#include "DOFVector.hpp"
#include "Mesh.hpp"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     *
     * \brief
     * Namespace of methods which read a value file in AMDiS format and copies
     * the data to a DOF vector.
     */
    namespace ValueReader
    {

      template<typename Container>
      void readValue(std::string /*filename*/,
                     Mesh* /*mesh*/,
                     Container& /*vec*/,
                     MacroInfo* /*macroFileInfo*/)
      {
        ERROR_EXIT("ValueReader not implemented for this container type!\n");
      }

      /// Copies the values of a value file to a DOF vector.
      void readValue(std::string filename,
                     Mesh* mesh,
                     DOFVector<double>* dofVector,
                     MacroInfo* macroFileInfo);

      inline void readValue(std::string filename,
                     Mesh* mesh,
                     DOFVector<double>& dofVector,
                     MacroInfo* macroFileInfo)
      {
        readValue(filename, mesh, &dofVector, macroFileInfo);
      }

    } // end namespace ValueReader
  }
} // end namespace io, AMDiS
