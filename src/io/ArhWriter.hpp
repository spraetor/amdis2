#pragma once

ERROR

#include "DOFVector.hpp"
#include "SystemVector.hpp"
#include "Mesh.hpp"
#include "MeshStructure.hpp"
#include "io/detail/ArhWriter.hpp"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Output
      * \brief Writer for the AMDiS ARH-format - version 1
      *
      * A collection of methods to write various container types to
      * ARH-files.
      **/
    namespace ArhWriter
    {

      /// Writes a \ref DOFVector to file. Using container pointer.
      inline void writeFile(DOFVector<double>* vec0,
                            std::string filename)
      {
        detail::write(filename, vec0->getFeSpace()->getMesh(), vec0);
      }


      /// Writes a \ref DOFVector to file. Using container reference.
      inline void writeFile(DOFVector<double>& vec0,
                            std::string filename)
      {
        detail::write(filename, vec0.getFeSpace()->getMesh(), &vec0);
      }


      /// Writes a \ref SystemVector to file. Using container pointer.
      inline void writeFile(SystemVector* vec,
                            std::string filename)
      {
        std::vector<DOFVector<double>*> vecs;
        for (int i = 0; i < vec->getSize(); i++)
          vecs.push_back(vec->getDOFVector(i));
        detail::write(filename, vecs[0]->getFeSpace()->getMesh(), vecs);
      }


      /// Writes a vector of \ref DOFVector to file.
      inline void writeFile(std::vector<DOFVector<double>*> vecs,
                            std::string filename)
      {
        detail::write(filename, vecs[0]->getFeSpace()->getMesh(), vecs);
      }

      // ________ below are obsolete functions, for backward compatibility _______


      inline void write(std::string filename, Mesh* mesh,
                        DOFVector<double>* vec0 = NULL,
                        DOFVector<double>* vec1 = NULL,
                        DOFVector<double>* vec2 = NULL)
      {
        FUNCNAME("ArhWriter::write()");
        WARNING("this function is obsolete.\n");
        detail::write(filename, mesh, vec0, vec1, vec2);
      }

      inline void write(std::string filename, Mesh* mesh,
                        std::vector<DOFVector<double>*> vecs,
                        bool writeParallel = true)
      {
        FUNCNAME("ArhWriter::write()");
        WARNING("this function is obsolete.\n");
        detail::write(filename, mesh, vecs, writeParallel);
      }

      inline void writeMacroElement(std::ofstream& file,
                                    MeshStructure& code,
                                    std::vector<std::vector<double>>& values,
                                    uint32_t elIndex)
      {
        FUNCNAME("ArhWriter::writeMacroElement()");
        WARNING("this function is obsolete.\n");
        detail::writeMacroElement(file, code, values, elIndex);
      }
    }// end namespace ArhWriter
  }
} // end namespace io, AMDiS
