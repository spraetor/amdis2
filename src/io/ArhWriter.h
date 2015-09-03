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



/** \file ArhWriter.h */

#ifndef AMDIS_ARH_WRITER_H
#define AMDIS_ARH_WRITER_H

#include "DOFVector.h"
#include "SystemVector.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "detail/ArhWriter.h"

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

#endif
