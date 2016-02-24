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



/** \file DofWriter.h */

#ifndef AMDIS_DOF_WRITER_H
#define AMDIS_DOF_WRITER_H

#include <fstream>
#include <vector>
#include "AMDiS_fwd.h"
#include "FiniteElemSpace.h"
#include "Mesh.h"
#include "SystemVector.h"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Output
     * \brief Writer for the AMDiS DOF-format
     *
     * A collection of methods to write various container types to
     * DOF-files.
     **/
    namespace DofWriter
    {

      /// Interface for general containers not implemented. Specializations below.
      template<typename Container>
      void writeFile(Container& vec, std::string filename)
      {
        ERROR_EXIT("DofWriter not implemented for this container type!\n");
      }

      /// Writes a vector of \ref DOFVector to file.
      void writeFile(std::vector<DOFVector<double>*>& vec,
                     std::string filename);


      /// Writes a \ref DOFVector to file. Using container pointer.
      inline
      void writeFile(DOFVector<double>* dofVector,
                     std::string filename)
      {
        std::vector<DOFVector<double>*> dofVectors;
        dofVectors.push_back(dofVector);
        writeFile(dofVectors, filename);
      }


      /// Writes a \ref DOFVector to file. Using container reference.
      inline
      void writeFile(DOFVector<double>& dofVector,
                     std::string filename)
      {
        writeFile(&dofVector, filename);
      }


      /// Writes a \ref SystemVector to file. Using container pointer.
      inline
      void writeFile(SystemVector* vecs,
                     std::string filename)
      {
        std::vector<DOFVector<double>*> dofVectors;
        for (int i = 0; i < vecs->getSize(); i++)
          dofVectors.push_back(vecs->getDOFVector(i));
        writeFile(dofVectors, filename);
      }


      /// Writes a \ref SystemVector to file. Using container reference.
      inline
      void writeFile(SystemVector& vecs,
                     std::string filename)
      {
        writeFile(&vecs, filename);
      }

    } // end namespace DofWriter
  } // end namespace io
} // end AMDiS

#endif
