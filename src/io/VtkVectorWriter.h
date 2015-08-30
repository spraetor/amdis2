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



/** \file VtkVectorWriter.h */

#ifndef AMDIS_VTKVECTORWRITER_H
#define AMDIS_VTKVECTORWRITER_H

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

#include "DataCollector.h"
#include "SystemVector.h"

#include "detail/VtkVectorWriter.h"
#include "detail/VtkWriter.h"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Output
     * \brief Writer for the ParaView VTU-format for vector data
     *
     * A collection of methods to write various container types to
     * VTU-files. The files generated are text-based XML-structured
     * files and can be read by ParaView. This writer reads vector containers
     * like \ref DOFVector <WorldVector> and generates vector components
     * in the VTU-file.
     **/
    namespace VtkVectorWriter
    {

      /** \brief Write a vector of \ref DataCollector to file.
       *
       * \param dcList vector of \ref DataCollector of arbitrary type T
       * \param filename name of the ouput file
       * \param writeParallel if true: generate container file for parallel VTU files
       * \param writeAs3dVector generate 3 components for the output to allow ParaView to display the data as 3d vector.
       **/
      template<typename T>
      void writeFile(std::vector<DataCollector<T>*>& dcList,
                     std::string filename,
                     bool writeParallel = true,
                     bool writeAs3dVector = false)
      {
        ::AMDiS::io::VtkVectorWriter::Aux<T> writer(&dcList, writeAs3dVector);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
        if (writeParallel)
        {
          int sPos = filename.find(".vtu");
          TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
          std::string name = filename.substr(0, sPos);

          if (MPI::COMM_WORLD.Get_rank() == 0)
            writer.writeParallelFile(name + ".pvtu", MPI::COMM_WORLD.Get_size(),
                                     name, ".vtu");

          filename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.vtu";
        }
#endif
        writer.writeFile(filename);
      }

      /// Write a \ref DOFVector to file. Using container pointer.
      template<typename T>
      void writeFile(DOFVector<T>* values,
                     std::string filename,
                     bool writeParallel = true,
                     bool writeAs3dVector = false)
      {
        std::vector<DataCollector<T>*> dcList(0);
        dcList.push_back(new DataCollector<T>(values->getFeSpace(), values));
        writeFile(dcList, filename, writeParallel, writeAs3dVector);
        delete dcList[0];
      }

      /// Write a \ref DOFVector to file. Using container reference.
      template<typename T>
      void writeFile(DOFVector<T>& values,
                     std::string filename,
                     bool writeParallel = true,
                     bool writeAs3dVector = false)
      {
        writeFile(&values, filename, writeParallel, writeAs3dVector);
      }

      /// Write a vector of \ref DOFVector to file.
      void writeFile(std::vector<DOFVector<double>*>& values,
                     std::string filename,
                     bool writeParallel = true,
                     bool writeAs3dVector = false);

      /// Write a \ref WorldVector of \ref DOFVector to file.
      void writeFile(WorldVector<DOFVector<double>*>& values,
                     std::string filename,
                     bool writeParallel = true,
                     bool writeAs3dVector = false);

      /// Write a \ref SystemVector to file. Using container pointer.
      void writeFile(SystemVector* values,
                     std::string filename,
                     bool writeParallel = true,
                     bool writeAs3dVector = false);

      /// Write a \ref SystemVector to file. Using container reference.
      inline
      void writeFile(SystemVector& values,
                     std::string filename,
                     bool writeParallel = true,
                     bool writeAs3dVector = false)
      {
        writeFile(&values, filename, writeParallel, writeAs3dVector);
      }

    } // end namespace VtkVectorWriter
  }
} // end namespace io, AMDiS

#endif // AMDIS_VTKVECTORWRITER_H
