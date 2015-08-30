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



/** \file VtkWriter.h */

#ifndef AMDIS_VTKWRITER_H
#define AMDIS_VTKWRITER_H

#include "DataCollector.h"
#include "DOFVector.h"
#include "FixVec.h"
#include "SystemVector.h"
#include "detail/VtkWriter.h"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Output
     * \brief Writer for the ParaView VTU-format
     *
     * A collection of methods to write various container types to
     * VTU-files. The files generated are text-based XML-structured
     * files and can be read by ParaView.
     **/

    namespace VtkWriter
    {
      /// Interface for general containers not implemented. Specializations below.
      template<typename Container>
      void writeFile(Container& vec, std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true)
      {
        ERROR_EXIT("VtkWriter not implemented for this container type!\n");
      }

      /// Write a \ref DOFVector to a file. Using a container pointer.
      void writeFile(DOFVector<double>* values,
                     std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true
                    );

      /// Write a \ref DOFVector to file. Using a container reference.
      inline
      void writeFile(DOFVector<double>& values,
                     std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true
                    )
      {
        writeFile(&values, filename, format, highPrecision, writeParallel);
      }

      /// Write a vector of \ref DOFVector to file.
      void writeFile(std::vector<DOFVector<double>*>& values,
                     std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true
                    );

      /// Write a \ref WorldVector of \ref DOFVector to file.
      void writeFile(WorldVector<DOFVector<double>*>& values,
                     std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true
                    );

      /// Write a \ref DOFVector <WorldVector> to file, by converting
      /// to a \ref WorldVector of \ref DOFVector and calling the appropriate
      /// method above. Using a container pointer.
      void writeFile(DOFVector<WorldVector<double>>* values,
                     std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true
                    );

      /// Write a \ref DOFVector <WorldVector> to file, by converting
      /// to a \ref WorldVector of \ref DOFVector and calling the appropriate
      /// method above. Using a container reference.
      inline
      void writeFile(DOFVector<WorldVector<double>>& values,
                     std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true
                    )
      {
        writeFile(&values, filename, format, highPrecision, writeParallel);
      }

      /// Write a \ref SystemVector to file. Using a container pointer.
      void writeFile(SystemVector* values,
                     std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true
                    );

      /// Write a \ref SystemVector to file. Using a container reference.
      inline
      void writeFile(SystemVector& values,
                     std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true
                    )
      {
        writeFile(&values, filename, format, highPrecision, writeParallel);
      }

      /// Write a vector of \ref DataCollector to file.
      void writeFile(std::vector<DataCollector<>*>& dcList,
                     std::string filename,
                     Vtuformat format = ASCII,
                     bool highPrecision = false,
                     bool writeParallel = true
                    );

    } // end namespace VtkWriter
  }
} // end namespace io, AMDiS

#endif // AMDIS_VTKWRITER_H
