#pragma once

#include "io/DataCollector.hpp"
#include "DOFVector.hpp"
#include "FixVec.hpp"
#include "SystemVector.hpp"
#include "io/detail/VtkWriter.hpp"

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
      void writeFile(Container& /*vec*/, std::string /*filename*/,
                     Vtuformat /*format*/ = ASCII,
                     bool /*highPrecision*/ = false,
                     bool /*writeParallel*/ = true)
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
