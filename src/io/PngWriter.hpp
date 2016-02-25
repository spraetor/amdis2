#pragma once

#ifdef HAVE_PNG

#include "DataCollector.hpp"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Output
     *
     * \brief
     * class which writes a container to a png-file.
     */
    class PngWriter
    {
    public:
      PngWriter(DataCollector<>* dc)
        : dataCollector(dc)
      {
        FUNCNAME("PngWriter::PngWriter()");

        TEST_EXIT(dataCollector->getFeSpace()->getBasisFcts()->getDegree() == 1)
        ("PNG Image output only for linear lagrange elements!\n");

        TEST_EXIT(dataCollector->getMesh()->getDim() == 2)
        ("PNG image output only for 2D!\n");
      }

      /// Writes a PNG image file.
      /// \param filename name of the file to write
      /// \param imageType 0..grayscale, 1..color
      int writeFile(std::string filename, int imageType);

    private:
      /// Datacollector with values for the output file.
      DataCollector<>* dataCollector;
    };

  }
} // end namespace io, AMDiS

#endif
