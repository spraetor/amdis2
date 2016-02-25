#pragma once

#include <cstring>

#include "DOFVector.hpp"
#include "SystemVector.hpp"

#include "io/Arh3Writer.hpp"
#include "io/DataCollector.hpp"
#include "io/DofWriter.hpp"
#include "io/GNUPlotWriter.hpp"
#include "io/MacroWriter.hpp"

#ifdef HAVE_PNG
#include "io/PngWriter.hpp"
#endif

#include "io/PovrayWriter.hpp"
#include "io/ValueWriter.hpp"
#include "io/VtkWriter.hpp"

#include "io/detail/ReaderWriter.hpp"


namespace AMDiS
{
  namespace io
  {

    /** \ingroup Output
     * \brief Writes a file in a format corresponding to the file-extension
     *
     * Filewriter interface that selects a writer depending on the file
     * extension of the filename.
     *
     * @param container One of DOFVector&, std::vector<DOFVector*>,
     * 			SystemVector&, Mesh&.
     * @param filename Filename of the output file, used to extract extension.
     **/
    template<typename Container>
    void writeFile(Container& container,
                   std::string filename)
    {
      std::string ext = filename.substr(filename.find_last_of("."));

      if (ext == ".1d" || ext == ".2d" || ext == ".3d")
      {
        DataCollector<double> dc(detail::getFeSpace(container, 0),
                                 detail::getDOFVector(container, 0));
        MacroWriter::writeMacro(&dc, filename.c_str());
      }
      else if (ext == ".arh")
      {
        Arh3Writer::writeFile(container, filename);
      }
      else if (ext == ".dat")
      {
        DataCollector<double> dc(detail::getFeSpace(container, 0),
                                 detail::getDOFVector(container, 0));
        ValueWriter::writeValues(&dc, filename);
      }
      else if (ext == ".dof")
      {
        DofWriter::writeFile(container, filename);
      }
      else if (ext == ".plt" || ext == ".gp" || ext == ".gnu" || ext == ".gpi" ||
               ext == ".gnuplot")
      {
        GNUPlotWriter::writeFile(container, filename);
      }
#ifdef HAVE_PNG
      else if (ext == ".png")
      {
        DataCollector<double> dc(detail::getFeSpace(container, 0),
                                 detail::getDOFVector(container, 0));
        PngWriter pngWriter(&dc);
        int color = 1;
        pngWriter.writeFile(filename, color);
      }
#endif
      else if (ext == ".pov")
      {
        DataCollector<double> dc(detail::getFeSpace(container, 0),
                                 detail::getDOFVector(container, 0));
        PovrayWriter povWriter(&dc);
        povWriter.writeFile(filename);
      }
      else if (ext == ".vtu")
      {
        VtkWriter::writeFile(container, filename);
      }
      else
      {
        ERROR_EXIT("File-extensions %s can not be assigned to a writer!\n", ext.c_str());
      }
    }


    /** \ingroup Output
     * \brief Wrapper for pointers to container types
     **/
    template<typename Container>
    void writeFile(Container* container,
                   std::string filename)
    {
      writeFile(*container, filename);
    }

  } // end namespace io
} // end namespace AMDiS
