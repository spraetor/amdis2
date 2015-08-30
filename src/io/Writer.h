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


/** \file Writer.h */

#ifndef AMDIS_WRITER_H
#define AMDIS_WRITER_H

#include <cstring>
#include "DOFVector.h"
#include "SystemVector.h"

#include "Arh3Writer.h"
#include "DataCollector.h"
#include "DofWriter.h"
#include "GNUPlotWriter.h"
#include "MacroWriter.h"

#ifdef HAVE_PNG
#include "PngWriter.h"
#endif

#include "PovrayWriter.h"
#include "ValueWriter.h"
#include "VtkWriter.h"

#include "detail/ReaderWriter.h"


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

#endif // AMDIS_WRITER_H
