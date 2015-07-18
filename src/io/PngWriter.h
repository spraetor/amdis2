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



/** \file PngWriter.h */

#ifndef AMDIS_PNGWRITER_H
#define AMDIS_PNGWRITER_H

#if defined HAVE_PNG

#include "DataCollector.h"

namespace AMDiS { namespace io {

  /** \ingroup Output
   *  
   * \brief
   * class which writes a container to a png-file.
   */ 
  class PngWriter
  {
  public:
    PngWriter(DataCollector<> *dc)
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
  
} } // end namespace io, AMDiS

#endif

#endif
