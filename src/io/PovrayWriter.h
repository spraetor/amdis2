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



/** \file FileWriter.h */

// Marcel Schiffel, 23.06.09

#ifndef AMDIS_POVRAYWRITER_H
#define AMDIS_POVRAYWRITER_H

#include <string>
#include <fstream>
#include <iostream>
#include "Global.h"
#include "DataCollector.h"
#include "FileWriter.h"

namespace AMDiS { namespace io {

  typedef struct {
    double minx;
    double maxx;
    double miny;
    double maxy;
    double minz;
    double maxz;
  } BoundingBox;


  /** \ingroup Output
   *  
   * \brief
   * Class which writes a container to a povray-file.
   */ 
  class PovrayWriter 
  {
  public:
    PovrayWriter(DataCollector<> *dc) : dataCollector(dc), bBox(NULL) { }

    ~PovrayWriter();
    
    /// writes a povray script for the current time step to the specified file.
    void writeFile(std::string filename);

    // provides the bounding box of the mesh (lazy evaluation, bounding box is computed only once)
    BoundingBox *getBoundingBox(std::ofstream&);

  private:
    DataCollector<> *dataCollector;
    BoundingBox *bBox;

    // TODO: remove!
    void tryMeshTraversal(std::ofstream&);

    void writeTestStuff(std::ofstream&, DataCollector<>&); // TODO: remove/rename
    void writeHeader(std::ofstream&);
    void writeIncludes(std::ofstream&);
    void writeCamera(std::ofstream&);
    void writeLight(std::ofstream&);
    void writeMesh2(std::ofstream&, DataCollector<>&);
    void writeVertexVectors(std::ofstream&, DataCollector<>&);
    void writeTextureList(std::ofstream&, DataCollector<>&);
    void writeFaceIndices(std::ofstream&, DataCollector<>&);
  };

} } // end namespace io, AMDiS

#endif
