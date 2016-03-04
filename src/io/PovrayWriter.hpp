// Marcel Schiffel, 23.06.09

#pragma once

#include <string>
#include <fstream>
#include <iostream>

#include "io/DataCollector.hpp"
#include "io/FileWriter.hpp"
#include "Global.hpp"

namespace AMDiS
{
  namespace io
  {

    typedef struct
    {
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
      PovrayWriter(DataCollector<>* dc) : dataCollector(dc), bBox(NULL) { }

      ~PovrayWriter();

      /// writes a povray script for the current time step to the specified file.
      void writeFile(std::string filename);

      // provides the bounding box of the mesh (lazy evaluation, bounding box is computed only once)
      BoundingBox* getBoundingBox(std::ofstream&);

    private:
      DataCollector<>* dataCollector;
      BoundingBox* bBox;

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

  }
} // end namespace io, AMDiS
