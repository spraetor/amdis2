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


#include <fstream>

#include "MacroWriter.h"
#include "DataCollector.h"
#include "Mesh.h"
#include "DOFAdmin.h"
#include "FiniteElemSpace.h"
#include "DOFVector.h"
#include "DOFIterator.h"
#include "ElInfo.h"
#include "SurfaceRegion_ED.h"
#include "ElementRegion_ED.h"
#include <string>
#include "Traverse.h"

namespace AMDiS
{
  namespace io
  {

    int MacroWriter::writeMacro(DataCollector<>* dc,
                                std::string name,
                                double time,
                                int level,
                                Flag traverseFlag,
                                bool (*writeElem)(ElInfo*))
    {
      FUNCNAME("MacroWroter::writeFile()");

      TEST_EXIT(dc)("no data collector\n");

      std::ofstream file;
      std::list<ElementInfo>* elements = dc->getElementInfos();
      DOFVector<std::list<VertexInfo>>* vertexInfos = dc->getVertexInfos();

      int dow = Global::getGeo(WORLD);
      int dim = dc->getMesh()->getDim();
      int nv = dc->getNumberVertices();
      int ne = dc->getNumberElements();
      int vertices = dc->getMesh()->getGeo(VERTEX);

      file.open(name.c_str());

      // === print file header ===
      file << "DIM: " << dim << std::endl;
      file << "DIM_OF_WORLD: " << dow << std::endl << std::endl;

      file << "number of vertices: " << nv << std::endl;
      file << "number of elements: " << ne << std::endl << std::endl;

      // === print vertex coords and remember global output indices ===
      file << "vertex coordinates:" << std::endl;

      DOFVector<std::list<VertexInfo>>::Iterator it(vertexInfos, USED_DOFS);
      int counter = 0;

      // for all DOFs
      for (it.reset(); !it.end(); ++it)
      {
        // for all vertex infos of this DOF
        std::list<VertexInfo>::iterator it2;
        for (it2 = it->begin(); it2 != it->end(); ++it2)
        {
          it2->outputIndex = counter++;
          for (int i = 0; i < dow; i++)
            file << std::scientific << it2->coords[i] << " ";
          file << std::endl;
        }
      }

      // === print element vertices ===
      file << std::endl << "element vertices:" << std::endl;

      // iterate the element list
      std::list<ElementInfo>::iterator elementIt;

      for (elementIt = elements->begin(); elementIt != elements->end(); ++elementIt)
      {
        // for all vertices
        for (int i = 0; i < vertices; i++)
          file << elementIt->vertexInfo[i]->outputIndex << " ";
        file << "\n";
      }

      // === print boundaries ===
      file << std::endl << "element boundaries:" << std::endl;

      for (elementIt = elements->begin(); elementIt != elements->end(); ++elementIt)
      {
        // for all vertices
        for (int i = 0; i < vertices; i++)
          file << elementIt->boundary[i] << " ";
        file << "\n";
      }

      // === print neighbours ===
      file << std::endl << "element neighbours:" << std::endl;

      for (elementIt = elements->begin(); elementIt != elements->end(); ++elementIt)
      {
        // for all vertices
        for (int i = 0; i < vertices; i++)
          file << elementIt->neighbour[i] << " ";
        file << "\n";
      }

      // === print boundary projections ===
      file << std::endl << "projections:" << std::endl;

      for (elementIt = elements->begin(); elementIt != elements->end(); ++elementIt)
      {
        // for all vertices
        for (int i = 0; i < vertices; i++)
        {
          if (elementIt->projection[i])
            file << elementIt->projection[i]->getID() << " ";
          else
            file << "0 ";
        }
        file << std::endl;
      }

      // === print element regions ===
      file << std::endl << "element region:" << std::endl;

      for (elementIt = elements->begin(); elementIt != elements->end(); ++elementIt)
        file << elementIt->elementRegion << std::endl;

      // === print surface regions ===
      file << std::endl << "surface region:" << std::endl;

      for (elementIt = elements->begin(); elementIt != elements->end(); ++elementIt)
      {
        for (int i = 0; i < vertices; i++)
          file << elementIt->surfaceRegions[i] << " ";

        file << std::endl;
      }

      // === print element types if necessary ===
      if (dim == 3)
      {
        file << std::endl << "element type:" << std::endl;

        for (elementIt = elements->begin(); elementIt != elements->end(); ++elementIt)
          file << (int)(elementIt->type) << " ";
      }

      file.close();

      return 0;
    }


    void MacroWriter::writePeriodicFile(DataCollector<>* dc, std::string filename)
    {
      FUNCNAME("MacroWriter::writePeriodicFile2");
      TEST_EXIT(dc)("no data collector\n");

      std::ofstream file;

      file.open(filename.c_str());

      file << "associations: " << dc->getNumberConnections() << std::endl;
      file << std::endl << "mode  bc  el1 - local vertices <->  el2 - local vertices"
           << std::endl;

      std::list<PeriodicInfo>::iterator periodicIt;

      // Iterate on all periodic connections
      for (periodicIt = dc->getPeriodicInfos()->begin();
           periodicIt != dc->getPeriodicInfos()->end(); ++periodicIt)
      {

        std::map<int, int>::iterator mapIt;

        // Write mode and type of periodic connection
        file << periodicIt->mode << " " << periodicIt->type << " ";

        // Write index of the first element
        file << periodicIt->outputIndex << " ";

        // Write local indices of the first element
        for (mapIt = periodicIt->vertexMap.begin(); mapIt != periodicIt->vertexMap.end();
             ++mapIt)
          file << mapIt->first << " ";

        // Write index of the second element
        file << periodicIt->neighIndex << " ";

        // Write local indices of the second element
        for (mapIt = periodicIt->vertexMap.begin(); mapIt != periodicIt->vertexMap.end();
             ++mapIt)
          file << mapIt->second << " ";

        file << std::endl;
      }


      file.close();
    }

  }
} // end namespace io, AMDiS
