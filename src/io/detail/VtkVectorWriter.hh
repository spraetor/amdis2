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

/** \file VtkVectorWriter.hh */

#ifndef AMDIS_VTKVECTORWRITER_DETAIL_HH
#define AMDIS_VTKVECTORWRITER_DETAIL_HH

#include <list>
#include <vector>

#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

#ifdef HAVE_COMPRESSION
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#endif

#include "DOFVector.h"
#include "io/DataCollector.h"
// #include "SurfaceRegion_ED.h"
// #include "ElementRegion_ED.h"

#include "io/detail/VtkWriter.h"

namespace AMDiS
{
  namespace io
  {
    namespace VtkVectorWriter
    {

      template<typename S>
      int Aux<S>::writeFile(std::string name)
      {
        FUNCNAME("Aux<S>::writeFile()");

        boost::iostreams::filtering_ostream file;
#ifdef HAVE_COMPRESSION
        switch (compress)
        {
        case GZIP:
          file.push(boost::iostreams::gzip_compressor());
          name.append(".gz");
          break;
        case BZIP2:
          file.push(boost::iostreams::bzip2_compressor());
          name.append(".bz2");
          break;
        default:
          break;
        }
#endif

        {
          std::ofstream swapfile(name.c_str(), std::ios::out | std::ios::trunc);
          TEST_EXIT(swapfile.is_open())
          ("Cannot open file %s for writing!\n", name.c_str());
          swapfile.close();
        }

        file.push(boost::iostreams::file_descriptor_sink(name, std::ios::trunc));
        writeFileToStream(file);

        return 0;
      }


      template<typename S>
      void Aux<S>::writeParallelFile(std::string name, int nRanks,
                                     std::string fnPrefix, std::string fnPostfix)
      {
        FUNCNAME("Aux<S>::writeParallelFile()");

        boost::iostreams::filtering_ostream file;
        {
          std::ofstream swapfile(name.c_str(), std::ios::out | std::ios::trunc);
          TEST_EXIT(swapfile.is_open())
          ("Cannot open file %s for writing!\n", name.c_str());
          swapfile.close();
        }
        file.push(boost::iostreams::file_descriptor_sink(name, std::ios::trunc));

        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"PUnstructuredGrid\">\n";
        file << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
        file << "    <PPoints>\n"
             << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>\n"
             << "    </PPoints>\n";
        file << "    <PCells>\n"
             << "      <PDataArray type=\"Int32\" Name=\"offsets\"/>\n"
             << "      <PDataArray type=\"UInt8\" Name=\"types\"/>\n"
             << "      <PDataArray type=\"Int32\" Name=\"connectivity\"/>\n"
             << "    </PCells>\n";
        file << "    <PPointData>\n";

        for (int i = 0; i < static_cast<int>(dataCollector->size()); i++)
        {
          S temp = (*((*dataCollector)[i]->getValues()))[0];
          int numComponent = num_rows(temp);
          std::string name = (*dataCollector)[i]->getValues()->getName();
          if (name.find_first_not_of(" \n\r\t") == std::string::npos)
            name = "value" + std::to_string(i);

          file << "      <PDataArray type=\"Float32\" Name=\""
               << name
               << "\" format=\"ascii\""
               << (numComponent > 1 ? " NumberOfComponents=\"" + std::to_string(numComponent) + "\"" : "")
               << "/>\n";
        }

        file << "    </PPointData>\n";

        for (int i = 0; i < nRanks; i++)
        {
          std::stringstream oss;
          oss << fnPrefix << "-p" << i << "-" << fnPostfix;
          boost::filesystem::path filepath(oss.str());
          file << "    <Piece Source=\""
               << boost::filesystem::basename(filepath)
               << boost::filesystem::extension(filepath) << "\"/>\n";

        }

        file << "  </PUnstructuredGrid>\n";
        file << "</VTKFile>\n";
      }


      template<typename S> template<typename T>
      void Aux<S>::writeFileToStream(T& file)
      {
        int nVertices = (*dataCollector)[0]->getNumberVertices();
        int nElements = (*dataCollector)[0]->getNumberElements();
        int vertices = (*dataCollector)[0]->getMesh()->getGeo(VERTEX);

        if ((dim == 2) && (degree == 2))
        {
          nVertices += (*dataCollector)[0]->getNumberInterpPoints();
          nElements *= 4;
        }
        else if ((dim == 2) && (degree == 3))
        {
          nVertices += (*dataCollector)[0]->getNumberInterpPoints();
          nElements *= 9;
        }
        else if ((dim == 2) && (degree == 4))
        {
          nVertices += (*dataCollector)[0]->getNumberInterpPoints();
          nElements *= 16;
        }
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "  <UnstructuredGrid>\n";
        file << "    <Piece NumberOfPoints=\"" << nVertices
             << "\" NumberOfCells=\"" << nElements << "\">\n";
        file << "      <Points>\n";
        file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

        writeVertexCoords(file);

        file << "        </DataArray>\n";
        file << "      </Points>\n";
        file << "      <Cells>\n";
        file << "        <DataArray type=\"Int32\" Name=\"offsets\">\n";

        for (int i = 0; i < nElements; i++)
          file << " " << (i + 1) * vertices << "\n";

        file << "        </DataArray>\n";
        file << "        <DataArray type=\"UInt8\" Name=\"types\">\n";

        for (int i = 0; i < nElements; i++)
        {
          switch (vertices)
          {
          case 2:
            file << " 3\n";
            break;
          case 3:
            file << " 5\n";
            break;
          case 4:
            file << " 10\n";
            break;
          default:
            break;
          }
        }

        file << "        </DataArray>\n";
        file << "        <DataArray type=\"Int32\" Name=\"connectivity\">\n";

        writeConnectivity(file);

        file << "        </DataArray>\n";
        file << "      </Cells>\n";
        file << "      <PointData>\n";

        for (int i = 0; i < static_cast<int>(dataCollector->size()); i++)
        {
          S temp = (*((*dataCollector)[i]->getValues()))[0];
          int numComponent = num_rows(temp);
          std::string name = (*dataCollector)[i]->getValues()->getName();
          if (name.find_first_not_of(" \n\r\t") == std::string::npos)
            name = "value" + std::to_string(i);

          file << "        <DataArray type=\"Float32\" Name=\""
               << name
               << "\" format=\"ascii\"" << (numComponent > 1 ? " NumberOfComponents=\"" + std::to_string(numComponent) + "\"" : "") << ">\n";

          writeVertexValues(file, i);

          file << "        </DataArray>\n";
        }

        file << "      </PointData>\n";
        file << "    </Piece>\n";
        file << "  </UnstructuredGrid>\n";
        file << "</VTKFile>\n";
      } //


      template<typename S> template<typename T>
      void Aux<S>::writeVertexCoords(T& file)
      {
        using ::AMDiS::io::VtkWriter::detail::writeCoord;

        DOFVector<std::list<VertexInfo>>* vertexInfos =
                                        (*dataCollector)[0]->getVertexInfos();
        DOFVector<std::list<VertexInfo>>::Iterator it(vertexInfos, USED_DOFS);
        int counter = 0;

        // For all DOFs of vertices, write the coordinates.
        for (it.reset(); !it.end(); ++it)
        {
          // for all vertex infos of this DOF
          for (std::list<VertexInfo>::iterator it2 = it->begin(); it2 != it->end(); ++it2)
          {
            it2->outputIndex = counter++;
            writeCoord(file, it2->coords);
          }
        }

        // For the second dim case, write also the interpolation points.
        if ((dim == 2) && (degree > 1))
        {
          DOFVector<std::list<WorldVector<double>>>* interpPointCoords = (*dataCollector)[0]->getInterpPointCoords();
          DOFVector<std::list<WorldVector<double>>>::Iterator pointIt(interpPointCoords, USED_DOFS);

          counter = 0;
          for (pointIt.reset(); !pointIt.end(); ++pointIt)
          {
            for (std::list<WorldVector<double>>::iterator it2 = pointIt->begin();
                 it2 != pointIt->end(); ++it2)
            {
              counter++;
              writeCoord(file, *it2);
            }
          }
        }
      }


      template<typename S> template<typename T>
      void Aux<S>::writeVertexValues(T& file, int componentNo)
      {
        DOFVector<int>* interpPointInd = (*dataCollector)[componentNo]->getInterpPointInd();
        DOFVector<S>* values = (*dataCollector)[componentNo]->getValues();
        DOFVector<std::list<WorldVector<double>>>* dofCoords = (*dataCollector)[componentNo]->getDofCoords();

        DOFVector<int>::Iterator intPointIt(interpPointInd, USED_DOFS);
        typename DOFVector<S>::Iterator valueIt(values, USED_DOFS);
        DOFVector<std::list<WorldVector<double>>>::Iterator coordIt(dofCoords, USED_DOFS);

        file << std::scientific;
        file.precision(15);

        // Write the values for all vertex DOFs.
        for (intPointIt.reset(), valueIt.reset(), coordIt.reset();
             !intPointIt.end();
             ++intPointIt, ++valueIt, ++coordIt)
        {

          if (*intPointIt == -2)
            for (int i = 0; i < static_cast<int>(coordIt->size()); i++)
            {
              file << " " ;
              print(*valueIt,file);
              file << "\n";
            }
        }

        // For the second dim case, write also the values of the interpolation points.
        if ((dim == 2) && (degree > 1))
        {
          DOFVector<std::list<WorldVector<double>>>::Iterator
          interpCoordIt((*dataCollector)[componentNo]->getInterpPointCoords(), USED_DOFS);

          for (intPointIt.reset(), valueIt.reset(), interpCoordIt.reset();
               !intPointIt.end();
               ++intPointIt, ++valueIt, ++interpCoordIt)
          {

            if (*intPointIt >= 0)
            {
              for (unsigned int i = 0; i < interpCoordIt->size(); i++)
              {
                file << " ";
                print(*valueIt,file);
                file << "\n";
              }
            }
          }
        }
      }

      template<typename S> template<typename T>
      void Aux<S>::writeConnectivity(T& file)
      {
        using ::AMDiS::io::VtkWriter::detail::writeConnectivity_dim2_degree2;
        using ::AMDiS::io::VtkWriter::detail::writeConnectivity_dim2_degree3;
        using ::AMDiS::io::VtkWriter::detail::writeConnectivity_dim2_degree4;

        // For the second dim case, and if higher order Lagrange elements are used,
        // write the connectivity by extra functions.
        if ((dim == 2) && (degree == 2))
        {
          writeConnectivity_dim2_degree2((*dataCollector)[0], file);
        }
        else if ((dim == 2) && (degree == 3))
        {
          writeConnectivity_dim2_degree3((*dataCollector)[0], file);
        }
        else if ((dim == 2) && (degree == 4))
        {
          writeConnectivity_dim2_degree4((*dataCollector)[0], file);
        }
        else
        {
          std::list<ElementInfo>* elements = (*dataCollector)[0]->getElementInfos();
          std::list<ElementInfo>::iterator elementIt;
          int vertices = (*dataCollector)[0]->getMesh()->getGeo(VERTEX);

          for (elementIt = elements->begin(); elementIt != elements->end(); ++elementIt)
          {
            // for all vertices
            for (int i = 0; i < vertices; i++)
            {
              file << " " << elementIt->vertexInfo[i]->outputIndex;
            }
            file << "\n";
          }
        }
      }

    } // end namespace VtkVectorWriter
  } // end namespace io
} // end namespace AMDiS

#endif // AMDIS_VTKVECTORWRITER_HH
