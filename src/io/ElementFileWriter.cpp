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


#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>

#include "io/VtkWriter.hpp"
#include "ElementFileWriter.h"
#include "BasisFunction.h"
#include "Initfile.h"
#include "Traverse.h"
#include "AdaptInfo.h"

namespace AMDiS
{
  namespace io
  {

    using namespace std;

    ElementFileWriter::ElementFileWriter(string name_,
                                         Mesh* mesh_,
                                         map<int, double>& mapvec)
      : name(name_),
        amdisMeshDatExt(".elem.mesh"),
        vtkExt(".vtu"),
        pvdExt(".pvd"),
        writeAMDiSFormat(0),
        writeVtkFormat(0),
        writeVtkVectorFormat(0),
        writeAs3dVector(true),
        writeParaViewAnimation(0),
        appendIndex(0),
        indexLength(5),
        indexDecimals(3),
        tsModulo(1),
        timestepNumber(-1),
        mesh(mesh_),
        vec(&mapvec),
        vecs(NULL)
    {
      if (name != "")
      {
        Parameters::get(name + "->output->filename", filename);
        Parameters::get(name + "->output->AMDiS format", writeAMDiSFormat);
        Parameters::get(name + "->output->AMDiS mesh-dat ext", amdisMeshDatExt);
        Parameters::get(name + "->output->ParaView format", writeVtkFormat);
        Parameters::get(name + "->output->ParaView animation", writeParaViewAnimation);
        Parameters::get(name + "->output->append index", appendIndex);
        Parameters::get(name + "->output->index length", indexLength);
        Parameters::get(name + "->output->index decimals", indexDecimals);
        Parameters::get(name + "->output->write every i-th timestep", tsModulo);
      }
    }


    ElementFileWriter::ElementFileWriter(string name_,
                                         Mesh* mesh_,
                                         map<int, vector<double>>& mapvec)
      : name(name_),
        amdisMeshDatExt(".elem.mesh"),
        vtkExt(".vtu"),
        pvdExt(".pvd"),
        writeAMDiSFormat(0),
        writeVtkFormat(0),
        writeVtkVectorFormat(0),
        writeAs3dVector(true),
        writeParaViewAnimation(0),
        appendIndex(0),
        indexLength(5),
        indexDecimals(3),
        tsModulo(1),
        timestepNumber(-1),
        mesh(mesh_),
        vec(NULL),
        vecs(&mapvec)
    {
      if (name != "")
      {
        Parameters::get(name + "->output->filename", filename);
        Parameters::get(name + "->output->AMDiS format", writeAMDiSFormat);
        Parameters::get(name + "->output->AMDiS mesh-dat ext", amdisMeshDatExt);
        Parameters::get(name + "->output->ParaView format", writeVtkFormat);
        Parameters::get(name + "->output->ParaView vector format", writeVtkVectorFormat);
        Parameters::get(name + "->output->write vector as 3d vector", writeAs3dVector);
        Parameters::get(name + "->output->ParaView animation", writeParaViewAnimation);
        Parameters::get(name + "->output->append index", appendIndex);
        Parameters::get(name + "->output->index length", indexLength);
        Parameters::get(name + "->output->index decimals", indexDecimals);
        Parameters::get(name + "->output->write every i-th timestep", tsModulo);
      }
    }


    void ElementFileWriter::writeFiles(AdaptInfo& adaptInfo, bool force,
                                       int level, Flag traverseFlag,
                                       bool (*writeElem)(ElInfo*))
    {
      FUNCNAME("ElementFileWriter::writeFiles()");

      timestepNumber++;
      timestepNumber %= tsModulo;

      if (timestepNumber != 0 && !force)
        return;

      string fn = filename;

      if (appendIndex)
      {
        TEST_EXIT(indexLength <= 99)("index lenght > 99\n");
        TEST_EXIT(indexDecimals <= 97)("index decimals > 97\n");
        TEST_EXIT(indexDecimals < indexLength)("index length <= index decimals\n");

        char formatStr[9];
        sprintf(formatStr, "%%0%d.%df", indexLength, indexDecimals);

        char timeStr[20];
        sprintf(timeStr, formatStr, adaptInfo.getTime());

        fn += timeStr;
      }

      if (writeAMDiSFormat)
      {
        TEST_EXIT(mesh)("no mesh\n");

        writeMeshDatValues(const_cast<char*>((fn + amdisMeshDatExt).c_str()),
                           adaptInfo.getTime());
        MSG("MeshDat file written to %s\n", (fn + amdisMeshDatExt).c_str());
      }

      if (writeVtkFormat || writeVtkVectorFormat)
      {
        TEST_EXIT(mesh)("no mesh\n");

        writeVtkValues(fn, vtkExt, -1, writeVtkVectorFormat == 1);
        MSG("VTK file written to %s\n", (fn + vtkExt).c_str());
      }

      if (writeParaViewAnimation)
        VtkWriter::detail::updateAnimationFile(adaptInfo,
                                               (fn + vtkExt),
                                               &paraViewAnimationFrames,
                                               (filename + pvdExt));
    }


    void ElementFileWriter::writeFile(map<int, double>& vec,
                                      Mesh* mesh,
                                      string filename,
                                      string postfix,
                                      int level,
                                      bool writeAsVector)
    {
      ElementFileWriter efw("", mesh, vec);
      efw.writeVtkValues(filename, postfix, level, writeAsVector);
    }


    void ElementFileWriter::writeFile(map<int, vector<double>>& vecs,
                                      Mesh* mesh,
                                      string filename,
                                      string postfix,
                                      int level,
                                      bool writeAsVector)
    {
      ElementFileWriter efw("", mesh, vecs);
      efw.writeVtkValues(filename, postfix, level, writeAsVector);
    }


    void ElementFileWriter::writeMeshDatValues(string filename, double time)
    {
      FUNCNAME("ElementFileWriter::writeMeshDatValues()");

      boost::iostreams::filtering_ostream file;
      {
        //boost::iostreams seems not to truncate the file
        ofstream swapfile(filename.c_str(), ios::out | ios::trunc);
        TEST_EXIT(swapfile.is_open())
        ("Cannot open file %s for writing!\n", name.c_str());
        swapfile.close();
      }
      file.push(boost::iostreams::file_descriptor_sink(filename, ios::trunc));

      int dim = mesh->getDim();
      double val;

      // === Write header. ===
      file << "mesh name: " << mesh->getName().c_str() << "\n\n";
      file << "time: " << time << "\n\n";
      file << "DIM: " << dim << "\n";
      file << "DIM_OF_WORLD: " << Global::getGeo(WORLD) << "\n\n";
      file << "number of vertices: " << (dim+1)*mesh->getNumberOfLeaves() << "\n";
      file << "number of elements: " << mesh->getNumberOfLeaves() << "\n\n";


      // === Write vertex coordinates (every vertex for every element). ===
      file << "vertex coordinates:\n";
      TraverseStack stack;

      ElInfo* elInfo = stack.traverseFirst(mesh,
                                           -1,
                                           Mesh::CALL_LEAF_EL |
                                           Mesh::FILL_COORDS);

      while (elInfo)
      {

        // Write coordinates of all element vertices.
        for (int i = 0; i <= dim; ++i)
        {
          for (int j = 0; j < dim; ++j)
          {
            file << elInfo->getCoord(i)[j] << " ";
          }
          file << "\n";
        }

        elInfo = stack.traverseNext(elInfo);
      }


      // === Write elements. ===
      int numLeaves = mesh->getNumberOfLeaves();
      int vertCntr = 0;
      file << "\n";
      file << "element vertices:\n";
      for (int i = 0; i < numLeaves; ++i)
      {
        for (int j = 0; j <= dim; ++j)
        {
          file << vertCntr << " ";
          ++vertCntr;
        }
        file << "\n";
      }


      // === Write values. ===

      // Write values header.
      file << "\n";
      file << "number of values: 1\n\n";
      file << "value description: " << name.c_str() << "\n";
      file << "number of interpolation points: 0" << "\n";
      file << "type: scalar" << "\n";
      file << "interpolation type: lagrange" << "\n";
      file << "interpolation degree: 1" << "\n";
      file << "end of description: " << name.c_str() << "\n\n";

      // Write values.
      file << "vertex values: " << name.c_str() << "\n";

      file.setf(ios::scientific,ios::floatfield);

      elInfo = stack.traverseFirst(mesh,
                                   -1,
                                   Mesh::CALL_LEAF_EL);

      while (elInfo)
      {
        // Get element value.
        val = (*vec)[elInfo->getElement()->getIndex()];

        // Write value for each vertex of each element.
        for (int i = 0; i <= dim; ++i)
        {
          file << val << "\n";
        }

        elInfo = stack.traverseNext(elInfo);
      }  // end of: mesh traverse


      // Write values trailor.
      file << "\n";
      file << "interpolation values: " << name.c_str() << "\n\n\n";
      file << "element interpolation points: " << name.c_str() << "\n";
    }


    void ElementFileWriter::writeVtkValues(string fname,
                                           string postfix,
                                           int level,
                                           bool writeAsVector)
    {
      FUNCNAME("ElementFileWriter::writeVtkValues()");

      TEST_EXIT((vec!=NULL || vecs!=NULL) && (vec==NULL || vecs==NULL))
      ("Ether vec or vecs must be given, not both and not nothing!");

#if HAVE_PARALLEL_DOMAIN_AMDIS
      string filename =
        fname  +
        "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-" +
        postfix;
#else
      string filename = fname + postfix;
#endif

      boost::iostreams::filtering_ostream file;
      {
        //boost::iostreams seems not to truncate the file
        ofstream swapfile(filename.c_str(), ios::out | ios::trunc);
        TEST_EXIT(swapfile.is_open())
        ("Cannot open file %s for writing!\n", name.c_str());
        swapfile.close();
      }
      file.push(boost::iostreams::file_descriptor_sink(filename, ios::trunc));

      int dim = mesh->getDim();
      int vertices = mesh->getGeo(VERTEX);
      int nVertices = 0;
      int nElements = 0;
      Flag traverseFlag = level == -1 ? Mesh::CALL_LEAF_EL : Mesh::CALL_EL_LEVEL;

      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(mesh, level, traverseFlag);
      while (elInfo)
      {
        nElements++;
        elInfo = stack.traverseNext(elInfo);
      }

      // === Write vertex coordinates (every vertex for every element). ===
      elInfo = stack.traverseFirst(mesh, level, traverseFlag | Mesh::FILL_COORDS);

      std::map<DegreeOfFreedom, std::pair<int, WorldVector<double>>> dof2Idx;
      std::vector<DegreeOfFreedom> idx2Dof;
      while (elInfo)
      {
        // Write coordinates of all element vertices.
        for (int i = 0; i <= dim; i++)
        {
          DegreeOfFreedom idx = elInfo->getElement()->getDof(i, 0);
          if (dof2Idx.find(idx) == dof2Idx.end())
          {
            dof2Idx[idx] = std::make_pair(nVertices, elInfo->getCoord(i));
            idx2Dof.push_back(idx);
            nVertices++;
          }
        }
        elInfo = stack.traverseNext(elInfo);
      }

      file << "<?xml version=\"1.0\"?>\n";
      file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
      file << "  <UnstructuredGrid>\n";
      file << "    <Piece NumberOfPoints=\"" << nVertices
           << "\" NumberOfCells=\"" <<  nElements << "\">\n";
      file << "      <Points>\n";
      file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

      for (int i = 0; i < nVertices; i++)
      {
        DegreeOfFreedom dof = idx2Dof[i];
        writeCoord(file, dof2Idx[dof].second);
      }

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
      elInfo = stack.traverseFirst(mesh, level, traverseFlag);
      while (elInfo)
      {
        // Write coordinates of all element vertices.
        for (int i = 0; i <= dim; i++)
        {
          DegreeOfFreedom dof = elInfo->getElement()->getDof(i, 0);
          file << (dof2Idx[dof].first) << " ";
        }
        file << "\n";
        elInfo = stack.traverseNext(elInfo);
      }
      file << "        </DataArray>\n";
      file << "      </Cells>\n";

      int dataLength = (vecs != NULL ? (*(vecs->begin())).second.size() : 1);
      int nComponents = (!writeAsVector || (vecs == NULL && vec != NULL) ? 1 : dataLength);
      int nDataArrays = (!writeAsVector && (vec == NULL && vecs != NULL) ? dataLength : 1);
      file << "      <CellData>\n";
      for (int i = 0; i < nDataArrays; i++)
      {
        file << "        <DataArray type=\"Float32\" Name=\"value"<<i<<"\" format=\"ascii\" NumberOfComponents=\""<<(writeAsVector ? std::max(3,nComponents) : nComponents)<<"\">\n";

        file.setf(ios::scientific,ios::floatfield);

        elInfo = stack.traverseFirst(mesh, level, traverseFlag);
        int vc = 0;
        while (elInfo)
        {
          // Get element value.
          int idx = elInfo->getElement()->getIndex();

          for (int j = 0; j < nComponents; j++)
          {
            double val = (vec != NULL ? (*vec)[idx] : (static_cast<int>((*vecs)[idx].size())==dataLength ? (*vecs)[idx][i*nComponents+j] : 0.0));

            // Write value for each vertex of each element.
            if (fabs(val) < 1.e-40)
              file << "0.0";
            else
              file << val;

            if (j < nComponents-1)
              file << " ";
          }
          if (writeAs3dVector && writeAsVector && vecs != NULL)
          {
            for (int j = nComponents; j < 3; j++)
              file << " 0.0";
          }
          file << "\n";

          vc++;
          elInfo = stack.traverseNext(elInfo);
        }

        file << "        </DataArray>\n";
      }

      file << "      </CellData>\n";
      file << "    </Piece>\n";
      file << "  </UnstructuredGrid>\n";
      file << "</VTKFile>\n";

#if HAVE_PARALLEL_DOMAIN_AMDIS
      if (MPI::COMM_WORLD.Get_rank() == 0)
      {
        vector<string> componentNames;
        componentNames.push_back("elvalue");
        VtkWriter::detail::writeParallelFile(fname + ".pvtu", MPI::COMM_WORLD.Get_size(),
                                             fname, ".vtu", componentNames);
      }
#endif
    }

  }
} // end namespace io, AMDiS


