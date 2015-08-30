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


#include "PovrayWriter.h"
#include "DOFVector.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include "Traverse.h"
#include "DataCollector.h"

using namespace std;

namespace AMDiS
{
  namespace io
  {

    /* destructor */
    PovrayWriter::~PovrayWriter()
    {
      // free the memory for the bounding box
      if (bBox)
      {
        delete bBox;
        bBox = NULL;
      }
    }


    /* some forward declarations */
    double computeWeight(double&, double&, double&);
    string getColorString(double&);


    void PovrayWriter::tryMeshTraversal(ofstream& out)
    {
      /* prepare some data structures */
      DOFVector<double>* values = dataCollector->getValues();
      const FiniteElemSpace* feSpace = dataCollector->getFeSpace();
      const BasisFunction* basFcts = feSpace->getBasisFcts();
      std::vector<DegreeOfFreedom> dofs(basFcts->getNumber());

      /* create an 'iterator' */
      TraverseStack stack;
      ElInfo* elInfo =  stack.traverseFirst(dataCollector->getMesh(), -1,
                                            Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);

      /* find min and max values */
      double min_value = 0;
      double max_value = 0;
      if (values->getSize() > 0)
        min_value = max_value = (*values)[0];

      for (int index = 1; index < values->getSize(); index++)
      {
        min_value = std::min(min_value, (*values)[index]);
        max_value = std::max(max_value, (*values)[index]);
      }

      /* map DOFs to values */
      std::map <DegreeOfFreedom, double> value_map;
      while (elInfo)
      {
        basFcts->getLocalIndices(elInfo->getElement(), feSpace->getAdmin(), dofs);

        int max = Global::getGeo(WORLD)==2 ? 3 : 4;

        for (int i = 0; i < max; i++)
        {
          double value = (*values)[dofs[i]];
          value_map[dofs[i]] = value;
        }

        elInfo = stack.traverseNext(elInfo);
      }

      out << "\ttexture_list {" << endl;
      out << "\t" << value_map.size() << "," << endl;
      std::map<DegreeOfFreedom, double>::iterator map_iter;
      for (map_iter = value_map.begin(); map_iter != value_map.end(); ++map_iter)
      {
        double value = map_iter->second;
        double color_weight =  computeWeight(min_value,max_value,value);
        out << "\t\ttexture{ pigment{ " << getColorString(color_weight) << "}}";
        std::map<DegreeOfFreedom, double>::iterator final_iter = value_map.end();
        --final_iter;
        if (map_iter != final_iter)
          out << ",";
        out << endl;
      }
      out << "\t}" << endl; // end of texture_list
    }

    // provides the bounding box of the mesh
    //
    // note: computation is done during the first call other, this might cause
    //       problems in simulations with changing geometries
    BoundingBox* PovrayWriter::getBoundingBox(ofstream& out)
    {
      if (bBox != NULL)
        return bBox;

      DOFVector<std::list<VertexInfo>>* vertexInfos = dataCollector->getVertexInfos();
      DOFVector<std::list<VertexInfo>>::Iterator it(vertexInfos, USED_DOFS);
      int dow = Global::getGeo(WORLD);
      bBox = new BoundingBox();

      it.reset();
      std::list<VertexInfo>::iterator it2 = it->begin();
      bBox->minx = it2->coords[0];
      bBox->maxx = it2->coords[0];
      bBox->miny = it2->coords[1];
      bBox->maxy = it2->coords[1];
      if (dow == 3)
      {
        bBox->minz = it2->coords[2];
        bBox->maxz = it2->coords[2];
      }
      else
      {
        bBox->minz = bBox->maxz = 0;
      }

      for (it.reset(); !it.end(); ++it)
      {
        // initialize mit and max values with coordinates of first vertex
        std::list<VertexInfo>::iterator it2 = it->begin();
        bBox->minx = std::min(bBox->minx, it2->coords[0]);
        bBox->maxx = std::max(bBox->maxx, it2->coords[0]);
        bBox->miny = std::min(bBox->miny, it2->coords[1]);
        bBox->maxy = std::max(bBox->maxy, it2->coords[1]);
        if (dow == 3)
        {
          bBox->minz = std::min(bBox->minz, it2->coords[2]);
          bBox->maxz = std::max(bBox->maxz, it2->coords[2]);
        }
      }

      return bBox;
    }

    // computes the weight of a value for color computation
    double computeWeight(double& min, double& max, double& value)
    {
      double interval_length = (max - min);
      return (value - min) / interval_length;
    }

    // determines the color for a given value/weight
    string getColorString(double& value)
    {
      value = std::max(0.0, value);
      value = std::min(value, 1.0);

      // rot: 1,0,0
      // blau: 0,0,1

      // currently blue->red gradient
      stringstream res;
      res << "rgb<" << value << ",0," << (1 - value) << ">";
      return res.str();
    }

    // writes a povray script for the current time step to the specified file
    void PovrayWriter::writeFile(std::string filename)
    {
      ofstream out;
      out.open(filename.c_str(), ios::out);

      // delegate all the work to other subroutines
      writeHeader(out);
      writeIncludes(out);
      writeCamera(out);
      writeLight(out);
      writeTestStuff(out, *dataCollector);
      writeMesh2(out, *dataCollector);

      out.close();
    }


    // writes a simple info header (comments only)
    void PovrayWriter::writeHeader(ofstream& out)
    {
      out << "// povray file created by AMDiS" << endl;
      out << endl;
    }


    // writes all neccessary include file declarations
    void PovrayWriter::writeIncludes(ofstream& out)
    {
      out << "#include \"colors.inc\"" << endl;
      out << "#include \"stones.inc\"" << endl;
      out << endl;
    }

    // computes an appropriate camera position and writes it to the povray file
    void PovrayWriter::writeCamera(ofstream& out)
    {
      /* TODO: compute adequate camera position (computation needs to be performed
         only once, same camera position will be used for all iterations) */

      /* TODO: check whether camera position and location have been specified in the
         init file, use them if present*/

      BoundingBox* box = getBoundingBox(out);
      int dow = Global::getGeo(WORLD);
      double centerx = (box->minx+box->maxx)/2.;
      double centery = (box->miny+box->maxy)/2.;
      double centerz = (box->minz+box->maxz)/2.;

      if (dow == 2)
      {
        out << "camera {" << endl;
        out << "\tlocation <" << centerx << ", " << centery << ", -3>" << endl;
        out << "\tlook_at  <" << centerx << ", " << centery << ",  0>" << endl;
        out << "}" << endl << endl;
      }
      else
      {
        out << "camera {" << endl;
        out << "\tlocation <"<< centerx <<", "<<centery<<", -3>" << endl;
        out << "\tlook_at  <"<< centerx <<", "<<centery<<",  "<<centerz<<">" << endl;
        out << "}" << endl << endl;
      }
    }


    // computes an appropriate light source and writes it to the povray file
    void PovrayWriter::writeLight(ofstream& out)
    {
      // TODO: compute adequate light position
      out << "light_source { <0, 0, -15> color White}" << endl;
    }


    // writes the data in the 'mesh2' format
    void PovrayWriter::writeMesh2(ofstream& out, DataCollector<>& dataCollector)
    {
      // initialization and preparations
      //Mesh *mesh = dataCollector.getMesh();

      // begin of mesh2 block
      out << "mesh2 {" << endl;

      // delegate work to other methods
      writeVertexVectors(out, dataCollector);
      //writeTextureList(out, dataCollector);
      tryMeshTraversal(out);
      writeFaceIndices(out, dataCollector);
      out << "\t" << "pigment {rgb 1}" << endl;
      // end of mesh2 block
      out << "}" << endl;
    }


    // writes all vertices of the model
    void PovrayWriter::writeVertexVectors(ofstream& out, DataCollector<>& dataCollector)
    {
      // initialization and preparations
      Mesh* mesh = dataCollector.getMesh();
      DOFVector<std::list<VertexInfo>>* vertexInfos = dataCollector.getVertexInfos();
      DOFVector<std::list<VertexInfo>>::Iterator it(vertexInfos, USED_DOFS);

      // begin of vertex_vectors block
      out << "\tvertex_vectors {" << endl;
      int nVertices = mesh->getNumberOfVertices();
      out << "\t\t" << nVertices << "," << endl;

      // For all DOFs of vertices, write the coordinates.
      for (it.reset(); !it.end(); ++it)
      {
        // for all vertex infos of this DOF
        for (std::list<VertexInfo>::iterator it2 = it->begin(); it2 != it->end(); ++it2)
        {
          out << endl;
          // TODO: right implementation using STL vector/iterator
          out  << "\t\t<";
          WorldVector<double> coords;
          coords = it2->coords;
          for (int i = 0; i < Global::getGeo(WORLD); i++)
          {
            if (i > 0)
              out << ", ";

            out << std::scientific << coords[i];
          }
          for (int i = Global::getGeo(WORLD); i < 3; i++)
            out << ", 0.0";
          out << ">,";
        }
      }
      // undo last comma
      long pos=out.tellp(); //tells pos the actual stream position
      out.seekp(pos - 1);     //sets stream position -1
      out <<  endl;         //overwrites the las comma with space
      // end of vertex_vectors block
      out << "\t}" << endl;
    }


    // writes the list of textures (basically one color for each vertex, used to draw a color gradient on each triangle)
    void PovrayWriter::writeTextureList(ofstream& out, DataCollector<>& dataCollector)
    {
      // prepare some data structures
      Mesh* mesh = dataCollector.getMesh();
      DOFVector<std::list<VertexInfo>>* vertexInfos = dataCollector.getVertexInfos();
      DOFVector<std::list<VertexInfo>>::Iterator it(vertexInfos, USED_DOFS);


      // begin of texture_list block
      out << "\ttexture_list {" << endl;

      // compute and write the number of textures/colors (one for each vertex)
      int nTextures = mesh->getNumberOfVertices();
      out << "\t\t" << nTextures << "," << endl;

      // For all DOFs of vertices, write the coordinates.
      for (it.reset(); !it.end(); ++it)
      {
        // for all vertex infos of this DOF
        for (std::list<VertexInfo>::iterator it2 = it->begin(); it2 != it->end(); ++it2)
        {
          // test: use y-coordinate to compute color
          double redValue = it2->coords[1];
          redValue = std::max(0., redValue);
          redValue = std::min(redValue, 1.);

          out << "\n\t\ttexture{ pigment{ rgb"<< getColorString(redValue) <<" } }" << ",";
        }
      }
      long pos = out.tellp();
      out.seekp(pos - 1);
      out << endl;
      // end of texture_list block
      out << "\t}" << endl;
    }


    // writes the connectivity
    void PovrayWriter::writeFaceIndices(ofstream& out, DataCollector<>& dataCollector)
    {
      // begin of face_indices block
      out << "\tface_indices {" << endl;

      std::list<ElementInfo>* elementInfos = dataCollector.getElementInfos();
      int vertices = dataCollector.getMesh()->getGeo(VERTEX);
      int meshDim = dataCollector.getMesh()->getDim();
      int counter = 0;

      if (meshDim == 2)
      {
        int nElements = dataCollector.getElementInfos()->size();
        out << "\t\t" << nElements << "," << endl;
        for (std::list<ElementInfo>::iterator iter = elementInfos->begin();
             iter != elementInfos->end(); iter++)
        {
          out << endl;
          out << "\t\t<";
          for (int i = 0; i < vertices; i++)
          {
            if (i > 0)
              out << ", ";

            out << iter->vertexInfo[i]->outputIndex;
          }
          out << ">";

          // write indices again for texture/color mapping
          for (int i = 0; i < vertices; i++)
          {
            out << ", ";
            out << iter->vertexInfo[i]->outputIndex;
          }

          out << ",";
        }
        //Delete last comma
        long pos = out.tellp();
        out.seekp(pos - 1);
        out << endl;
      }
      else if (meshDim == 3)
      {
        int nElements = dataCollector.getElementInfos()->size() * 4;
        out << "\t\t" << nElements << "," << endl;
        for (std::list<ElementInfo>::iterator iter = elementInfos->begin();
             iter != elementInfos->end(); iter++)
        {
          for (int j = 0; j < vertices; j++)
          {
            out << "/* " << counter++ << " */";
            out << "\t\t<";
            for (int i = 0; i < vertices; i++)
            {
              if (i==j)
                continue;

              if (!((j > 0 && i == 0) || (j == 0 && i == 1)))
                out << ", ";

              out << iter->vertexInfo[i]->outputIndex;
            }
            out << ">";
            // write indices again for texture/color mapping
            for (int i = 0; i < vertices; i++)
            {
              if (i == j)
                continue;
              out << ", ";
              out << iter->vertexInfo[i]->outputIndex;
            }
            out << "," << endl;
          }
        }
      }

      // end of face_indices block
      out << "\n\t}" << endl; //
    }



    /// writes some additional information for debugging purposes
    void PovrayWriter::writeTestStuff(ofstream& out, DataCollector<>& dataCollector)
    {
      Mesh* mesh = dataCollector.getMesh();

      out << "/*" << endl;
      out << "TestData:" << endl;

      int dow = Global::getGeo(WORLD);
      out << " world dimension:\t" << dow << endl;

      out << " dimension:\t" << mesh->getDim() << endl;
      out << " numberOfNodes:\t" << mesh->getNumberOfNodes() << endl;
      out << " numberOfEdges:\t" << mesh->getNumberOfEdges() << endl;
      out << " numberOfVertices:\t" << mesh->getNumberOfVertices() << endl;
      out << " numberOfFaces:\t" << mesh->getNumberOfFaces() << endl;
      out << " numberOfElements:\t" << mesh->getNumberOfElements() << endl;
      out << "*/" << endl << endl;
    }

  }
} // end of namespace io, AMDiS
