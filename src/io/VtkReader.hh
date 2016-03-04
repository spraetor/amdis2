#pragma once

// need some extension-methods/libs (pugixml, nanoflann)
#ifdef HAVE_EXTENSIONS

#include <boost/filesystem.hpp>

#include <pugixml.hpp>

// extension {
#include "kdtree_nanoflann.h"
#include "VectorOperations.h"
// }

#include "Mesh.hpp"
#include "io/detail/VtkReader.hpp"

namespace AMDiS
{
  namespace io
  {
    namespace VtkReader
    {

      template<typename T>
      void readByName(std::string filename,
                      std::vector<DOFVector<T>*> dofVectors,
                      std::vector<std::string> componentNames)
      {
        using namespace std;
        using namespace pugi;
        using namespace AMDiS::io::VtkReader;

        FUNCNAME("VtuReader::readValue()");

        TEST_EXIT(filename != "")("Filename not specified!\n");
        TEST_EXIT(dofVectors.size() > 0)("no DOF vectors specified\n");
        TEST_EXIT(componentNames.size() > 0)("no componentName specified\n");

        TEST_EXIT(boost::filesystem::exists(filename))((filename + " does not exist!\n").c_str());

        xml_document vtu;
        TEST_EXIT(vtu.load_file(filename.c_str()))("Could not load vtu file! Error in xml structure.\n");

        xml_node VTKFile = vtu.child("VTKFile");
        string zlib = VTKFile.attribute("compressor").value();
        TEST_EXIT(zlib == "" || zlib == "vtkZLibDataCompressor")("The only supported compressor is Zlib.\n Should not happen.\n ");

        xml_node UnstructuredGrid = VTKFile.child("UnstructuredGrid");
        xml_node Piece = UnstructuredGrid.child("Piece");
        xml_node Points = Piece.child("Points");
        xml_node PointsDataArray = Points.child("DataArray");
        xml_node AppendedData = VTKFile.child("AppendedData");

        typedef vector<WorldVector<double>> PL;
        PL pointList;
        vector<vector<T>> valueList(dofVectors.size());

        if(!AppendedData)
        {

          xml_node PointData = Piece.child("PointData");
          string points = Points.child_value("DataArray");
          detail::string2pointList(points, pointList);

          T test;

          for (xml_node DataArray = PointData.child("DataArray"); DataArray; DataArray = DataArray.next_sibling("DataArray"))
          {
            string Name = DataArray.attribute("Name").value();
            for (size_t i = 0; i < componentNames.size(); i++)
            {
              if (Name == componentNames[i])
              {
                string values = DataArray.last_child().value();
                int nComponents = -1;
                if (DataArray.attribute("NumberOfComponents"))
                  nComponents = DataArray.attribute("NumberOfComponents").as_int();
                TEST_EXIT(nComponents == -1 || static_cast<int>(size(test)) <= nComponents)
                ("Can not read values in DOFVector with given value type. Too many components!\n");
                detail::string2valueList(values, valueList[i], size(test), nComponents);
                break;
              }
            }
          }
        }
        else
        {
          string encoding = AppendedData.attribute("encoding").value();
          TEST_EXIT(encoding == "base64")
          ("Currently the encoding of AppendedData only supports base64. But it's easy to extend to raw.\n");

          string appendedData = AppendedData.last_child().value();
          int start = appendedData.find("_");
          appendedData = appendedData.substr(start + 1, appendedData.length() - 1);

          vector<int> offsetVec;
          int pointsIndex = 0, index = 0;
          vector<pair<xml_node, int>> pointDataIndex;
          for (xml_node Parent = Piece.first_child(); Parent; Parent = Parent.next_sibling())
          {
            for(xml_node Data = Parent.child("DataArray"); Data; Data = Data.next_sibling("DataArray"), index++)
            {
              string parentName = Data.parent().name();
              if(parentName == "Points")
                pointsIndex = index;
              if(parentName == "PointData")
                pointDataIndex.push_back(make_pair(Data, index));

              offsetVec.push_back(Data.attribute("offset").as_int());
            }
          }

          string format = PointsDataArray.attribute("format").value();
          string type = PointsDataArray.attribute("type").value();
          TEST_EXIT(format == "appended")("The format of DataArray is not appended in appended mode. Should not happen.\n");
          TEST_EXIT(type == "Float32" || type == "Float64")("Currently only supports Points with type Float32 and Float64.\n");

          int len = ((unsigned)(pointsIndex + 1) == offsetVec.size()) ?
                    appendedData.length() - offsetVec[pointsIndex] :
                    offsetVec[pointsIndex + 1] - offsetVec[pointsIndex];

          string points = appendedData.substr(offsetVec[pointsIndex], len);
          detail::binary2pointList(points, type, (zlib != ""), pointList);

          T test;

          for (size_t i = 0 ; i < pointDataIndex.size(); i++)
          {
            string Name = pointDataIndex[i].first.attribute("Name").value();
            for (size_t j = 0; j < componentNames.size(); j++)
            {
              if(Name == componentNames[j])
              {
                format = pointDataIndex[i].first.attribute("format").value();
                type = pointDataIndex[i].first.attribute("type").value();
                TEST_EXIT(format == "appended")("The format of DataArray is not appended in appended mode. Should not happen.\n");
                TEST_EXIT(type == "Float32" || type == "Float64")("Currently only supports PointData with type Float32 and Float64.\n");

                len = ((unsigned)(pointDataIndex[i].second + 1) == offsetVec.size())?
                      appendedData.length() - offsetVec[pointDataIndex[i].second]:
                      offsetVec[pointDataIndex[i].second + 1] - offsetVec[pointDataIndex[i].second];

                string values = appendedData.substr(offsetVec[pointDataIndex[i].second], len);

                int nComponents = -1;
                if (pointDataIndex[i].first.attribute("NumberOfComponents"))
                  nComponents = pointDataIndex[i].first.attribute("NumberOfComponents").as_int();
                TEST_EXIT(nComponents == -1 || static_cast<int>(size(test)) <= nComponents)
                ("Can not read values in DOFVector with given value type. Too many components!\n");
                detail::binary2valueList(values, type, (zlib != ""), valueList[j], size(test), nComponents);
                break;
              }
            }
          }
        }
        for (size_t i = 0; i < dofVectors.size(); i++)
        {
          TEST_EXIT(valueList[i].size() != 0)
          ("No values found for component name '%s'!\n", componentNames[i].c_str());
          TEST_EXIT(pointList.size() == valueList[i].size())
          ("Coordinates and values[%d] do not match!\n", i);
        }

        DOFVector<WorldVector<double>> coords(dofVectors[0]->getFeSpace(), "coords");
        dofVectors[0]->getFeSpace()->getMesh()->getDofIndexCoords(coords);

        extensions::KD_Tree tree(Global::getGeo(WORLD), pointList, 10);
        tree.index->buildIndex();

        DOFIterator<WorldVector<double>> it_coords(&coords, USED_DOFS);
        vector<DOFIterator<T>*> it_results;
        for (size_t i = 0; i < dofVectors.size(); i++)
        {
          it_results.push_back(new DOFIterator<T>(dofVectors[i], USED_DOFS));
          it_results[i]->reset();
        }
        for (it_coords.reset(); !it_coords.end(); ++it_coords)
        {
          size_t idx = detail::getNearestIndex(tree, *it_coords);

          for (size_t i = 0; i < dofVectors.size(); i++)
          {
            if (valueList[i].size() > idx)
              *(*it_results[i]) = valueList[i][idx];
            (*it_results[i])++;
          }
        }


        MSG("VTU file read from: %s\n", filename.c_str());
      }

    } // end namespace VtuReader
  } // end namespace io
} // end namespace

#endif // HAVE_EXTENSIONS
