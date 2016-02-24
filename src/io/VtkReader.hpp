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


/** \file VtkReader.h */

#ifndef AMDIS_VTKREADER_H
#define AMDIS_VTKREADER_H

// need some extension-methods/libs (pugixml, nanoflann)
#ifdef HAVE_EXTENSIONS

#include <cstring>
#include "DOFVector.h"
#include "SystemVector.h"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     * \brief Reader for the ParaView VTU-format
     *
     * A collection of methods to read to various container types from
     * VTU-files.
     **/
    namespace VtkReader
    {
      /// Read  various components of file with given \ref componentNames to a
      /// vector of \ref DOFVector.
      template<typename T>
      void readByName(std::string filename,
                      std::vector<DOFVector<T>*> dofVectors,
                      std::vector<std::string> componentNames);


      /// Read a component of file with given \ref componentName to a \ref DOFVector.
      template<typename T>
      void readByName(std::string filename,
                      DOFVector<T>* dofVector,
                      std::string componentName)
      {
        std::vector<DOFVector<T>*> dofVectors;
        dofVectors.push_back(dofVector);
        std::vector<std::string> componentNames;
        componentNames.push_back(componentName);

        readByName(filename, dofVectors, componentNames);
      }


      /// Interface for general containers not implemented. Specializations below.
      template<typename Container>
      void readFile(std::string filename, Container& vec)
      {
        ERROR_EXIT("VtkReader not implemented for this container type!\n");
      }


      /// Read a file to a \ref DOFVector. Using container pointer.
      template<typename T>
      void readFile(std::string filename,
                    DOFVector<T>* dofVector)
      {
        readByName(filename, dofVector, dofVector->getName());
      }


      /// Read a file to a \ref DOFVector. Using container reference.
      template<typename T>
      void readFile(std::string filename,
                    DOFVector<T>& dofVector)
      {
        readFile(filename, &dofVector);
      }


      /// Read a file to a vector of \ref DOFVector.
      template<typename T>
      void readFile(std::string filename,
                    std::vector<DOFVector<T>*> dofVectors)
      {
        std::vector<std::string> componentNames(dofVectors.size());
        for (size_t i = 0; i < dofVectors.size(); i++)
          componentNames[i] = dofVectors[i]->getName();
        readByName(filename, dofVectors, componentNames);
      }


      /// Read a file to a \ref SystemVector. Using container pointer.
      inline
      void readFile(std::string filename,
                    SystemVector* vecs)
      {
        std::vector<DOFVector<double>*> dofVectors(vecs->getSize());
        std::vector<std::string> componentNames(vecs->getSize());
        for (int i = 0; i < vecs->getSize(); i++)
        {
          dofVectors[i] = vecs->getDOFVector(i);
          componentNames[i] = vecs->getDOFVector(i)->getName();
        }
        readByName(filename, dofVectors, componentNames);
      }

      /// Read a file to a \ref SystemVector. Using container reference.
      inline
      void readFile(std::string filename,
                    SystemVector& vecs)
      {
        readFile(filename, &vecs);
      }

    } // end namespace VtkReader
  } // end namespace io
} // end namespace AMDiS

#include "VtkReader.hh"

#endif // HAVE_EXTENSIONS

#endif // AMDIS_VTKREADER_H
