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


/** \file ReaderWriter.h */

#ifndef AMDIS_READER_DETAIL_H
#define AMDIS_READER_DETAIL_H

#include <cstring>
#include "DOFVector.h"
#include "SystemVector.h"

namespace AMDiS
{
  namespace io
  {
    namespace detail
    {

      template<typename T>
      Mesh* getMesh(DOFVector<T>& container)
      {
        return container.getFeSpace()->getMesh();
      }


      template<typename T>
      Mesh* getMesh(std::vector<DOFVector<T>*>& container)
      {
        return container[0]->getFeSpace()->getMesh();
      }


      inline
      Mesh* getMesh(SystemVector& container)
      {
        return container.getDOFVector(0)->getFeSpace()->getMesh();
      }


      inline
      Mesh* getMesh(Mesh& container)
      {
        return &container;
      }



      template<typename T>
      const FiniteElemSpace* getFeSpace(DOFVector<T>& container, int /*i*/ = 0)
      {
        return container.getFeSpace();
      }


      template<typename T>
      const FiniteElemSpace* getFeSpace(std::vector<DOFVector<T>*>& container, int i = 0)
      {
        return container[i]->getFeSpace();
      }


      inline
      const FiniteElemSpace* getFeSpace(SystemVector& container, int i = 0)
      {
        return container.getDOFVector(i)->getFeSpace();
      }


      inline
      const FiniteElemSpace* getFeSpace(Mesh& /*container*/, int /*i*/ = 0)
      {
        ERROR_EXIT("Can not extract feSpace from Mesh!\n");
        return NULL;
      }


      template<typename T>
      DOFVector<T>* getDOFVector(DOFVector<T>& container, int /*i*/ = 0)
      {
        return &container;
      }


      template<typename T>
      DOFVector<T>* getDOFVector(std::vector<DOFVector<T>*>& container, int i = 0)
      {
        return container[i];
      }


      inline
      DOFVector<double>* getDOFVector(SystemVector& container, int i = 0)
      {
        return container.getDOFVector(i);
      }


      inline
      DOFVector<double>* getDOFVector(Mesh& /*container*/, int /*i*/ = 0)
      {
        ERROR_EXIT("Can not extract DOFVector from Mesh!\n");
        return NULL;
      }

    } // end namespace detail
  } // end namespace io
} // end namespace AMDiS

#endif // AMDIS_READER_DETAIL_H
