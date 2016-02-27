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



/** \file VtkVectorWriter.h */

#ifndef AMDIS_VTKVECTORWRITER_DETAIL_H
#define AMDIS_VTKVECTORWRITER_DETAIL_H

#include "io/DataCollector.hpp"
#include "SystemVector.hpp"

#ifdef HAVE_COMPRESSION
#include "io/FileCompression.hpp"
#endif

namespace AMDiS
{
  namespace io
  {
    namespace VtkVectorWriter
    {

      template<typename S>
      struct Aux
      {
        Aux(std::vector<DataCollector<S>*>* dc, bool writeAs3dVector_=false)
          : dataCollector(dc)
#ifdef HAVE_COMPRESSION
          , compress(NONE)
#endif
          , writeAs3dVector(writeAs3dVector_)
          , three(3)
        {
          degree = (*dataCollector)[0]->getFeSpace()->getBasisFcts()->getDegree();
          dim = (*dataCollector)[0]->getMesh()->getDim();
        }

        /// Writes a ParaView-VTK file.
        int writeFile(std::string name);

        /// Writes a pvtu file, which contains the links to all the rank files.
        void writeParallelFile(std::string name, int nRanks,
                               std::string fnPrefix, std::string fnPostfix);


#ifdef HAVE_COMPRESSION
        /// Set a compressing method for file output.
        void setCompression(FileCompression c)
        {
          compress = c;
        }
#endif

        void setWriteAs3dVector(bool w)
        {
          writeAs3dVector = w;
        }

      protected:
        /// Writes the VTK file to an arbitrary stream.
        template<typename T>
        void writeFileToStream(T& file);

        /// Writes all coordinates of vertices and interpolation points to an output file.
        template<typename T>
        void writeVertexCoords(T& file);

        /// Writes all values of vertices and interpolation point to an output file.
        template<typename T>
        void writeVertexValues(T& file, int componentNo);

        /// Writes the connectivity of all simplices to an output file.
        template<typename OutputStream>
        void writeConnectivity(OutputStream& file);


        std::ostream& print(const WorldVector<double>& l, std::ostream& o)
        {
          using size_type = Size_t< WorldVector<double> >;
          for (size_type i = 0; i < l.getSize(); i++)
          {
            o << (fabs(l[i]) < 1e-40 ? 0.0 : l[i]) << " ";
          }
          for (size_type i = l.getSize(); i < 3; i++)
            o << 0.0 << " ";
          return o;
        }

        std::ostream& print(const std::vector<double>& l, std::ostream& o)
        {
          using size_type = Size_t< WorldVector<double> >;
          for (size_type i = 0; i < l.size(); i++)
          {
            o << (fabs(l[i]) < 1e-40 ? 0.0 : l[i]) << " ";
          }
          if (writeAs3dVector)
          {
            for (size_type i = l.size(); i < 3; i++)
              o << 0.0 << " ";
          }
          return o;
        }

        template<typename T>
        std::ostream& print(const T& value, std::ostream& o)
        {
          o << (fabs(value) < 1e-40 ? 0.0 : value);
          return o;
        }

        // num_rows for scalar types
        template<typename T>
        size_t num_rows(T& /*v*/)
        {
          return 1;
        }
        // for WorldVectors print 3 components
        template<typename T>
        size_t num_rows(WorldVector<T>& /*v*/)
        {
          return three; //(writeAs3dVector ? three : static_cast<size_t>(v.getSize()));
        }
        // num_rows for vector types
        template<typename T>
        size_t num_rows(std::vector<T>& v)
        {
          size_t v_size = v.size();
          return (writeAs3dVector ? std::max(three,v_size) : v_size);
        }
        template<typename T>
        size_t num_rows(const std::vector<T>& v)
        {
          size_t v_size = v.size();
          return (writeAs3dVector ? std::max(three,v_size) : v_size);
        }
        template<typename T>
        size_t num_rows(mtl::dense_vector<T>& v)
        {
          size_t v_size = mtl::num_rows(v);
          return (writeAs3dVector ? std::max(three,v_size) : v_size);
        }
        template<typename T>
        size_t num_rows(const mtl::dense_vector<T>& v)
        {
          size_t v_size = mtl::num_rows(v);
          return (writeAs3dVector ? std::max(three,v_size) : v_size);
        }

      private:
        /// List of DataCollectors, for each component of the problem one.
        std::vector<DataCollector<S>*>* dataCollector;

#ifdef HAVE_COMPRESSION
        /** \brief
        * Defines if the file has to be compressed for ouput, and with which
        * kind of compress method.
        */
        FileCompression compress;
#endif

        /// Degree of the basis function of the problem.
        int degree;

        /// Dimension of the geometry.
        int dim;

        bool writeAs3dVector;
        size_t three;
      };

    } // end namespace VtkVectorWriter
  } // end namespace io
} // end namespace AMDiS

#include "VtkVectorWriter.hh"

#endif
