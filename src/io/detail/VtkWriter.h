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



/** \file VtkVector.h */

#ifndef AMDIS_VTKWRITER_DETAIL_H
#define AMDIS_VTKWRITER_DETAIL_H

#include <string>
#include <vector>
#include <sstream>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>

#include <AdaptInfo.h>
#include <io/DataCollector.h>
#include <utility/string.hpp>

#define AMDIS_ZLIB_BLOCK_SIZE 32768

namespace AMDiS
{
  namespace io
  {

    namespace VtkWriter
    {

      typedef enum
      {
        ASCII = 0,
        APPENDED = 1,
        APPENDED_COMPRESSED = 2
      } Vtuformat;

    }
  }

  // conversion function used in Initfile parser
  namespace detail
  {
    /// convert special enums
    inline void convert(const std::string valStr, ::AMDiS::io::VtkWriter::Vtuformat& value)
    {
      using namespace ::AMDiS::io::VtkWriter;
      std::string swapStr = to_upper(valStr);

      value = static_cast<Vtuformat>(ASCII);
      if (swapStr == "ASCII")
        value = static_cast<Vtuformat>(ASCII);
      else if (swapStr == "APPENDED")
        value = static_cast<Vtuformat>(APPENDED);
      else if (swapStr == "APPENDED_COMPRESSED")
        value = static_cast<Vtuformat>(APPENDED_COMPRESSED);
      else
      {
        throw std::runtime_error("Unknown ParaView mode.");
      }
    }
  } // end namespace detail

  namespace io
  {

    namespace VtkWriter
    {

      using namespace std;

      // A class definition of a VtkWriter, to accept Datacollectors in
      // the constructor.
      class Aux
      {
      public:

        class BinaryStream: public stringstream
        {
        public:
          explicit BinaryStream(bool hp = false, ios_base::openmode mode = ios_base::in | ios_base::out | ios_base::binary)
            : stringstream(mode | ios_base::binary)
          {
            highPrecision = hp;
          }

          bool isHighPrecision()
          {
            return highPrecision;
          }

          int getSize()
          {
            seekp(0, ios::end);
            return boost::lexical_cast<int>(tellp());
          }

        private:
          /// True: write data using FLOAT64; False: write data using FLOAT32.
          bool highPrecision;
        };

        Aux(std::vector<DataCollector<>*>* dc,
            std::vector<std::string>& names_,
            Vtuformat f = ASCII,
            bool hp = false,
            bool writeAsVector_ = false)
          : dataCollector(dc),
            componentNames(names_),
            bstream(hp),
            writeAsVector(writeAsVector_)
        {
#ifndef HAVE_COMPRESSION
          FUNCNAME("VtkWriter::Aux::Aux()");
          if(f == APPENDED_COMPRESSED)
          {
            f = APPENDED;
            WARNING("HAVE_COMPRESSION OFF. So vtuformat is automatically changed from APPENDED_COMPRESSED to APPENDED.\n");
          }
#endif
          degree = (*dataCollector)[0]->getFeSpace()->getBasisFcts()->getDegree();
          dim = (*dataCollector)[0]->getMesh()->getDim();
          format = f;
        }

        /// Writes a ParaView-VTK file.
        int writeFile(std::string name);

        void setFormat(Vtuformat f)
        {
          format = f;
        }

      protected:
        /// Writes the VTK file to an arbitrary stream.
        template<typename T>
        void writeFileToStream(T& file);

        /// Writes tje VTL file to an arbitrary stream using appended mode.
        template<typename T>
        void writeFileToStreamAppended(T& file);

        /// Writes all coordinates of vertices and interpolation points to an
        /// output file.
        template<typename T>
        void writeVertexCoords(T& file);

        /// Writes all values of vertices and interpolation point to an output file.
        template<typename T>
        void writeVertexValues(T& file, int componentNo);

        /// Writes all values of vertices and interpolation point as vector to an output file.
        template<typename T>
        void writeVertexValues(T& file);

        /// Writes the connectivity of all simplices to an output file.
        template<typename OutputStream>
        void writeConnectivity(OutputStream& file);

        template<typename Stream>
        Stream& print(const WorldVector<double>& l, Stream& o)
        {
          for (size_t i = 0; i < static_cast<size_t>(l.getSize()); i++)
          {
            o << (fabs(l[i]) < 1e-40 ? 0.0 : l[i]) << " ";
          }
          for (int i = l.getSize(); i < 3; i++)
            o << 0.0 << " ";
          return o;
        }

        template<typename Stream>
        Stream& print(const std::vector<double>& l, Stream& o)
        {
          for (size_t i = 0; i < l.size(); i++)
          {
            o << (fabs(l[i]) < 1e-40 ? 0.0 : l[i]) << " ";
          }
          for (int i = l.size(); i < 3; i++)
            o << 0.0 << " ";
          return o;
        }

        template<typename T, typename Stream>
        Stream& print(const T& value, Stream& o)
        {
          o << (fabs(value) < 1e-40 ? 0.0 : value);
          return o;
        }

      private:
        std::string getStreamData();

        /// List of DataCollectors, for each component of the problem one.
        std::vector<DataCollector<>*>* dataCollector;

        std::vector<std::string> componentNames;

        /// Degree of the basis function of the problem.
        int degree;

        /// Dimension of the geometry.
        int dim;

        /// 0:ascii mode; 1:append mode; 2:appended_compressed mode.
        Vtuformat format;

        /// The binary stream for filtering data. Only used in append and appended_compressed mode.
        BinaryStream bstream;

        /// write data-vectors as vector-valued component
        bool writeAsVector;
      };

      template <class T>
      Aux::BinaryStream& operator<< (Aux::BinaryStream& stream, const T& /*param*/)
      {
        return stream;
      }

      template <>
      inline Aux::BinaryStream& operator<< (Aux::BinaryStream& stream, const double& param)
      {
        double* param_ = const_cast<double*>(&param);

        if(stream.isHighPrecision())
        {
          stream.write(reinterpret_cast<char*>(param_), 8);
        }
        else
        {
          float param__ = static_cast<float>(*param_);
          stream.write(reinterpret_cast<char*>(&param__), 4);
        }
        return stream;
      }

      template <>
      inline Aux::BinaryStream& operator<< (Aux::BinaryStream& stream, const uint64_t& param)
      {
        uint64_t* param_ = const_cast<uint64_t*>(&param);
        stream.write(reinterpret_cast<char*>(param_), 8);
        return stream;
      }

      template <>
      inline Aux::BinaryStream& operator<< (Aux::BinaryStream& stream, const int& param)
      {
        int* param_ = const_cast<int*>(&param);
        stream.write(reinterpret_cast<char*>(param_), 4);
        return stream;
      }

      template <>
      inline Aux::BinaryStream& operator<< (Aux::BinaryStream& stream, const uint32_t& param)
      {
        uint32_t* param_ = const_cast<uint32_t*>(&param);
        stream.write(reinterpret_cast<char*>(param_), 4);
        return stream;
      }

      template <>
      inline Aux::BinaryStream& operator<< (Aux::BinaryStream& stream, const uint8_t& param)
      {
        uint8_t* param_ = const_cast<uint8_t*>(&param);
        stream.write(reinterpret_cast<char*>(param_), 1);
        return stream;
      }

      namespace detail
      {
        typedef
        boost::archive::iterators::base64_from_binary<// convert binary values to base64 characters
        boost::archive::iterators::transform_width<// retrieve 6 bit integers from a sequence of 8 bit bytes
        std::string::const_iterator, 6, 8
        >
        > binary_base64;

        std::string base64Encode(std::string text);

        std::string extract_relative_path(std::string valueFilename, std::string animationFilename);

        /// Adds a new entry to a ParaView animation file.
        int updateAnimationFile(AdaptInfo& adaptInfo,
                                std::string valueFilename,
                                std::vector<std::pair<double, std::string>>* paraViewAnimationFrames,
                                std::string animationFilename);

        /// Writes a pvtu file, which contains the links to all the rank files.
        void writeParallelFile(std::string name, int nRanks,
                               std::string fnPrefix, std::string fnPostfix,
                               const std::vector<std::string>& componentNames,
                               ::AMDiS::io::VtkWriter::Vtuformat format = ::AMDiS::io::VtkWriter::ASCII,
                               bool highPrecision = false,
                               bool writeAsVector = false,
                               bool createSubDir = false
                              );

        /// Writes a pvtu file which contains the links to the rank files in @subNames
        void writeParallelFile(std::string name, int nRanks,
                               std::vector<std::string>& subNames,
                               const std::vector<std::string>& componentNames,
                               ::AMDiS::io::VtkWriter::Vtuformat format = ::AMDiS::io::VtkWriter::ASCII,
                               bool highPrecision = false,
                               bool writeAsVector = false,
                               bool createSubDir = false
                              );

        /// Writes the connectivity for the case dim = 2 and degree = 2 to an output file.
        template<typename S, typename OutputStream>
        void writeConnectivity_dim2_degree2(DataCollector<S>* dataCollector, OutputStream& file);

        /// Writes the connectivity for the case dim = 2 and degree = 3 to an output file.
        template<typename S, typename OutputStream>
        void writeConnectivity_dim2_degree3(DataCollector<S>* dataCollector, OutputStream& file);

        /// Writes the connectivity for the case dim = 2 and degree = 4 to an output file.
        template<typename S, typename OutputStream>
        void writeConnectivity_dim2_degree4(DataCollector<S>* dataCollector, OutputStream& file);

        /// Writes a world coordinate to a given file.
        template<typename OutputStream>
        void writeCoord(OutputStream& file, WorldVector<double> coord)
        {
          file << std::scientific;
          file.precision(15);

          for (int i = 0; i < Global::getGeo(WORLD); i++)
            file << " " << coord[i];
          for (int i = Global::getGeo(WORLD); i < 3; i++)
            file << " "<<0.0;

          file << "\n";
        }

      } // end namespace detail
    } // end namespace VtkWriter
  }
} // end namespace io, AMDiS

#include "VtkWriter.hh"

#endif // AMDIS_VTKWRITER_DETAIL_H
