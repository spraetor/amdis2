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


#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/convenience.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/copy.hpp>
#ifdef HAVE_COMPRESSION
#include <boost/iostreams/filter/zlib.hpp>
#endif
#include "AdaptInfo.h"
#include "VtkWriter.h"

namespace AMDiS
{
  namespace io
  {

    namespace VtkWriter
    {

      using namespace std;

      int Aux::writeFile(string name)
      {
        FUNCNAME("VtkWriter::writeFile()");

        boost::iostreams::filtering_ostream file;

        {
          ofstream swapfile(name.c_str(), ios::out | ios::trunc);
          TEST_EXIT(swapfile.is_open())
          ("Cannot open file %s for writing!\n", name.c_str());
          swapfile.close();
        }

        file.push(boost::iostreams::file_descriptor_sink(name, ios::trunc));

        if(format == APPENDED || format == APPENDED_COMPRESSED)
          writeFileToStreamAppended(file);
        else
        {
          TEST_EXIT(format == ASCII)("Unknown ParaView mode.\n");
          writeFileToStream(file);
        }
        return 0;
      }

      std::string Aux::getStreamData()
      {
        using namespace std;

        string result(""), header(""), body("");

        switch(format)
        {
        case APPENDED:
        {
          BinaryStream hstream;
          hstream << bstream.getSize();
          body = bstream.str();
          header = hstream.str();
          result = detail::base64Encode(header + body);
          break;
        }
        case APPENDED_COMPRESSED:
        {
#ifdef HAVE_COMPRESSION
          BinaryStream hstream;
          int nBlocks = bstream.getSize() / AMDIS_ZLIB_BLOCK_SIZE + 1;
          int finalsize = bstream.getSize() % AMDIS_ZLIB_BLOCK_SIZE;
          hstream << nBlocks << AMDIS_ZLIB_BLOCK_SIZE << finalsize;

          string data = bstream.str();

          for(int i = 0; i < nBlocks; i++)
          {
            string subData = (nBlocks - 1 == i) ? data.substr(i * AMDIS_ZLIB_BLOCK_SIZE, finalsize) : data.substr(i * AMDIS_ZLIB_BLOCK_SIZE, AMDIS_ZLIB_BLOCK_SIZE);
            stringstream tmp1, tmp2;
            tmp1 << subData;
            boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
            in.push(boost::iostreams::zlib_compressor());
            in.push(tmp1);

            boost::iostreams::copy(in, tmp2);
            subData = tmp2.str();
            int size = subData.length();
            hstream << size;
            body += subData;
          }
          header = hstream.str();
          result = detail::base64Encode(header) + detail::base64Encode(body);
#endif
          break;
        }
        case ASCII:
        default:
          break;
        }
        bstream.str("");
        return result;
      }

      // __________________________________________________________________________ //

      namespace detail
      {
        std::string base64Encode(std::string text)
        {
          unsigned int writePaddChars = (3-text.length()%3)%3;
          std::string base64(binary_base64(text.begin()),binary_base64(text.end()));
          base64.append(writePaddChars,'=');
          return base64;
        }

        std::string extract_relative_path(std::string valueFilename, std::string animationFilename)
        {
          using namespace boost::filesystem;
          path vtu_path = valueFilename;
          path pvd_path = animationFilename;

          path::iterator it_vtu, it_pvd;
#if (BOOST_VERSION < 104800)
          path vtu_path0 = complete(vtu_path);
          vtu_path0.remove_filename();
          path pvd_path0 = complete(pvd_path);
          pvd_path0.remove_filename();
#else
          path vtu_path0 = absolute(vtu_path);
          vtu_path0.remove_filename();
          path pvd_path0 = absolute(pvd_path);
          pvd_path0.remove_filename();
#endif
          // find matching root directories
          for (it_vtu = vtu_path0.begin(), it_pvd = pvd_path0.begin();
               it_vtu != vtu_path0.end() && it_pvd != pvd_path0.end() && *it_vtu == *it_pvd;
               it_vtu++, it_pvd++) {}

          // create relative path
          path new_vtu_path;
          for (; it_pvd != pvd_path0.end(); it_pvd++)
            new_vtu_path /= "..";
          for (; it_vtu != vtu_path0.end(); it_vtu++)
            new_vtu_path /= *it_vtu;
          new_vtu_path /= vtu_path.filename();

          return new_vtu_path.string();
        }

        int updateAnimationFile(AdaptInfo& adaptInfo,
                                string valueFilename,
                                vector<pair<double, string>>* paraViewAnimationFrames,
                                string animationFilename)
        {
          FUNCNAME("updateAnimationFile()");

          paraViewAnimationFrames->push_back(
            make_pair(adaptInfo.getTime(), extract_relative_path(valueFilename, animationFilename)));

          boost::iostreams::filtering_ostream file;
          {
            ofstream swapfile(animationFilename.c_str(),
                              ios::out | ios::trunc);
            TEST_EXIT(swapfile.is_open())
            ("Cannot open file %s for writing!\n", animationFilename.c_str());
            swapfile.close();
          }
          file.push(boost::iostreams::file_descriptor_sink(animationFilename,
                    ios::trunc));

          file << "<?xml version=\"1.0\"?>\n";
          file << "<VTKFile type=\"Collection\" version=\"0.1\" >"  << "\n";
          file << "<Collection>\n";

          for (vector<pair<double, string>>::iterator it = paraViewAnimationFrames->begin();
               it < paraViewAnimationFrames->end(); ++it)
          {
            file << "<DataSet timestep=\"" << it->first
                 << "\" part=\"0\" file=\"" << it->second << "\"/>\n";
          }

          file << "</Collection>\n";
          file << "</VTKFile>\n";

          return 0;
        }




        void writeParallelFile(string name, int nRanks,
                               string fnPrefix, string fnPostfix,
                               const vector<string>& componentNames,
                               ::AMDiS::io::VtkWriter::Vtuformat format,
                               bool highPrecision,
                               bool writeAsVector,
                               bool createSubDir
                              )
        {
          vector<string> fileNames(nRanks);
          for (int i = 0; i < nRanks; i++)
            fileNames[i] = fnPrefix + "-p" + std::to_string(i) + "-" + fnPostfix;

          writeParallelFile(name, nRanks, fileNames, componentNames, format, highPrecision, writeAsVector, createSubDir);
        }

        void writeParallelFile(string name, int nRanks,
                               vector<string>& subNames,
                               const vector<string>& componentNames,
                               ::AMDiS::io::VtkWriter::Vtuformat format,
                               bool highPrecision,
                               bool writeAsVector,
                               bool createSubDir
                              )
        {
          FUNCNAME("writeParallelFile()");
#ifndef HAVE_COMPRESSION
          if(format == APPENDED_COMPRESSED)
            format = APPENDED;
#endif

          boost::iostreams::filtering_ostream file;
          {
            ofstream swapfile(name.c_str(), ios::out | ios::trunc);
            TEST_EXIT(swapfile.is_open())
            ("Cannot open file %s for writing!\n", name.c_str());
            swapfile.close();
          }
          file.push(boost::iostreams::file_descriptor_sink(name, ios::trunc));

          file << "<?xml version=\"1.0\"?>\n";

          if(format == ::AMDiS::io::VtkWriter::APPENDED_COMPRESSED)
            file << "<VTKFile type=\"PUnstructuredGrid\" compressor=\"vtkZLibDataCompressor\" >\n";
          else
            file << "<VTKFile type=\"PUnstructuredGrid\" >\n";

          file << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
          file << "    <PPoints>\n";

          if(highPrecision && format != ::AMDiS::io::VtkWriter::ASCII)
            file << "      <PDataArray type=\"Float64\"";
          else
            file << "      <PDataArray type=\"Float32\"";

          if(format == ::AMDiS::io::VtkWriter::ASCII)
            file << " NumberOfComponents=\"3\" format=\"ascii\"/>\n";
          else
            file << " NumberOfComponents=\"3\" format=\"appended\"/>\n";

          file << "    </PPoints>\n"
               << "    <PCells>\n"
               << "      <PDataArray type=\"Int32\" Name=\"offsets\"/>\n"
               << "      <PDataArray type=\"UInt8\" Name=\"types\"/>\n"
               << "      <PDataArray type=\"Int32\" Name=\"connectivity\"/>\n"
               << "    </PCells>\n";
          file << "    <PPointData>\n";

          int nValues = static_cast<int>(componentNames.size());
          nValues = writeAsVector ? std::min(nValues,1) : nValues;
          for (int i = 0; i < nValues; i++)
          {
            if(highPrecision && format != ::AMDiS::io::VtkWriter::ASCII)
              file << "      <PDataArray type=\"Float64\" Name=\"";
            else
              file << "      <PDataArray type=\"Float32\" Name=\"";

            file << componentNames[i];

            if(format == ::AMDiS::io::VtkWriter::ASCII)
              file << "\" format=\"ascii\"";
            else
              file << "\" format=\"appended\"";
            if (writeAsVector && componentNames.size() > 1)
              file << " NumberOfComponents=\"" << std::max(3,static_cast<int>(componentNames.size())) << "\"";
            file << "/>\n";
          }


          file << "    </PPointData>\n";

          for (int i = 0; i < nRanks; i++)
          {
            boost::filesystem::path filepath(subNames[i]);
            file << "    <Piece Source=\"";
            if (createSubDir)
              file << "./data/";
            file << boost::filesystem::basename(filepath)
                 << boost::filesystem::extension(filepath) << "\"/>\n";
          }

          file << "  </PUnstructuredGrid>\n";
          file << "</VTKFile>\n";
        }

      } // end namespace detail
    } // end namespace VtkWriter
  }
} // end namespace io, AMDiS
