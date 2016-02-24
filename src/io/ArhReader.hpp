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



/** \file ArhReader.h */

#ifndef AMDIS_ARH_READER_H
#define AMDIS_ARH_READER_H

#include "Mesh.h"
#include "DOFVector.h"
#include "Global.h"
#include "SystemVector.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#define WRITE_PARALLEL true
#else
#define WRITE_PARALLEL false
#endif

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
      * \brief Reader for the AMDiS ARH-format - version 1
      *
      * A collection of methods to read to various container types from
      * ARH-files.
      **/
    namespace ArhReader
    {

      /// Read a file to 1-3 \ref DOFVector. Using container pointers.
      void readFile(std::string filename,
                    DOFVector<double>* vec0 = NULL,
                    DOFVector<double>* vec1 = NULL,
                    DOFVector<double>* vec2 = NULL,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1);


      /// Read a file to a \ref DOFVector. Using container reference.
      inline
      void readFile(std::string filename,
                    DOFVector<double>& vec0,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1)
      {
        readFile(filename, &vec0, NULL, NULL, writeParallel, nProcs);
      }



      /// Read a file to a vector of \ref DOFVector.
      void readFile(std::string filename,
                    std::vector<DOFVector<double>*> vecs,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1);


      /// Read a file to a \ref SystemVector. Using container pointer.
      inline
      void readFile(std::string filename,
                    SystemVector* sysVec,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1)
      {
        std::vector<DOFVector<double>*> vecs;
        for (int i = 0; i < sysVec->getSize(); i++)
          vecs.push_back(sysVec->getDOFVector(i));
        readFile(filename, vecs, writeParallel, nProcs);
      }


      /// Read a file to a \ref SystemVector. Using container reference.
      inline
      void readFile(std::string filename,
                    SystemVector& sysVec,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1)
      {
        readFile(filename, &sysVec, writeParallel, nProcs);
      }


      /// Read a file to a \ref Mesh.
      void readFile(std::string filename,
                    Mesh* mesh,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1);


      /** \brief first read a meta ARH file and get \ref elInRank. And then uses the elInRank map
        * to find all the arh files that contains the dof value vecs needs and sets data into vecs by order.
        * \param filename the name of meta ARH file.
        * \param vec0 the first \ref DOFVector you want to read to.
        * \param vec1 the second \ref DOFVector you want to read to.
        * \param vec2 the third \ref DOFVector you want to read to.
        */
      void readMeta(std::string filename,
                    DOFVector<double>* vec0 = NULL,
                    DOFVector<double>* vec1 = NULL,
                    DOFVector<double>* vec2 = NULL);


      /// read a meta ARH file and the corresponding parallel ARH files to a vector
      /// of \ref DOFVector.
      void readMeta(std::string filename, std::vector<DOFVector<double>*> vecs);


      /// read only the meta data from a meta ARH file
      int readMetaData(std::string filename,
                       std::map<int, int>& elInRank,
                       std::map<int, int>& elCodeSize,
                       std::string& arhPrefix);


      /// read only the meta data from a meta ARH file
      inline int readMetaData(std::string filename,
                              std::map<int, int>& elInRank,
                              std::map<int, int>& elCodeSize)
      {
        std::string tmp;
        return readMetaData(filename, elInRank, elCodeSize, tmp);
      }


      /// Returns just the number of subdomains a meta ARH file is defined for.
      int readMetaData(std::string filename);


      /// Reader like readFile, but reads to data from the block \ref data.
      void readFromMemoryBlock(std::vector<char>& data, Mesh* mesh,
                               DOFVector<double>* vec0 = NULL,
                               DOFVector<double>* vec1 = NULL,
                               DOFVector<double>* vec2 = NULL,
                               bool writeParallel = WRITE_PARALLEL,
                               int nProcs = -1);


      /// Reader like readFile, but reads to data from the block \ref data.
      void readFromMemoryBlock(std::vector<char>& data, Mesh* mesh,
                               std::vector<DOFVector<double>*> vecs,
                               bool writeParallel = WRITE_PARALLEL,
                               int nProcs = -1);


      /// read the header information of an ARH file to extract the number of components
      int readNumOfValueVectors(std::string filename);


      // ________ below are obsolete functions, for backward compatibility _______


      void readFile(std::string filename,
                    Mesh* mesh,
                    std::vector<DOFVector<double>*> vecs);

      void read(std::string filename,
                Mesh* mesh,
                DOFVector<double>* vec0 = NULL,
                DOFVector<double>* vec1 = NULL,
                DOFVector<double>* vec2 = NULL,
                bool writeParallel = WRITE_PARALLEL,
                int nProcs = -1);

      void read(std::string filename,
                Mesh* mesh,
                std::vector<DOFVector<double>*> vecs,
                bool writeParallel = WRITE_PARALLEL,
                int nProcs = -1);

      int getNumValueVectors(std::string filename);

      void setDofValues(int macroElIndex, Mesh* mesh,
                        std::vector<double>& values, DOFVector<double>* vec);

      void readBlock(std::vector<char>& data,
                     Mesh* mesh,
                     std::vector<DOFVector<double>*> vecs);

      void readMeta(std::string filename,
                    Mesh* mesh,
                    std::vector<DOFVector<double>*> vecs);

      void readMeta(std::string filename,
                    Mesh* mesh,
                    DOFVector<double>* vec0 = NULL,
                    DOFVector<double>* vec1 = NULL,
                    DOFVector<double>* vec2 = NULL);

    } // end namespace ArhReader
  }
} // end namespace io, AMDiS

#endif
