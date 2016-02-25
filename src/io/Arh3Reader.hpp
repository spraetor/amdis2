#pragma once

#include "AMDiS_fwd.hpp"
#include "Global.hpp"

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
      * \brief Reader for the AMDiS ARH-format - version 3
      *
      * A collection of methods to read to various container types from
      * ARH-files.
      **/
    namespace Arh3Reader
    {
      const uint8_t MAJOR = 3;
      const uint8_t MINOR = 0;

      /**
       * \brief Read MeshStructure, refine the mesh and read dof values to sysVec by order.
       * When reading dof values, it puts the data into the correspond DOFVector
       * (according to their order in the sysVec and in the file, which is the same as the old ArhReader).
       * Normally, just use the default value of \ref writeParallel and \ref nProcs.
       *
       * You should notice these:
       * 1. DOFVectors in SystemVector are not allowed to have identical name.
       * 2. The length of DOFVectors in SystemVector is less than the length of values in the file.
       * 3. NULL of DOFVector is allowed in SystemVector. But all the non-null DOFvector have the same fespace
       * (number of Dofs per position) as the correspond value in the file.
       *
       * \param writeParallel
       * the default value in Sequential Model is false.
       * the default value in Parallel Model is true.
       * if the value is set to false, no matter in sequential or parallel mode, the processor(s) will
       *   read the file of the same given name, without adding thier process numbers: for example: "-p1-".
       * if the value is set to true, in Parallel Model, the behavior depends on \ref nProcs.
       *   But in Sequential Model, an error will be thrown.
       * \param nProcs
       * the default value is -1. But it only affects when \ref writeParallel is set to true.
       * If nProcs is -1, and it's in Parallel Model, every processor reads their own filename (with processor
       *   numbers). But in Sequential Model, an error will be thrown.
       * And you can also set nProcs to a value which is smaller than the number of Processors. Then every
       *   processor will read all the files with names from -p0- to -p[nProcs]-.
       */
      void readFile(std::string filename,
                    SystemVector* sysVec,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1);

      inline
      void readFile(std::string filename,
                    SystemVector& sysVec,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1)
      {
        readFile(filename, &sysVec, writeParallel, nProcs);
      }


      /// Read MeshStructure, refine the mesh and read dof values to vec0, vec1 and vec2 by order.
      /// the behavior is equal to readFile(string filename, SystemVector* sysVec).
      void readFile(std::string filename,
                    DOFVector<double>* vec0 = NULL,
                    DOFVector<double>* vec1 = NULL,
                    DOFVector<double>* vec2 = NULL,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1);

      inline
      void readFile(std::string filename,
                    DOFVector<double>& vec0,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1)
      {
        readFile(filename, &vec0, NULL, NULL, writeParallel, nProcs);
      }

      /// Read MeshStructure, refine the mesh and read dof values to vecs by order.
      /// the behavior is equal to readFile(string filename, SystemVector* sysVec).
      void readFile(std::string filename,
                    std::vector<DOFVector<double>*> vecs,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1);

      /// Read MeshStructure and refine the mesh from arh file.
      void readFile(std::string filename,
                    Mesh* mesh,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1);

      inline
      void readFile(std::string filename,
                    Mesh& mesh,
                    bool writeParallel = WRITE_PARALLEL,
                    int nProcs = -1)
      {
        readFile(filename, &mesh, writeParallel, nProcs);
      }

      /**
       * \brief Read MeshStructure, refine the mesh and read dof values to sysVec by name.
       * It traverses all the values in the file to get the first matching value with the same name
       * and fespace, and then reads the data. If it doesn't find any matching, an error will be thrown.
       * Normally, just use the default value of \ref writeParallel and \ref nProcs. See readFile for
       * detail information.
       *
       * You should notice these:
       * 1. DOFVectors in SystemVector are not allowed to have identical name.
       * 2. The length of SystemVector is less than the length of values in the file.
       * 3. There is at least one value in the file has the same name and fespace as the DOFVector in SystemVector.
       */
      void readByName(std::string filename,
                      SystemVector* sysVec,
                      bool writeParallel = WRITE_PARALLEL,
                      int nProcs = -1);

      /// Read MeshStructure, refine the mesh and read dof values to vecs by name.
      /// the behavior is equal to readByName(string filename, SystemVector* sysVec).
      void readByName(std::string filename,
                      std::vector<DOFVector<double>*> vecs,
                      bool writeParallel = WRITE_PARALLEL,
                      int nProcs = -1);

      /// Read MeshStructure, refine the mesh and read dof values to vec by name.
      /// the behavior is equal to readByName(string filename, SystemVector* sysVec).
      void readByName(std::string filename,
                      DOFVector<double>& vec,
                      bool writeParallel = WRITE_PARALLEL,
                      int nProcs = -1);

      /// Read MeshStructure, refine the mesh and read dof values to vec by name.
      /// the behavior is equal to readByName(string filename, SystemVector* sysVec).
      void readByName(std::string filename,
                      DOFVector<double>* vec,
                      bool writeParallel = WRITE_PARALLEL,
                      int nProcs = -1);


      /// the behavior is equal to readMeta(string filename, vector(DOFvector*) vecs).
      void readMeta(std::string filename,
                    DOFVector<double>* vec0 = NULL,
                    DOFVector<double>* vec1 = NULL,
                    DOFVector<double>* vec2 = NULL);

      /**
       * \brief first read a meta ARH file and get \ref elInRank. And then uses the elInRank map
       * to find all the arh files that contains the dof value vecs needs and sets data into vecs by order.
       * \param filename the name of meta ARH file.
       * \param vecs the vector of DOFVectors which you want to get value.
       */
      void readMeta(std::string filename,
                    std::vector<DOFVector<double>*> vecs);

      /**
       * \brief Read meta data from a meta ARH file and put the information into \ref elInRank
       * and \ref elCodeSize and \ref arhPrefix.
       * \param elInRank map of macro index and the rank it belongs to.
       * \param elCodeSize map of macro index and the code size of MeshStructure of the macro.
       * \param arhPrefix the prefix of arh file which the meta data file comes from.
       * \return the number of subdomains a meta ARH file is defined for.
       */
      int readMetaData(std::string filename,
                       std::map<int, int>& elInRank,
                       std::map<int, int>& elCodeSize,
                       std::string& arhPrefix);

      ///Read meta data from a meta ARH file and put the information into \ref elInRank
      ///and \ref elCodeSize. And return the number of subdomains a meta ARH file is defined for.
      inline int readMetaData(std::string filename,
                              std::map<int, int>& elInRank,
                              std::map<int, int>& elCodeSize)
      {
        std::string tmp;
        return readMetaData(filename, elInRank, elCodeSize, tmp);
      }

      /// Same as readMetaData  but collects inform^ation from a set of ARH-files
      int readMetaFromArh(std::string filename,
                          std::map<int, int>& elInRank,
                          std::map<int, int>& elCodeSize);

      /// Only returns just the number of subdomains a meta ARH file is defined for.
      int readMetaData(std::string filename);


      /// Returns the number of value vectors in the file.
      int readNumOfValueVectors(std::string filename, bool writeParallel = WRITE_PARALLEL);

      /// Returns the Header size of the file.
      int readHeaderSize(std::string filename, bool writeParallel = WRITE_PARALLEL);

      /// If the current version Arh3Reader can read the file, return true, else false.
      bool isReadable(std::string filename, bool writeParallel = WRITE_PARALLEL);

    } // end namespace Arh3Reader
  }
} // end namespace io, AMDiS
