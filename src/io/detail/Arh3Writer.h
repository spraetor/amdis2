#ifndef AMDIS_ARH_WRITER3_DETAIL_H
#define AMDIS_ARH_WRITER3_DETAIL_H

#include "Global.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "DOFVector.h"
#include "SystemVector.h"
#include "boost/assign.hpp"

namespace AMDiS
{
  namespace io
  {

    namespace Arh3Writer
    {
      typedef enum
      {
        NONE = 0,
        ZLIB = 1,
        BZIP2 = 2
      } Cpsformat;

      typedef enum {SI08, SI16, SI32, SI64, UI08, UI16, UI32, UI64, SF32, SF64} Valformat;

      using namespace boost::assign;

      static const std::map<std::string, Valformat> dataformatMap = map_list_of("SI08", SI08)
          ("SI16", SI16)
          ("SI32", SI32)
          ("SI64", SI64)
          ("UI08", UI08)
          ("UI16", UI16)
          ("UI32", UI32)
          ("UI64", UI64)
          ("SF32", SF32)
          ("SF64", SF64);

      namespace detail
      {
        //Maybe remove this later
        void write(std::string filename,
                   DOFVector<double>* vec0 = NULL,
                   DOFVector<double>* vec1 = NULL,
                   DOFVector<double>* vec2 = NULL,
                   bool writeParallel = true,
                   Cpsformat cps = NONE,
                   std::string dataformat = "SF64",
                   std::string filenameType = "cont");

        /**
         * \ingroup Output
         *
         * \brief
         * checks the arguments, eventually splits the list of DOFVectors, to write
         * these to different files.
         *
         * The behavior is as follows:
         *
         * If mesh is given, then all the DOFVectors in vecs must belong to the
         * same mesh as the given one.
         *
         * Else vecs can contain DOFVector belong to different meshes and they
         * will be seperated to diff files. If multiple output files generated,
         * they are named like "filename.meshname.arh".
         *
         * Note: NULL pointer in vecs is not allowed for writeFile.
         *       If mesh is NULL and vecs is empty, a warning is showed and returned.
         *       Identical name in DOFVectors is not allowed.
        */
        void write(std::string filename,
                   Mesh* mesh,
                   std::vector<DOFVector<double>*> vecs,
                   bool writeParallel = true,
                   Cpsformat cps = NONE,
                   std::string dataformat = "SF64",
                   std::string filenameType = "cont");

        void writeAux(std::string filename, Mesh* mesh,
                      std::vector<DOFVector<double>*> vecs,
                      bool writeParallel, Cpsformat cps,
                      std::string dataformat);

        ///\return the size of the macro block in file
        std::pair<int, int> writeMacroElement(std::ofstream& file,
                                              MeshStructure& code,
                                              std::vector<std::vector<double>>& values,
                                              std::map<const FiniteElemSpace*,
                                              std::vector<int>>& feSpaces,
                                              Cpsformat cps, std::string dataformat);

        int writeValues(std::stringstream& file,
                        std::string dataformat,
                        std::vector<double>& values);

        template<typename T>
        int writeValues(std::stringstream& file, std::vector<double>& values);

        int writeHeader(std::ofstream& file,
                        Mesh* mesh,
                        std::vector<DOFVector<double>*> vecs,
                        std::map<const FiniteElemSpace*,
                        std::vector<int>>& feSpaces,
                        Cpsformat cps, std::string dataformat);

        ///internal method, don't call
        void setMacrosPos(std::ofstream& file, int headerLen,
                          std::vector<std::pair<int, int>>& macroSize);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
        void writeParallelFile(std::string filename,
                               Mesh* mesh,
                               std::string filenameType);
#endif
      }//end namespace detail
    } // end namespace Arh3Writer
  }
} // end namespace io, AMDiS

#endif
