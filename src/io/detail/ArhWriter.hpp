#ifndef AMDIS_ARH_WRITER_DETAIL_H
#define AMDIS_ARH_WRITER_DETAIL_H

#include "Global.h"
#include "DOFVector.h"
#include "SystemVector.h"
#include "MeshStructure.h"

namespace AMDiS
{
  namespace io
  {

    namespace ArhWriter
    {
      namespace detail
      {
        void write(std::string filename, Mesh* mesh,
                   DOFVector<double>* vec0 = NULL,
                   DOFVector<double>* vec1 = NULL,
                   DOFVector<double>* vec2 = NULL);

        void write(std::string filename, Mesh* mesh,
                   std::vector<DOFVector<double>*> vecs,
                   bool writeParallel = true);

        void writeMacroElement(std::ofstream& file,
                               MeshStructure& code,
                               std::vector<std::vector<double>>& values,
                               uint32_t elIndex);

      }//end namespace detail
    } // end namespace ArhWriter
  }
} // end namespace io, AMDiS

#endif
