#ifndef AMDIS_ARH_READER_DETAIL_H
#define AMDIS_ARH_READER_DETAIL_H

#include "Mesh.h"
#include "DOFVector.h"


namespace AMDiS { namespace io {

  namespace ArhReader
  {
    namespace detail
    {
      
      void setDofValues(int macroElIndex, Mesh *mesh,
			     std::vector<double>& values, DOFVector<double>* vec);

      void readBlock(std::vector<char> &data,
			  Mesh *mesh,
			  std::vector<DOFVector<double>*> vecs);
      
      void read(std::string filename,
			 Mesh *mesh,
			 std::vector<DOFVector<double>*> vecs);
      
      void readFile(std::string filename, 
		     Mesh *mesh,
		     std::vector<DOFVector<double>*> vecs,
		     bool writeParallel,
		     int nProcs = -1);
      
    }//end namespace detail
  } // end namespace ArhReader
} } // end namespace io, AMDiS

#endif
