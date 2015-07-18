#ifndef AMDIS_ARH_READER2_DETAIL_H
#define AMDIS_ARH_READER2_DETAIL_H

#include "Mesh.h"
#include "MeshStructure.h"
#include "DOFVector.h"
#include "SystemVector.h"

namespace AMDiS { namespace io {

  namespace Arh2Reader 
  {
    namespace detail
    {
      
     /*
     * \ingroup Output
     * 
     * \brief
     * read the first part of Arh2 format file, and check the following:
     * 1. the type of file is equal to "arh2".
     * 2. the major version of Arh2Reader is equal to the one in the file.
     * 3. the minor version of Arh2Reader is bigger than the one in the file.
     * return value: minor version
     */
      uint8_t firstRead(std::ifstream& file);
      
      void setDofValues(int macroElIndex, Mesh *mesh,
			     std::vector<std::vector<double> >& values, 
                             std::vector<DOFVector<double>*>& vecs,
                             std::vector<std::vector<int> >& feSpaces);
      
       /*
	* \ingroup Output
	*
	* \brief 
	* when calling this function, one has to ensure that the following
	* conditions need to be fulfilled:
	* 
	* 0. DOFVectors in vecs are not allowed to have identical name.
	* 1. the length of vecs is less than the length of values in the file. 
	* 
	* In normal way(byName = false), it reads data and put it into the correspond 
	* vec(according to their order like old ArhReader). 
	* 
	* Condition 2 will be checked:
	* 2. all the non-null DOFvector in vecs have the same fespace
	* (number of Dofs per position) as the correspond value. 
	* 
	* If reading by name, it traverses all the values in the file to get the
	* first matching value with the same name and fespace. If doesn't find, error.
	* 
	* Condition 3 will be checked:
	* 3. There is at least one value in the file has the same name and fespace.
	* 
	* Note: Null DOFVector pointer in vecs is allowed.
      */
      void read(std::string filename,
                         Mesh *mesh,
                         std::vector<DOFVector<double>*> vecs,
                         bool byName = false);
      
      void readFile(std::string filename, 
		     Mesh *mesh,
		     std::vector<DOFVector<double>*> vecs,
		     bool writeParallel,
		     int nProcs = -1,
                     bool byName = false);
      
      /// read meta data from a single ARH-file 
      void readMetaFromSgArh(std::string filename, int nProc, 
			   std::vector< std::set<std::pair<int, int> > >& data);

      
      int readNumOfMacrosFromSgArh(std::string filename, int nProc = -1);

    }//end namespace detail
  } // end namespace Arh2Reader
} } // end namespace io, AMDiS

#endif