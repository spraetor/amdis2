#ifndef AMDIS_ARH_WRITER2_DETAIL_H
#define AMDIS_ARH_WRITER2_DETAIL_H

#include "Global.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "DOFVector.h"
#include "SystemVector.h"

namespace AMDiS { namespace io {

  namespace Arh2Writer
  {
    namespace detail
    {
      //Maybe remove this later
      void write(std::string filename,                                                
                      DOFVector<double>* vec0 = NULL,                                                  
                      DOFVector<double>* vec1 = NULL,
                      DOFVector<double>* vec2 = NULL,
		      bool writeParallel = true);
    
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
		  bool writeParallel = true);

       void writeAux(std::string filename, Mesh *mesh,
		       std::vector<DOFVector<double>*> vecs,
		       bool writeParallel);

       ///\return the size of the macro block in file
       int writeMacroElement(std::ofstream &file,
                                  MeshStructure &code,
                                  std::vector<std::vector<double> >& values,
                                  std::map<const FiniteElemSpace*, std::vector<int> >& feSpaces);

       int writeHeader(std::ofstream& file,
                            Mesh* mesh,
                            std::vector<DOFVector<double>*> vecs,
                            std::map<const FiniteElemSpace*, std::vector<int> >& feSpaces);    
       
       ///internal method, don't call
       void setMacrosPos(std::ofstream& file,
	                    std::vector<int>& macrosPosVec);
 
    }//end namespace detail
  } // end namespace Arh2Writer
} } // end namespace io, AMDiS

#endif