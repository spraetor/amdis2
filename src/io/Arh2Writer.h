#ifndef AMDIS_ARH_WRITER2_H
#define AMDIS_ARH_WRITER2_H

#include "DOFVector.h"
#include "SystemVector.h"
#include "detail/Arh2Writer.h"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Output
      * \brief Writer for the AMDiS ARH-format - version 2
      *
      * A collection of methods to write various container types to
      * ARH-files.
      **/
    namespace Arh2Writer
    {
      const uint8_t MAJOR = 2;
      const uint8_t MINOR = 1;

      /**
        * \brief write the meshstructure and the dof values of DOFVectors in sysVec
        * to arh files.
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
        * Note:
        * 1. NULL pointer in vecs is not allowed for writeFile.
        * 2. Identical name in DOFVectors is not allowed.
        *
        * \param writeParallel
        * the default value is true.
        * In Sequential Model, the value is ignored...
        * In Parallel Model, you can set the value to false if you want to write
        * arh file without adding the processor number. But if you do so, it might
        * highly cause conflict when processors on the same machine try to access
        * the identical file.
      */
      inline void writeFile(SystemVector* sysVec,
                            std::string filename,
                            bool writeParallel = true)
      {
        std::vector<DOFVector<double>*> vecs;
        for (int i = 0; i < sysVec->getSize(); i++)
          vecs.push_back(sysVec->getDOFVector(i));
        detail::write(filename, NULL, vecs, writeParallel);
      }

      inline void writeFile(SystemVector& sysVec,
                            std::string filename,
                            bool writeParallel = true)
      {
        writeFile(&sysVec, filename, writeParallel);
      }

      /// write the meshstructure and the dof values of DOFVectors in vec0
      /// the behavior is equal to writeFile(SystemVector* sysVec, string filename).
      inline void writeFile(DOFVector<double>* vec0,
                            std::string filename,
                            bool writeParallel = true)
      {
        std::vector<DOFVector<double>*> vecs;
        vecs.push_back(vec0);
        detail::write(filename, NULL, vecs, writeParallel);
      }

      /// write the meshstructure and the dof values of DOFVectors in vec0
      /// the behavior is equal to writeFile(SystemVector* sysVec, string filename).
      inline void writeFile(DOFVector<double>& vec0,
                            std::string filename,
                            bool writeParallel = true)
      {
        writeFile(&vec0, filename, writeParallel);
      }

      /// write the meshstructure and the dof values of DOFVectors in vecs
      /// the behavior is equal to writeFile(SystemVector* sysVec, string filename).
      inline void writeFile(std::vector<DOFVector<double>*> vecs,
                            std::string filename,
                            bool writeParallel = true)
      {
        detail::write(filename, NULL, vecs, writeParallel);
      }

      /// write the meshstructure of the mesh to arh file.
      inline void writeFile(Mesh* mesh,
                            std::string filename,
                            bool writeParallel = true)
      {
        std::vector<DOFVector<double>*> vecs;
        detail::write(filename, mesh, vecs, writeParallel);
      }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
      void writeMetaData(Mesh* mesh, std::string filename);
#endif


    } // end namespace Arh2Writer
  }
} // end namespace io, AMDiS

#endif