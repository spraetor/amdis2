#pragma once

#ifdef HAVE_PNG

#include <cstring>

#include "DOFVector.hpp"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     *
     * \brief
     * Namespace which provieds readers to read a png file and gets values for each pixel
     */
    namespace PngReader
    {

      /// Interface for general containers not implemented. Specializations below.
      template<typename Container>
      void readFile(std::string filename, Container& vec)
      {
        ERROR_EXIT("PngReader not implemented for this container type!\n");
      }

      /// Copies the values of a png file to a DOF vector. Using container pointer.
      void readFile(std::string filename, DOFVector<double>* dofVector);


      /// Copies the values of a png file to a DOF vector. Using container reference.
      inline
      void readFile(std::string filename, DOFVector<double>& dofVector)
      {
        readFile(filename, &dofVector);
      }

    }

  }
} // end namespace io, AMDiS

#endif // HAVE_PNG
