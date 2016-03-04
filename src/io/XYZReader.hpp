#pragma once

#ifdef HAVE_EXTENSIONS

#include <cstring>

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     *
     * \brief
     * Namespace of methods which read a XYZ-file
     */
    namespace XYZReader
    {
      template<typename Container>
      void readFile(std::string filename, Container& data)
      {
        ERROR_EXIT("Can not read xyz-file to given container!\n");
      }

      /// Copies the values of a value file to a DOF vector.
      void readFile(std::string filename,
                    std::pair<std::vector<WorldVector<double>>,
                    std::vector<std::vector<double>>>& data);

    } // end namespace XYZReader
  }
} // end namespace io, AMDiS

#endif
