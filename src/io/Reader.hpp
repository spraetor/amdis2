#pragma once

#include <cstring>

#include "DOFVector.hpp"
#include "SystemVector.hpp"

// #include "io/ArhReader.hpp"
// #include "io/Arh2Reader.hpp"
#include "io/Arh3Reader.hpp"
#include "io/MacroReader.hpp"
#include "io/ValueReader.hpp"
#include "io/XYZReader.hpp"

#ifdef HAVE_PNG
#include "io/PngReader.hpp"
#endif

#ifdef HAVE_EXTENSIONS
#include "io/VtkReader.hpp"
#endif

#include "io/detail/ReaderWriter.hpp"


namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     * \brief reader interface that selects a reader depending on the file extension of the filename.
     *
     * @param filename Filename of the input file, used to extract extension.
     * @param container One of DOFVector&, std::vector<DOFVector*>,
     * 			SystemVector&, Mesh&.
     **/
    template<typename Container>
    void readFile(std::string filename,
                  Container& container)
    {
      std::string ext = filename.substr(filename.find_last_of("."));

      Mesh* mesh = detail::getMesh(container);
      if (ext == ".1d" || ext == ".2d" || ext == ".3d")
      {
        std::string periodicFilename = "";
        int check = 1;
        Parameters::get(mesh->getName() + "->periodic file", periodicFilename);
        Parameters::get(mesh->getName() + "->check", check);
        MacroReader::readMacro(filename, mesh, periodicFilename, check);
      }
      else if (ext == ".arh")
      {
//         if (Arh2Reader::isReadable(filename))
//           Arh2Reader::readFile(filename, container);
//         else if (Arh3Reader::isReadable(filename))
        TEST_EXIT(Arh3Reader::isReadable(filename))("File is not in ARH3 format. Please convert!\n");
        Arh3Reader::readFile(filename, container);
//         else
//           ArhReader::readFile(filename, container);
      }
      else if (ext == ".dat")
      {
        bool preserveMacroFileInfo = false;
        Parameters::get(mesh->getName() + "->preserve macroFileInfo",
                        preserveMacroFileInfo);
        TEST_EXIT(mesh->getMacroFileInfo() != NULL && preserveMacroFileInfo)
        ("Yout have to set the flag 'mesh_name->preserve macroFileInfo' "
         "to 'true' (%d), in order to read .dat-files",
         int(preserveMacroFileInfo));

        ValueReader::readValue(filename, mesh, container,
                               mesh->getMacroFileInfo());
      }
#ifdef HAVE_PNG
      else if (ext == ".png")
      {
        PngReader::readFile(filename, container);
      }
#endif

#ifdef HAVE_EXTENSIONS
      else if (ext == ".vtu")
      {
        VtkReader::readFile(filename, container);
      }
#endif
      else
      {
        ERROR_EXIT("File-extensions %s can not be assigned to a reader!\n", ext.c_str());
      }
    }


    /** \ingroup Output
     * \brief Wrapper for pointers to container types
     **/
    template<typename Container>
    void readFile(std::string filename,
                  Container* container)
    {
      readFile(filename, *container);
    }


    /** \ingroup Output
     * \brief specialization for read meshes
     **/
    inline
    void readFile(std::string filename,
                  Mesh* mesh)
    {

      std::string ext = filename.substr(filename.find_last_of("."));
      if (ext == ".1d" || ext == ".2d" || ext == ".3d")
      {
        std::string periodicFilename = "";
        int check = 1;
        Parameters::get(mesh->getName() + "->periodic file", periodicFilename);
        Parameters::get(mesh->getName() + "->check", check);
        MacroReader::readMacro(filename, mesh, periodicFilename, check);
      }
      else if (ext == ".arh")
      {
//         if (Arh2Reader::isReadable(filename))
//           Arh2Reader::readFile(filename, mesh);
        TEST_EXIT(Arh3Reader::isReadable(filename))("File is not in ARH3 format. Please convert!\n");
        Arh3Reader::readFile(filename, mesh);
//         else
//           ArhReader::readFile(filename, mesh);
      }
      else
      {
        ERROR_EXIT("File-extensions %s can not be assigned to a reader!\n", ext.c_str());
      }
    }
  } // end namespace io
} // end namespace AMDiS
